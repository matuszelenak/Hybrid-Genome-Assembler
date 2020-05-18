#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string/join.hpp>
#include <fmt/format.h>
#include <filesystem>
#include <map>
#include <queue>
#include <set>

#include "../common/SequenceRecordIterator.h"
#include "../common/Utils.h"

#include "../lib/BloomFilter.h"


#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;

typedef std::pair<uint32_t, bool> Endpoint;
typedef std::pair<uint32_t, uint32_t> Interval;


struct ReadMetaData {
    ReadID id;
    uint32_t length;
    CategoryID category_id;
    uint32_t start;
    uint32_t end;
    tsl::robin_map<KmerID, uint32_t> kmer_positions;
};


class ReadComponent {
public:
    ReadComponent(ReadMetaData &initial_read_meta, std::vector<KmerID> &discriminative_kmer_ids){
        this->id = initial_read_meta.id;
        this->discriminative_kmer_ids = discriminative_kmer_ids;
        this->contained_read_ids = {initial_read_meta.id};
        this->categories = {initial_read_meta.category_id};
    };

    ComponentID id = 0;
    std::vector<ReadID> contained_read_ids;
    std::vector<KmerID> discriminative_kmer_ids;
    std::set<CategoryID> categories;

    [[nodiscard]] uint64_t size() const{ return contained_read_ids.size(); };

    std::string consistency(tsl::robin_map<ReadID, ReadMetaData> &metas){
        std::map<CategoryID, int> category_counts = {{0, 0}, {1, 0}};
        for (auto read_id : contained_read_ids) {
            category_counts.insert({metas[read_id].category_id, 0}).first->second += 1;
        }
        std::vector<std::string> category_count_vector;
        std::transform(
                category_counts.begin(),
                category_counts.end(),
                std::back_inserter(category_count_vector),
                [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
        );
        return boost::algorithm::join(category_count_vector, "/");
    }

    std::map<CategoryID, std::vector<Interval>> category_intervals(tsl::robin_map<ReadID, ReadMetaData> &metas){
        std::map<CategoryID, std::vector<Endpoint>> category_endpoints;
        for (auto read_id: contained_read_ids){
            category_endpoints.insert({metas[read_id].category_id, {}}).first->second.emplace_back(metas[read_id].start, true);
            category_endpoints.insert({metas[read_id].category_id, {}}).first->second.emplace_back(metas[read_id].end, false);
        }

        std::map<CategoryID, std::vector<Interval>> result;
        for (auto endpoints : category_endpoints){
            std::sort(endpoints.second.begin(), endpoints.second.end());

            result.insert({endpoints.first, {}});
            uint32_t current_interval_start;
            int opened_intervals = 0;
            for (auto endpoint : endpoints.second) {
                if (endpoint.second) {
                    opened_intervals++;
                    if (opened_intervals == 1) {
                        current_interval_start = endpoint.first;
                    }
                } else {
                    opened_intervals--;
                    if (opened_intervals == 0) {
                        result[endpoints.first].push_back({current_interval_start, endpoint.first});
                    }
                }
            }
        }
        return result;
    }

    std::string to_string(tsl::robin_map<ReadID, ReadMetaData> &metas){
        std::vector<std::string> category_intervals_strings;
        for (auto intervals : category_intervals(metas)){
            std::vector<std::string> interval_strings;
            std::transform(
                    intervals.second.begin(),
                    intervals.second.end(),
                    std::back_inserter(interval_strings),
                    [](Interval &i) -> std::string { return fmt::format("({},{})", i.first, i.second); }
            );
            category_intervals_strings.push_back(boost::algorithm::join(interval_strings, ";"));
        }

        return fmt::format("#{} : {} [{}]", this->id, this->consistency(metas), boost::algorithm::join(category_intervals_strings, " / "));
    }
};

typedef tsl::robin_map<ComponentID, ReadComponent *> ComponentIndex;
typedef std::unordered_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;

typedef std::vector<std::vector<ComponentID>> KmerComponentIndex;
typedef tsl::robin_map<KmerID, std::vector<ComponentID>> IndexRemovalMap;
typedef std::vector<ComponentID> Component;

typedef std::vector<std::pair<ComponentID, ComponentID> > SpanningTree;


struct ComponentConnection {
    ComponentID component_x_id;
    ComponentID component_y_id;
    ConnectionScore score;
    bool is_good; //DEBUG

    bool operator<(const ComponentConnection &conn) const {
        return (this->score < conn.score);
    }
};

struct ReadClusteringConfig {
    int scaffold_component_min_size = 30;
    int scaffold_component_max_size = -1;
    double scaffold_forming_fraction = 0.15;
    ConnectionScore scaffold_forming_score = 0;
    ConnectionScore enrichment_connections_min_score = 20;
    ConnectionScore tail_amplification_min_score = 40;
    int threads = 1;
    int spectral_dims = 16;
    bool force_spectral = false;
};

class ReadClusteringEngine {
protected:
    ReadClusteringConfig config;
    SequenceRecordIterator *reader;

    ComponentIndex component_index;
    KmerComponentIndex kmer_component_index;

    tsl::robin_map<ReadID, ReadMetaData> read_metas;

    // DEBUG
    bool debug = false;

    std::vector<Interval> get_read_coverage_intervals(std::vector<ReadID> &read_ids);

    void export_spanning_tree(SpanningTree &tree, std::pair<std::vector<ReadID>, std::vector<ReadID>> &tails);

    void plot_read_overlaps(std::vector<ComponentConnection> &connections);

    void print_components(std::vector<ComponentID> &ids);

    void plot_read_coverage(std::vector<ComponentID> &components);
    // END DEBUG

    int construct_indices(std::unordered_set<Kmer> &discriminative_kmers, int k);

    std::vector<ComponentConnection> get_all_connections(ConnectionScore min_score);

    std::vector<ComponentConnection> get_connections(Component &component, ConnectionScore min_score);

    std::vector<ComponentID> merge_components(std::vector<Component> &components);

    int approximate_read_overlap(ReadID x, ReadID y);

    std::pair<std::vector<ReadID>, std::vector<ReadID>> get_spanning_tree_tails(SpanningTree &tree);

    std::vector<ComponentConnection> get_core_component_connections(std::vector<std::pair<Component, SpanningTree>> &components_with_trees);

    std::vector<KmerID> accumulate_kmer_ids(std::vector<ComponentID> &component);

    std::vector<ComponentID> amplify_component(Component &component, int min_score);

public:
    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator, ReadClusteringConfig config);

    std::vector<ComponentID> run_clustering(std::unordered_set<Kmer> &discriminative_kmers, int k);

    void export_components(std::vector<ComponentID> &component_ids, std::string &directory_path);
};


#endif //SRC_READCLUSTERINGENGINE_H
