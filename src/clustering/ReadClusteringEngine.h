#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <experimental/filesystem>
#include <map>
#include <queue>

#include "../common/SequenceRecordIterator.h"
#include "../common/Utils.h"

#include "../lib/BloomFilter.h"

#include "ReadComponent.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ComponentID, ReadComponent*> ComponentIndex;

typedef std::vector<std::vector<ComponentID>> KmerComponentIndex;
typedef tsl::robin_map<KmerID, std::vector<ComponentID>> IndexRemovalMap;
typedef std::vector<ComponentID > Component;

typedef std::vector<std::pair<ComponentID, ComponentID> > SpanningTree;


struct ComponentConnection{
    ComponentID component_x_id;
    ComponentID component_y_id;
    ConnectionScore score;
    bool is_good;

    bool operator < (const ComponentConnection& conn) const
    {
        return (this->score < conn.score);
    }
};

void plot_connection_quality(std::vector<ComponentConnection> &connections);

class ReadClusteringEngine {
protected:
    SequenceRecordIterator* reader;
    bloom::BloomFilter<Kmer>* kmers = nullptr;
    int k = 0;

    ComponentIndex component_index;
    KmerComponentIndex kmer_component_index;

    tsl::robin_map<ReadID, uint32_t> read_lengths;
    tsl::robin_map<ReadID, tsl::robin_map<KmerID, int>> kmer_positions;

    // DEBUG
    tsl::robin_map<ReadID, std::pair<uint32_t, uint32_t> > read_intervals;

    int construct_indices();

    std::vector<ComponentConnection> get_connections(Component &component, ConnectionScore min_score);

    std::vector<ComponentID> merge_components(std::vector<Component> &components);

    std::vector<std::pair<Component, SpanningTree>> union_find(std::vector<ComponentConnection> &connections, std::set<ComponentID> &restricted, int min_component_size);

    std::pair<std::vector<ReadID>, std::vector<ReadID>> get_spanning_tree_boundary_reads(SpanningTree &tree);
    std::vector<ComponentConnection> get_core_component_connections(std::vector<std::pair<Component, SpanningTree>> &components_with_trees);

    std::map<ComponentID, std::vector<std::string>> assemble_clusters(std::vector<ComponentID> &cluster_ids);

    std::vector<Interval> get_read_coverage_intervals(std::vector<ReadID> &read_ids);

    std::vector<ComponentConnection> filter_connections(std::vector<ComponentConnection> &original, const std::function<bool(ComponentConnection &)> &func) {
        std::vector<ComponentConnection> result;
        for (auto conn : original) {
            if (func(conn)) result.push_back(conn);
        }
        return result;
    }

    void remove_merged_components(){
        for(auto it = component_index.begin(); it != component_index.end();){
            if (it->second->size() == 0){
                free(it->second);
                component_index.erase(it++);
            } else ++it;
        }
    }

public:
    void export_spanning_tree(SpanningTree &tree, std::pair<std::vector<ReadID>, std::vector<ReadID>> &tails);
    void plot_read_overlaps(std::vector<ComponentConnection> &connections);
    void print_components(std::vector<ComponentID> &ids);

    std::map<ComponentID, std::string> export_components(std::vector<ComponentID> &component_ids, std::experimental::filesystem::path &directory_path);

    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator);
    ~ReadClusteringEngine();

    void set_kmers(tsl::robin_set<Kmer> &_kmers, int k);
    void set_kmers(bloom::BloomFilter<Kmer>* _kmers, int k);
    std::vector<ComponentID> run_clustering();

    int approximate_read_overlap(ReadID x, ReadID y);

    void plot_overlap_differences(std::vector<ComponentConnection> &connections);

    std::vector<KmerID> accumulate_kmer_ids(std::vector<ComponentID> &component);

    std::vector<ComponentID> amplify_component(Component &component, int min_score);
};


#endif //SRC_READCLUSTERINGENGINE_H
