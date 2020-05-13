#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <filesystem>
#include <map>
#include <queue>

#include "../common/SequenceRecordIterator.h"
#include "../common/Utils.h"

#include "../lib/BloomFilter.h"

#include "ReadComponent.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef std::unordered_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef tsl::robin_map<ComponentID, ReadComponent *> ComponentIndex;

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
    int core_component_min_size = 30;
    double core_forming_connections = 0.2;
    ConnectionScore enrichment_connections_min_score = 20;
    ConnectionScore tail_amplification_min_score = 40;
};

class ReadClusteringEngine {
protected:
    ReadClusteringConfig config;
    SequenceRecordIterator *reader;

    ComponentIndex component_index;
    KmerComponentIndex kmer_component_index;

    tsl::robin_map<ReadID, uint32_t> read_lengths;
    tsl::robin_map<ReadID, tsl::robin_map<KmerID, int>> kmer_positions;

    // DEBUG
    bool debug = false;
    tsl::robin_map<ReadID, std::pair<uint32_t, uint32_t> > read_intervals;

    std::vector<Interval> get_read_coverage_intervals(std::vector<ReadID> &read_ids);

    void export_spanning_tree(SpanningTree &tree, std::pair<std::vector<ReadID>, std::vector<ReadID>> &tails);

    void plot_read_overlaps(std::vector<ComponentConnection> &connections);

    void plot_overlap_differences(std::vector<ComponentConnection> &connections);

    void print_components(std::vector<ComponentID> &ids);
    // END DEBUG

    int construct_indices(std::unordered_set<Kmer> &discriminative_kmers, int k);

    std::vector<ComponentConnection> get_all_connections();

    std::vector<ComponentConnection> get_connections(Component &component, ConnectionScore min_score);

    std::vector<ComponentID> merge_components(std::vector<Component> &components);

    std::vector<std::pair<Component, SpanningTree>> union_find(std::vector<ComponentConnection> &connections, std::set<ComponentID> &restricted, int min_component_size);

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
