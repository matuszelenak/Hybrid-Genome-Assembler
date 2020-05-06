#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <experimental/filesystem>
#include <map>
#include <queue>

#include "../common/SequenceRecordIterator.h"
#include "../common/Utils.h"

#include "../lib/BloomFilter.h"

#include "GenomeReadCluster.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

typedef std::vector<std::vector<ClusterID>> KmerClusterIndex;
typedef tsl::robin_map<KmerID, std::vector<ClusterID>> IndexRemovalMap;
typedef std::vector<ClusterID > IDComponent;

typedef std::vector<std::pair<ClusterID, ClusterID> > SpanningTree;
typedef tsl::robin_map<ClusterID, tsl::robin_set<ClusterID> > NeighborMap;


struct ClusterConnection{
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    ConnectionScore score;
    bool is_good;

    bool operator < (const ClusterConnection& conn) const
    {
        return (this->score < conn.score);
    }
};

void plot_connection_quality(std::vector<ClusterConnection> &connections);

class ReadClusteringEngine {
protected:
    SequenceRecordIterator* reader;
    bloom::BloomFilter<Kmer>* kmers = nullptr;
    int k = 0;

    ClusterIndex cluster_index;
    KmerClusterIndex kmer_cluster_index;

    // DEBUG
    tsl::robin_map<ReadID, std::pair<uint32_t, uint32_t> > read_intervals;

    void construct_indices_thread(KmerIndex &kmer_index);
    int construct_indices();

    void get_connections_thread(ConnectionScore min_score, ConcurrentQueue<ClusterID> &cluster_id_queue, std::vector<ClusterConnection> &accumulator);
    std::vector<ClusterConnection> get_connections(std::vector<ClusterID> &cluster_ids, ConnectionScore min_score);

    void kmer_cluster_index_update(ConcurrentQueue<IndexRemovalMap::value_type> &removal_list_queue);
    void merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal, std::vector<ClusterID> &merged_ids);
    std::vector<ClusterID> merge_clusters(std::vector<IDComponent> &components);

    std::vector<std::pair<IDComponent, SpanningTree>> union_find(std::vector<ClusterConnection> &connections, std::set<ClusterID> &restricted, int min_component_size);

    std::pair<std::vector<ReadID>, std::vector<ReadID>> get_spanning_tree_boundary_reads(SpanningTree &tree);
    std::pair<std::vector<KmerID>, std::vector<KmerID>> get_core_cluster_kmer_tails(SpanningTree &tree);
    std::vector<ClusterConnection> get_core_cluster_connections(std::vector<std::pair<IDComponent, SpanningTree>> &components_with_trees);

    std::map<ClusterID, std::vector<std::string>> assemble_clusters(std::vector<ClusterID> &cluster_ids);

    std::vector<Interval> get_read_coverage_intervals(std::vector<ReadID> &read_ids);

    std::vector<ClusterConnection> filter_connections(std::vector<ClusterConnection> &original, const std::function<bool(ClusterConnection &)> &func) {
        std::vector<ClusterConnection> result;
        for (auto conn : original) {
            if (func(conn)) result.push_back(conn);
        }
        return result;
    }

    void remove_merged_clusters(){
        for(auto it = cluster_index.begin();it != cluster_index.end();){
            if (it->second->size() == 0){
                free(it->second);
                cluster_index.erase(it++);
            } else ++it;
        }
    }

public:
    void export_spanning_tree(SpanningTree &tree);
    void plot_read_overlaps(std::vector<ClusterConnection> &connections);
    void print_clusters(std::vector<ClusterID> &ids);

    std::map<ClusterID, std::string> export_clusters(std::vector<ClusterID> &cluster_ids, std::experimental::filesystem::path &directory_path);

    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator);
    ~ReadClusteringEngine();

    void set_kmers(tsl::robin_set<Kmer> &_kmers, int k);
    void set_kmers(bloom::BloomFilter<Kmer>* _kmers, int k);
    std::vector<ClusterID> run_clustering();
};


#endif //SRC_READCLUSTERINGENGINE_H
