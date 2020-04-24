#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <experimental/filesystem>
#include <map>
#include <queue>

#include "SequenceRecordIterator.h"
#include "GenomeReadCluster.h"
#include "../lib/BloomFilter.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

typedef std::vector<std::vector<ClusterID>> KmerClusterIndex;
typedef tsl::robin_map<KmerID, std::vector<ClusterID>> IndexRemovalMap;
typedef std::vector<ClusterID > IDComponent;


struct ClusterConnection{
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    ConnectionScore score;
    bool is_good;

    bool operator < (const ClusterConnection& conn) const
    {
        return (this->score < conn.score);
    }
    bool operator > (const ClusterConnection& conn) const
    {
        return (this->score > conn.score);
    }
};


class ReadClusteringEngine {
protected:
    SequenceRecordIterator* reader;
    bloom::BloomFilter<Kmer>* kmers;
    int k;

    void construct_indices_thread();
    int construct_indices();

    ClusterIndex cluster_index;
    KmerIndex kmer_index;
    KmerClusterIndex kmer_cluster_index;
    tsl::robin_map<std::string, CategoryID> read_category_map;

    void get_connections_thread(std::vector<ClusterID> &cluster_indices, std::vector<ClusterConnection> &accumulator, tsl::robin_set<ClusterID> &restricted, int &index);
    std::vector<ClusterConnection> get_connections(std::vector<ClusterID> &cluster_ids, tsl::robin_set<ClusterID> &restricted);

    void kmer_cluster_index_update(IndexRemovalMap::iterator &removal_it, IndexRemovalMap::iterator &removal_end);
    void merge_clusters_thread(std::queue<IDComponent> &component_queue, IndexRemovalMap &for_removal);
    int merge_clusters(const tsl::robin_map<ClusterID, IDComponent> &components);

    int clustering_round();
    int clustering_round(std::vector<ClusterID> &cluster_ids, tsl::robin_set<ClusterID> &restricted);
public:
    void run_clustering();
    std::map<ClusterID, std::string> export_clusters(std::vector<ClusterID> &cluster_ids, std::experimental::filesystem::path &directory_path);

    ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers);

    void print_clusters(int first_n);

    std::string cluster_consistency(GenomeReadCluster *cluster);

    void construct_read_category_map();

    void assemble_clusters(std::vector<ClusterID> &cluster_ids);

    ClusterConnection get_connection(GenomeReadCluster *x, GenomeReadCluster *y);
};

void plot_connection_quality(std::vector<ClusterConnection> &connections);


#endif //SRC_READCLUSTERINGENGINE_H
