#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <experimental/filesystem>
#include <map>
#include <queue>

#include "SequenceRecordIterator.h"
#include "GenomeReadCluster.h"
#include "Utils.h"
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
    Platform platform;

    void construct_indices_thread(KmerIndex &kmer_index);
    int construct_indices();

    ClusterIndex cluster_index;
    KmerClusterIndex kmer_cluster_index;
    tsl::robin_map<ReadID, CategoryID> read_category_map;
    std::vector<Kmer> kmer_id_to_kmer;

    void get_connections_thread(ConnectionScore min_score, ConcurrentQueue<ClusterID> &cluster_id_queue, std::vector<ClusterConnection> &accumulator);
    std::vector<ClusterConnection> get_all_connections(ConnectionScore min_score);
    std::vector<ClusterConnection> get_connections(std::vector<ClusterID> &cluster_ids, ConnectionScore min_score);

    void kmer_cluster_index_update(IndexRemovalMap::iterator &removal_it, IndexRemovalMap::iterator &removal_end);
    void merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal, std::vector<ClusterID> &merge_result_ids);
    std::vector<ClusterID> merge_clusters(std::vector<IDComponent> &components);

    std::vector<IDComponent> union_find(std::vector<ClusterConnection> &connections, tsl::robin_set<ClusterIDPair> &restricted);

    std::vector<ClusterID> filter_clusters(const std::function<bool(GenomeReadCluster*)>& func){
        std::vector<ClusterID> result;
        for (auto cluster_pair : cluster_index){
            if (func(cluster_pair.second)){
                result.push_back(cluster_pair.first);
            }
        }
        return result;
    };
public:
    void run_clustering();
    std::map<ClusterID, std::string> export_clusters(std::vector<ClusterID> &cluster_ids, std::experimental::filesystem::path &directory_path);

    ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers, Platform platform);

    void print_clusters(int first_n);

    std::string cluster_consistency(GenomeReadCluster *cluster);

    void construct_read_category_map();

    void assemble_clusters(std::vector<ClusterID> &cluster_ids);
};

void plot_connection_quality(std::vector<ClusterConnection> &connections);


#endif //SRC_READCLUSTERINGENGINE_H
