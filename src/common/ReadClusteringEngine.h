#include <tsl/robin_map.h>
#include <map>
#include <queue>

#include "SequenceRecordIterator.h"
#include "GenomeReadCluster.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

typedef std::vector<std::set<ClusterID> > KmerClusterIndex;
typedef std::vector<ClusterID > IDComponent;


struct ClusterConnection{
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    ConnectionScore score;

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
    ClusterIndex cluster_index;
    KmerIndex kmer_index;
    KmerClusterIndex kmer_cluster_index;

    void get_connections_thread(std::vector<ClusterID> &cluster_indices, std::pair<int, int> range, std::vector<ClusterConnection> &accumulator);
    std::vector<ClusterConnection> get_connections();

    void kmer_cluster_index_thread(ClusterID cluster_from, ClusterID cluster_to);
    void construct_kmer_cluster_index();

    void cluster_index_thread(SequenceRecordIterator &reader, int k, std::vector<GenomeReadCluster*> &thread_result);
    void construct_cluster_index(SequenceRecordIterator &reader, std::set<Kmer> &characteristic_kmers, int k);

    void merge_clusters_thread(std::queue<IDComponent> &component_queue);
    int clustering_round();
public:
    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator, std::set<Kmer> &characteristic_kmers, int k);

    virtual void run_clustering();
};


#endif //SRC_READCLUSTERINGENGINE_H
