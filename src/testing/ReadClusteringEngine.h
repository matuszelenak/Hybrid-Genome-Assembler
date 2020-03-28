#include <tsl/robin_map.h>

#include "../common/SequenceRecordIterator.h"
#include "KmerAnalysis.h"
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
private:
    ClusterIndex cluster_index;
    KmerIndex kmer_index;
    KmerClusterIndex kmer_cluster_index;

    std::vector<ClusterConnection> get_connections();
    KmerClusterIndex get_index();
    ClusterIndex get_initial_read_clusters(SequenceRecordIterator &reader, KmerOccurrences &characteristic_kmers, int k);
    int clustering_round();
public:
    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerOccurrences &characteristic_kmers, int k);
    void run_clustering();
};


struct InClusterKmerInfo {
    uint32_t sum_of_qualities = 0;
    uint32_t total_count = 0;

    Quality avg_quality(){
        return (Quality)(this->sum_of_qualities / this->total_count);
    };
};


#endif //SRC_READCLUSTERINGENGINE_H
