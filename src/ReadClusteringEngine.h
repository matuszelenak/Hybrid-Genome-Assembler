#include "SequenceRecordIterator.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

//typedef std::map<ClusterID, Quality> InClusterKmerQuality;
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
    int clustering_round();
    KmerClusterIndex get_index();
    ClusterIndex get_initial_read_clusters(SequenceRecordIterator &reader, KmerOccurrences &characteristic_kmers, int k);
public:
    explicit ReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerOccurrences &characteristic_kmers, int k);
};


#endif //SRC_READCLUSTERINGENGINE_H
