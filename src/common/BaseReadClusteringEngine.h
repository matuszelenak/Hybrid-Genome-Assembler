#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <map>
#include <queue>

#include "SequenceRecordIterator.h"
#include "GenomeReadCluster.h"
#include "KmerCountingBloomFilter.h"

#ifndef SRC_BASEREADCLUSTERINGENGINE_H
#define SRC_BASEREADCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

typedef std::vector<std::set<ClusterID> > KmerClusterIndex;
typedef std::vector<ClusterID > IDComponent;

typedef tsl::robin_set<std::pair<ClusterID, ClusterID>, boost::hash<std::pair<ClusterID, ClusterID>>> ProcessedClusters;


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


class BaseReadClusteringEngine {
protected:
    SequenceRecordIterator* reader;
    KmerCountingBloomFilter* filter;
    int k, lower_coverage, upper_coverage;

    void construct_indices_thread();
    int construct_indices();

    ClusterIndex cluster_index;
    KmerIndex kmer_index;
    KmerClusterIndex kmer_cluster_index;

    void get_connections_thread(std::vector<ClusterID> &cluster_indices, std::pair<int, int> range, std::vector<ClusterConnection> &accumulator, ProcessedClusters &processed);
    std::vector<ClusterConnection> get_connections();

    void merge_clusters_thread(std::queue<IDComponent> &component_queue);
    int merge_clusters(const tsl::robin_map<ClusterID, IDComponent> &components);

    int clustering_round();
public:
    virtual void run_clustering() = 0;
    int export_clusters(int min_size);

    BaseReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, int lower_coverage, int upper_coverage);
};

void plot_connection_quality(std::vector<ClusterConnection> &connections);


#endif //SRC_BASEREADCLUSTERINGENGINE_H
