#include <vector>

#include "DNAStructures.h"
#include "SequenceRecordIterator.h"

#ifndef SRC_READCLUSTERING_H
#define SRC_READCLUSTERING_H

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

struct ClusterConnection{
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    uint64_t score;
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

void run_clustering(ClusterIndex &clusters);
ClusterIndex get_initial_read_clusters(SequenceRecordIterator &reader, int k, KmerOccurrences &characteristic_kmers);

#endif //SRC_READCLUSTERING_H
