#include <vector>

#include "DNAStructures.h"
#include "ReadDataLoader.h"

#ifndef SRC_READCLUSTERING_H
#define SRC_READCLUSTERING_H

struct ClusterConnection{
    uint64_t cluster_x_id;
    uint64_t cluster_y_id;
    uint64_t score;

    bool operator < (const ClusterConnection& conn) const
    {
        return (this->score < conn.score);
    }
    bool operator > (const ClusterConnection& conn) const
    {
        return (this->score > conn.score);
    }
};

void run_clustering(std::vector<GenomeReadCluster> &clusters);
std::vector<GenomeReadCluster> get_initial_read_clusters(ReadDataLoader &reader, int k, KmerOccurrences &characteristic_kmers);

#endif //SRC_READCLUSTERING_H
