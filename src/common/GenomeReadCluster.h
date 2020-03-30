#include <set>

#include "Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H


class GenomeReadCluster {
public:
    explicit GenomeReadCluster(ClusterID id, std::string &read_header, std::set<KmerID> &characteristic_kmers, CategoryID category);

    ClusterID reference_id = 0;
    std::vector<std::string> read_headers;
    std::set<KmerID> characteristic_kmer_ids;

    void absorb(GenomeReadCluster &cluster);
    uint64_t size();

    std::set<CategoryID> categories;
};


#endif //SRC_GENOMEREADCLUSTER_H
