#include <set>

#include "../common/Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H


struct InClusterReadData {
    std::string header;
    CategoryID category_id;
};


class GenomeReadCluster {
public:
    explicit GenomeReadCluster(ClusterID id, InClusterReadData &read_data, std::set<KmerID> &characteristic_kmers);

    ClusterID reference_id = 0;
    std::vector<InClusterReadData> reads;
    std::set<KmerID> characteristic_kmer_ids;

    void absorb(GenomeReadCluster &cluster);
    std::set<CategoryID> categories;
    std::string consistency();
    uint64_t size();
};


#endif //SRC_GENOMEREADCLUSTER_H
