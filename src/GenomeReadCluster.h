#include <vector>
#include <set>

#include "structs.h"

#ifndef SRC_GENOMECLUSTER_H
#define SRC_GENOMECLUSTER_H


class GenomeReadCluster {
private:
    std::set<Category> _categories;

public:
    explicit GenomeReadCluster(uint64_t id, GenomeReadData &read_data);

    uint64_t reference_id;
    std::vector<GenomeReadData> reads;
    std::set<Kmer> characteristic_kmers;

    void absorb(GenomeReadCluster &cluster);
    std::set<Category> categories();
    std::string consistency();
    uint64_t size();
};


#endif //SRC_GENOMECLUSTER_H
