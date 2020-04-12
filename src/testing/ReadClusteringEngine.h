#include "../common/BaseReadClusteringEngine.h"

#ifndef SRC_TESTREADCLUSTERINGENGINE_H
#define SRC_TESTREADCLUSTERINGENGINE_H


class ReadClusteringEngine : public BaseReadClusteringEngine {
protected:
    tsl::robin_map<std::string, CategoryID> read_category_map;
    std::set<Kmer> characteristic_kmers;

    void construct_read_category_map();
    std::string cluster_consistency(GenomeReadCluster* cluster);

public:
    void print_clusters(int first_n);
    void run_clustering() override;

    using BaseReadClusteringEngine::BaseReadClusteringEngine;
};


#endif //SRC_TESTREADCLUSTERINGENGINE_H
