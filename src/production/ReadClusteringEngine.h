#include "../common/BaseReadClusteringEngine.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H


class ReadClusteringEngine: public BaseReadClusteringEngine {
public:
    void run_clustering() override;
    using BaseReadClusteringEngine::BaseReadClusteringEngine;
};


#endif //SRC_READCLUSTERINGENGINE_H
