#include "../common/BaseReadClusteringEngine.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H


class ReadClusteringEngine: public BaseReadClusteringEngine {
protected:
    KmerCountingBloomFilter* filter;
    int k, lower_coverage, upper_coverage;

    void construct_indices_thread();
    void construct_indices();
public:
    void run_clustering() override;

    ReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, int lower_coverage, int upper_coverage);
};


#endif //SRC_READCLUSTERINGENGINE_H
