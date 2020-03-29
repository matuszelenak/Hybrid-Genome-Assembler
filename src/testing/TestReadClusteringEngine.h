#include "../common/ReadClusteringEngine.h"

#ifndef SRC_TESTREADCLUSTERINGENGINE_H
#define SRC_TESTREADCLUSTERINGENGINE_H


class TestReadClusteringEngine : public ReadClusteringEngine{
protected:
    tsl::robin_map<std::string, CategoryID> read_category_map;

    void construct_read_category_map(SequenceRecordIterator &read_iterator);
    std::string cluster_consistency(GenomeReadCluster* cluster);
public:
    TestReadClusteringEngine(SequenceRecordIterator &read_iterator, std::set<Kmer> &characteristic_kmers, int k);
    void print_clusters(int first_n);
    void run_clustering() override;
};


#endif //SRC_TESTREADCLUSTERINGENGINE_H
