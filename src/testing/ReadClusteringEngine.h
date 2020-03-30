#include "../common/BaseReadClusteringEngine.h"

#ifndef SRC_TESTREADCLUSTERINGENGINE_H
#define SRC_TESTREADCLUSTERINGENGINE_H


class ReadClusteringEngine : public BaseReadClusteringEngine {
protected:
    tsl::robin_map<std::string, CategoryID> read_category_map;
    std::set<Kmer> characteristic_kmers;

    void construct_read_category_map();
    std::string cluster_consistency(GenomeReadCluster* cluster);

    void kmer_cluster_index_thread(ClusterID cluster_from, ClusterID cluster_to);
    void construct_kmer_cluster_index();

    void cluster_index_thread(std::vector<GenomeReadCluster*> &thread_result);
    void construct_cluster_index();
public:
    ReadClusteringEngine(SequenceRecordIterator &read_iterator, std::set<Kmer> &characteristic_kmers, int k);
    void print_clusters(int first_n);
    void run_clustering() override;
};


#endif //SRC_TESTREADCLUSTERINGENGINE_H
