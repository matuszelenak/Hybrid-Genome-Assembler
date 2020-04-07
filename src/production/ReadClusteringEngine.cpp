#include <iostream>
#include "ReadClusteringEngine.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"

std::mutex index_merge;


void ReadClusteringEngine::construct_indices_thread(){
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        std::set<Kmer> in_read_characteristic;

        while (it.next_kmer()) {
            KmerCount count = filter->get_count(it.current_kmer);
            if (lower_coverage <= count && count <= upper_coverage){
                in_read_characteristic.insert(it.current_kmer);
            }
        }

        if (in_read_characteristic.size() >= 5){
            std::set<KmerID> in_read_characteristic_ids;

            index_merge.lock();

            ClusterID new_cluster_id = cluster_index.size();

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (Kmer kmer : in_read_characteristic){
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer, new_kmer_id));
                if (insert_result.second){
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_characteristic_ids.insert(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].insert(new_cluster_id);
            }

            cluster_index.insert(ClusterIndex::value_type(new_cluster_id, new GenomeReadCluster(new_cluster_id, read->header, in_read_characteristic_ids, read->category_id)));

            index_merge.unlock();
        }
    }
}


int ReadClusteringEngine::construct_indices() {
    reader->reset();
    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::construct_indices_thread, this);
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return 0;
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, int lower_coverage, int upper_coverage) {
    this->k = k;
    this->reader = &read_iterator;
    this->filter = &bf;
    this->lower_coverage = lower_coverage;
    this->upper_coverage = upper_coverage;

    auto r = timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
}

void ReadClusteringEngine::run_clustering(){
    while (clustering_round() > 0){}
    std::cout << fmt::format("Exported {} clusters\n", timeMeasureMemberFunc(&ReadClusteringEngine::export_clusters, this, "Exporting clusters")(10));
}