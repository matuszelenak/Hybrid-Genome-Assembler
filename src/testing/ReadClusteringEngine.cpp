#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <boost/thread.hpp>
#include <fmt/format.h>

#include "ReadClusteringEngine.h"
#include "../common/KmerIterator.h"

namespace algo = boost::algorithm;

std::mutex index_mut;

std::string ReadClusteringEngine::cluster_consistency(GenomeReadCluster *cluster) {
    std::map<CategoryID, int> category_counts;
    for (const auto &read_header : cluster->read_headers) {
        category_counts.insert(std::map<CategoryID, int>::value_type(read_category_map[read_header], 0)).first->second += 1;
    }
    std::vector<std::string> category_count_vector;
    std::transform(
            category_counts.begin(),
            category_counts.end(),
            std::back_inserter(category_count_vector),
            [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
    );
    return algo::join(category_count_vector, "/");
}


void ReadClusteringEngine::print_clusters(int first_n) {
    std::vector<GenomeReadCluster *> cluster_pointers;
    std::transform(
            cluster_index.begin(),
            cluster_index.end(),
            std::back_inserter(cluster_pointers),
            [](std::pair<const ClusterID, GenomeReadCluster *> &p) -> GenomeReadCluster * { return p.second; }
    );

    sort(cluster_pointers.rbegin(), cluster_pointers.rend(), [](GenomeReadCluster *x, GenomeReadCluster *y) -> bool { return x->size() < y->size(); });

    int iterate_first = first_n;
    if (first_n == -1) iterate_first = cluster_pointers.size();
    for (auto cluster : cluster_pointers) {
        std::cout << cluster_consistency(cluster) << " ";

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << std::endl;
}

void ReadClusteringEngine::run_clustering() {
    size_t num_of_components = cluster_index.size();
    std::cout << fmt::format("{} clusters", num_of_components) << std::endl;
    int merge_operations;
    while ((merge_operations = clustering_round()) > 0) {
        std::cout << fmt::format("Performed {} merge operations\n", merge_operations);
        print_clusters(200);
    }
    print_clusters(1000000);
}


void ReadClusteringEngine::construct_read_category_map() {
    this->reader->reset();
    std::optional<GenomeReadData> read;
    while ((read = this->reader->get_next_record()) != std::nullopt) {
        read_category_map.insert(tsl::robin_map<std::string, CategoryID>::value_type(read->header, read->category_id));
    }
}


void ReadClusteringEngine::cluster_index_thread(std::vector<GenomeReadCluster*> &thread_result){
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        std::set<KmerID> in_read_characteristic;

        while (it.next_kmer()) {
            if (kmer_index.contains(it.current_kmer)){
                in_read_characteristic.insert(kmer_index.at(it.current_kmer));
            }
        }

        if (!in_read_characteristic.empty()){
            thread_result.push_back(new GenomeReadCluster(0, read->header, in_read_characteristic, read->category_id));
        }
    }
}


void ReadClusteringEngine::construct_cluster_index() {
    // We will no longer be interested in the semantics of a kmer, so we can just convert it to an ID
    KmerID kmer_id = 0;
    for (Kmer kmer : characteristic_kmers){
        kmer_index[kmer] = kmer_id++;
    }

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];
    std::vector<std::vector<GenomeReadCluster*> >partial_results(num_threads, std::vector<GenomeReadCluster*>());
    reader->reset();
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::cluster_index_thread, this, std::ref(partial_results[i]));
    }

    for (int i = 0; i < num_threads; ++i) t[i].join();

    ClusterID cluster_counter = 0;
    for (const auto& partial_clusters: partial_results){
        for (auto cluster : partial_clusters){
            cluster->reference_id = cluster_counter;
            cluster_index[cluster_counter] = cluster;
            cluster_counter++;
        }
    }
    std::cout << fmt::format("Created {} clusters\n", cluster_counter);
}

void ReadClusteringEngine::kmer_cluster_index_thread(ClusterID cluster_from, ClusterID cluster_to){
    for (ClusterID cluster_id = cluster_from; cluster_id < cluster_to; cluster_id++){
        GenomeReadCluster* cluster = cluster_index[cluster_id];
        std::map<KmerID, std::set<ClusterID > > partial;
        for (KmerID kmer_id : cluster->characteristic_kmer_ids){
            partial.insert(std::map<KmerID, std::set<ClusterID > >::value_type(kmer_id, {})).first->second.insert(cluster_id);
        }

        index_mut.lock();
        for (auto it = begin(partial); it != end(partial); it++){
            kmer_cluster_index[it->first].insert(it->second.begin(), it->second.end());
        }
        index_mut.unlock();
    }
}


void ReadClusteringEngine::construct_kmer_cluster_index(){
    kmer_cluster_index = KmerClusterIndex(kmer_index.size(), std::set<ClusterID>());

    unsigned int num_threads = std::thread::hardware_concurrency();
    int clusters_per_thread = (int)ceil(cluster_index.size() / (double)num_threads);
    std::thread t[num_threads];
    std::cout << "Launching index computing threads\n";
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::kmer_cluster_index_thread, this, i * clusters_per_thread, std::min((i + 1) * clusters_per_thread, (int) cluster_index.size()));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();
}



ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, std::set<Kmer> &characteristic_kmers, int k){
    this->reader = &read_iterator;
    this->k = k;
    this->characteristic_kmers = characteristic_kmers;

    construct_cluster_index();
    construct_kmer_cluster_index();
    construct_read_category_map();
}