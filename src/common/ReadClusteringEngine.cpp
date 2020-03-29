#include <iostream>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <queue>

#include "KmerIterator.h"
#include "ReadClusteringEngine.h"

std::mutex index_mut;
std::mutex connections_mut;
std::mutex merge_mut;
std::mutex component_erase_mut;
unsigned int num_threads = std::thread::hardware_concurrency();

uint64_t MIN_SCORE = 5;


void ReadClusteringEngine::get_connections_thread(std::vector<ClusterID> &cluster_indices, std::pair<int, int> range, std::vector<ClusterConnection> &accumulator){
    for (int i = range.first; i < range.second; i++){
        ClusterID pivot_id = cluster_indices[i];
        std::map<ClusterID , ConnectionScore > shared_kmer_counts;
        for (KmerID kmer_id : cluster_index[pivot_id]->characteristic_kmer_ids){
            for (ClusterID candidate_id : kmer_cluster_index[kmer_id]){
                shared_kmer_counts.insert( std::map<ClusterID , ConnectionScore >::value_type(candidate_id, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        connections_mut.lock();
        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < MIN_SCORE) continue;
            accumulator.push_back({iter->first, pivot_id, iter->second});
        }
        connections_mut.unlock();
    }
}


std::vector<ClusterConnection> ReadClusteringEngine::get_connections(){
    std::vector<ClusterConnection> connections;

    std::vector<ClusterID> cluster_ids;
    std::transform(cluster_index.begin(), cluster_index.end(), std::back_inserter(cluster_ids), [](ClusterIndex::value_type &val) -> ClusterID { return val.first; });

    std::thread t[num_threads];
    std::cout << "Launching connection computing threads\n";
    int clusters_per_thread = (int)ceil(cluster_ids.size() / (double)num_threads);
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(
                &ReadClusteringEngine::get_connections_thread,
                this,
                std::ref(cluster_ids),
                std::make_pair(i * clusters_per_thread, std::min((i + 1) * clusters_per_thread, (int)cluster_ids.size())),
                std::ref(connections)
                );
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    std::cout << fmt::format("{} connections\n", connections.size());
    std::sort(connections.rbegin(), connections.rend());
    return connections;
}


ClusterID get_parent(ClusterID cluster_id, tsl::robin_map<ClusterID, ClusterID> &parents){
    if (cluster_id == parents[cluster_id]){
        return cluster_id;
    }
    parents[cluster_id] = get_parent(parents[cluster_id], parents);
    return parents[cluster_id];
}


void ReadClusteringEngine::merge_clusters_thread(std::queue<IDComponent> &component_queue){
    while (true){
        merge_mut.lock();
        if (component_queue.empty()){
            merge_mut.unlock();
            return;
        }
        IDComponent component = component_queue.front();
        component_queue.pop();
        merge_mut.unlock();

        //TODO faster merging using minheap?
        GenomeReadCluster* survivor = cluster_index[component[0]];
        for (int i = 1; i < component.size(); i++){
            survivor->absorb(*cluster_index[component[i]]);
        }

        component_erase_mut.lock();
        for (int i = 1; i < component.size(); i++){

            for (KmerID kmer_id : cluster_index[component[i]]->characteristic_kmer_ids){
                kmer_cluster_index[kmer_id].erase(component[i]);
                kmer_cluster_index[kmer_id].insert(survivor->reference_id);
            }
            cluster_index.erase(component[i]);
        }
        component_erase_mut.unlock();
    }
}


int ReadClusteringEngine::clustering_round(){
    tsl::robin_map<ClusterID, ClusterID>parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++){
        parents[cluster_iter->second->reference_id] = cluster_iter->second->reference_id;
        components[cluster_iter->second->reference_id] = {cluster_iter->second->reference_id};
    }

    std::vector<ClusterConnection> cluster_connections = get_connections();
    int merge_operations = 0;
    for (ClusterConnection &conn : cluster_connections){
        ClusterID parent_x = get_parent(conn.cluster_x_id, parents);
        ClusterID parent_y = get_parent(conn.cluster_y_id, parents);
        if (parent_x == parent_y) continue;

        ClusterID bigger, smaller;
        if (components[parent_x].size() > components[parent_y].size()){
            bigger = parent_x;
            smaller = parent_y;
        } else {
            bigger = parent_y;
            smaller = parent_x;
        }

        for (ClusterID cluster_id : components[smaller]){
            parents[cluster_id] = bigger;
        }
        components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
        components.erase(smaller);

        merge_operations++;
    }

    // Now do the actual merging of clusters in threads
    std::queue<IDComponent> component_queue;
    for (auto component_it = begin(components); component_it != end(components); component_it++){
        if (component_it->second.size() > 1){
            component_queue.push(component_it->second);
        }
    }

    std::thread t[num_threads];
    std::cout << "Launching merging threads\n";
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();


    return merge_operations;
}

void ReadClusteringEngine::run_clustering(){
    size_t num_of_components = cluster_index.size();
    std::cout << fmt::format("{} clusters", num_of_components) << std::endl;
    int merge_operations;
    while ((merge_operations = clustering_round()) > 0){
        std::cout << fmt::format("Performed {} merge operations\n", merge_operations);
    }
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

    int clusters_per_thread = (int)ceil(cluster_index.size() / (double)num_threads);
    std::thread t[num_threads];
    std::cout << "Launching index computing threads\n";
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::kmer_cluster_index_thread, this, i * clusters_per_thread, std::min((i + 1) * clusters_per_thread, (int) cluster_index.size()));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();
}


void ReadClusteringEngine::cluster_index_thread(SequenceRecordIterator &reader, int k, std::vector<GenomeReadCluster*> &thread_result){
    std::optional<GenomeReadData> read;
    while ((read = reader.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        std::set<KmerID> in_read_characteristic;

        while (it.next_kmer()) {
            if (kmer_index.contains(it.current_kmer)){
                in_read_characteristic.insert(kmer_index.at(it.current_kmer));
            }
        }

        if (!in_read_characteristic.empty()){
            thread_result.push_back(new GenomeReadCluster(0, read->header, in_read_characteristic));
        }
    }
}


void ReadClusteringEngine::construct_cluster_index(SequenceRecordIterator &reader, std::set<Kmer> &characteristic_kmers, int k) {
    // We will no longer be interested in the semantics of a kmer, so we can just convert it to an ID
    KmerID kmer_id = 0;
    for (Kmer kmer : characteristic_kmers){
        kmer_index[kmer] = kmer_id++;
    }

    std::cout << "Launching the calculation of initial clusters\n";
    std::thread t[num_threads];
    std::vector<std::vector<GenomeReadCluster*> >partial_results(num_threads, std::vector<GenomeReadCluster*>());
    reader.reset();
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&ReadClusteringEngine::cluster_index_thread, this, std::ref(reader), k, std::ref(partial_results[i]));
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

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, std::set<Kmer> &characteristic_kmers, int k) {
    construct_cluster_index(read_iterator, characteristic_kmers, k);
    construct_kmer_cluster_index();
}
