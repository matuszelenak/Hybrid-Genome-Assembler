#include <iostream>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <queue>

#include "KmerIterator.h"
#include "BaseReadClusteringEngine.h"

std::mutex connections_mut;
std::mutex merge_mut;
std::mutex component_erase_mut;

uint64_t MIN_SCORE = 5;


void BaseReadClusteringEngine::get_connections_thread(std::vector<ClusterID> &cluster_indices, std::pair<int, int> range, std::vector<ClusterConnection> &accumulator){
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
            accumulator.push_back({iter->first, pivot_id, iter->second, cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories});
        }
        connections_mut.unlock();
    }
}


std::vector<ClusterConnection> BaseReadClusteringEngine::get_connections(){
    std::vector<ClusterConnection> connections;

    std::vector<ClusterID> cluster_ids;
    std::transform(cluster_index.begin(), cluster_index.end(), std::back_inserter(cluster_ids), [](ClusterIndex::value_type &val) -> ClusterID { return val.first; });

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];
    int clusters_per_thread = (int)ceil(cluster_ids.size() / (double)num_threads);
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(
                &BaseReadClusteringEngine::get_connections_thread,
                this,
                std::ref(cluster_ids),
                std::make_pair(i * clusters_per_thread, std::min((i + 1) * clusters_per_thread, (int)cluster_ids.size())),
                std::ref(connections)
                );
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

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


void BaseReadClusteringEngine::merge_clusters_thread(std::queue<IDComponent> &component_queue){
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


int BaseReadClusteringEngine::clustering_round(){
    tsl::robin_map<ClusterID, ClusterID>parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++){
        parents[cluster_iter->second->reference_id] = cluster_iter->second->reference_id;
        components[cluster_iter->second->reference_id] = {cluster_iter->second->reference_id};
    }

    std::vector<ClusterConnection> cluster_connections = get_connections();
    int merge_operations = 0;
    for (ClusterConnection &conn : cluster_connections){
        if (!conn.is_good) continue;

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

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&BaseReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return merge_operations;
}
