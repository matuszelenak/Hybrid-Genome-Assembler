#include <unordered_map>
#include <vector>
#include <set>
#include <iostream>
#include <fmt/format.h>

#include "ReadClustering.h"
#include "KmerIterator.h"

typedef std::map<uint64_t, Quality> InClusterKmerQuality;
typedef std::unordered_map<Kmer, InClusterKmerQuality> KmerIndex;


KmerIndex get_kmer_index(std::vector<GenomeReadCluster> &clusters) {
    KmerIndex index;
    for (auto &cluster : clusters) {
        std::cout << fmt::format("Index for cluster #{}\n", cluster.reference_id);
        std::map<Kmer, InClusterKmerInfo>::iterator it;
        for (it = cluster.characteristic_kmers.begin(); it != cluster.characteristic_kmers.end(); it++) {
            if (index.contains(it->first)) {
                index[it->first][cluster.reference_id] = it->second.avg_quality(cluster.size());
            } else {
                index[it->first] = {{cluster.reference_id, it->second.avg_quality(cluster.size())}};
            }
        }
    }

    for (auto iter = begin(index); iter != end(index); iter++){
        std::cout << fmt::format("Kmer {} ({}): ", iter->first, iter->second.size());
        for (auto k_iter = begin(iter->second); k_iter != end(iter->second); k_iter++){
            std::cout << fmt::format("#{}", k_iter->first);
        }
        std::cout << std::endl;
    }
    return index;
}

std::vector<ClusterConnection> get_cluster_connection(std::vector<GenomeReadCluster> &clusters, KmerIndex &index) {
    std::vector<ClusterConnection> connections;

    for (auto &pivot_cluster : clusters) {
        std::unordered_map<uint64_t, uint64_t> shared_kmer_counts;

        for (auto it = begin(pivot_cluster.characteristic_kmers); it != end(pivot_cluster.characteristic_kmers); it++) {
            for (auto candidates_it = begin(index[it->first]); candidates_it != end(index[it->first]); candidates_it++){
                // TODO metric based on qualities (candidates_it->second)
                if (!shared_kmer_counts.contains(candidates_it->first)) {
                    shared_kmer_counts[candidates_it->first] = 1;
                } else {
                    shared_kmer_counts[candidates_it->first] += 1;
                }
            }
        }

        shared_kmer_counts.erase(pivot_cluster.reference_id);

        std::unordered_map<uint64_t, uint64_t>::iterator iter;
        for (iter = shared_kmer_counts.begin(); iter != shared_kmer_counts.end(); iter++) {
            connections.push_back({iter->second, iter->first, pivot_cluster.reference_id});
        }
    }

    std::cout << 2 << std::endl;
    std::cout << fmt::format("{} connections\n", connections.size());

    std::sort(connections.rbegin(), connections.rend());

    auto until = connections.cbegin();
    std::advance(until, 100);
    for (auto it = connections.cbegin(); it != until; it++){
        std::cout << fmt::format("Score {} X {} Y {}\n", it->score, it->cluster_x_id, it->cluster_y_id);
    }

    return connections;
}

std::vector<GenomeReadCluster> clustering_round(std::vector<GenomeReadCluster> &clusters, KmerIndex &index) {
    std::unordered_map<uint64_t, GenomeReadCluster*> id_to_cluster;
    for (auto cluster: clusters) { id_to_cluster[cluster.reference_id] = &cluster; }

    std::cout << 1 << std::endl;

    std::vector<ClusterConnection> cluster_connections = get_cluster_connection(clusters, index);
    for (ClusterConnection &connection : cluster_connections) {
        if (!(id_to_cluster.contains(connection.cluster_x_id) && id_to_cluster.contains(connection.cluster_y_id))) continue;

        GenomeReadCluster cluster_x = *id_to_cluster[connection.cluster_x_id];
        GenomeReadCluster cluster_y = *id_to_cluster[connection.cluster_y_id];

        uint64_t bigger_id, smaller_id;
        if (cluster_x.size() > cluster_y.size()) {
            bigger_id = cluster_x.reference_id;
            smaller_id = cluster_y.reference_id;
        } else {
            bigger_id = cluster_y.reference_id;
            smaller_id = cluster_x.reference_id;
        }

        GenomeReadCluster bigger = *id_to_cluster[bigger_id];
        GenomeReadCluster smaller = *id_to_cluster[smaller_id];

        bigger.absorb(smaller);

        for (auto it = begin(smaller.characteristic_kmers); it != end(smaller.characteristic_kmers); it++){
            index[it->first].erase(smaller_id);
            index[it->first][bigger_id] = bigger.characteristic_kmers[it->first].avg_quality(bigger.size());
        }

        id_to_cluster.erase(smaller_id);
    }

    std::unordered_map<uint64_t, GenomeReadCluster*>::iterator it;
    std::vector<GenomeReadCluster> result;
    for (it = id_to_cluster.begin(); it != id_to_cluster.end(); it++){
        result.push_back(*it->second);
    }
    return result;
}


void print_clusters(std::vector<GenomeReadCluster> &clusters){
    auto sort_by_size = [](GenomeReadCluster &x, GenomeReadCluster &y){
        return x.size() < y.size();
    };

    int iterate_first = 200;
    sort(clusters.rend(), clusters.rbegin(), sort_by_size);
    for (auto cluster : clusters){
        std::cout << cluster.consistency() << " ";

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << std::endl;
}


void run_clustering(std::vector<GenomeReadCluster> &clusters){
    KmerIndex index = get_kmer_index(clusters);
    size_t num_of_components = clusters.size();
    std::cout << fmt::format("{} clusters", num_of_components) << std::endl;
    while (true){
        clusters = clustering_round(clusters, index);
        print_clusters(clusters);

        if (clusters.size() == num_of_components) break;

        num_of_components = clusters.size();
    }
}

std::vector<GenomeReadCluster> get_initial_read_clusters(ReadDataLoader &reader, int k, KmerOccurrences &characteristic_kmers) {
    std::vector<GenomeReadCluster> clusters;

    std::optional<GenomeReadData> read;
    uint64_t cluster_id = 0;
    while ((read = reader.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k);
        std::optional<std::pair<Kmer, KmerQuality>> kmer_info;

        std::map<Kmer, InClusterKmerInfo> in_read_characteristic;

        while ((kmer_info = it.get_next_kmer()) != std::nullopt) {
            if (characteristic_kmers.contains(kmer_info->first)){
                in_read_characteristic[kmer_info->first] = {static_cast<uint64_t>(kmer_info->second.avg_quality)};
            }
        }

        if (!in_read_characteristic.empty()){
            InClusterReadData data = {read->header, read->category_flag};
            clusters.emplace_back(cluster_id++, data, in_read_characteristic);
        }
    }

    std::cout << fmt::format("{} reads converted to clusters", clusters.size()) << std::endl;

    return clusters;
}
