#include <vector>
#include <set>
#include <iostream>
#include <fmt/format.h>
#include <tsl/robin_map.h>
#include <boost/algorithm/string/join.hpp>

#include "ReadClustering.h"
#include "KmerIterator.h"
#include "Utils.h"

typedef std::map<ClusterID, Quality> InClusterKmerQuality;
typedef tsl::robin_map<Kmer, InClusterKmerQuality> KmerIndex;

namespace algo = boost::algorithm;


void plot_connection_quality(std::vector<ClusterConnection> &connections){
    std::map<bool, std::map<ConnectionScore , int>>score_to_matching;
    for (auto &conn : connections){
        score_to_matching[conn.is_good].insert(std::map<ConnectionScore , int>::value_type(conn.score, 0)).first->second += 1;
    }
    std::vector<std::string>good_connections, bad_connections;
    for (auto it = begin(score_to_matching[false]); it != end(score_to_matching[false]); it++){
        bad_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    for (auto it = begin(score_to_matching[true]); it != end(score_to_matching[true]); it++){
        good_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    std::string hist_input = fmt::format("{{{}}}\n{{{}}}", algo::join(good_connections, ", "), algo::join(bad_connections, ", "));
    run_command_with_input("python python_scripts/plot_connections_histogram.py", hist_input);
}


KmerIndex get_kmer_index(ClusterIndex &clusters) {
    KmerIndex index;
    for (auto cluster_iter = begin(clusters); cluster_iter != end(clusters); cluster_iter++) {
        for (auto it = begin(cluster_iter->second->characteristic_kmers); it != end(cluster_iter->second->characteristic_kmers); it++){
            KmerIndex::iterator iter = index.insert(KmerIndex::value_type(it->first, {})).first;
            iter.value()[cluster_iter->first] = it->second.avg_quality(cluster_iter->second->size());
        }
    }

    return index;
}

std::vector<ClusterConnection> get_cluster_connections(ClusterIndex &clusters, KmerIndex &index) {
    std::vector<ClusterConnection> connections;

    for (auto pivot_it = begin(clusters); pivot_it != end(clusters); pivot_it++) {
        std::map<ClusterID , ConnectionScore > shared_kmer_counts;

        for (auto it = begin(pivot_it->second->characteristic_kmers); it != end(pivot_it->second->characteristic_kmers); it++) {
            for (auto candidates_it = begin(index[it->first]); candidates_it != end(index[it->first]); candidates_it++){
                shared_kmer_counts.insert( std::map<ClusterID , ConnectionScore >::value_type(candidates_it->first, 0)).first->second += 1;
            }
        }

        shared_kmer_counts.erase(pivot_it->first);

        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            bool is_good = clusters[pivot_it->first]->categories == clusters[iter->first]->categories;
            connections.push_back({iter->first, pivot_it->first, iter->second, is_good});
        }
    }

    std::cout << fmt::format("{} connections\n", connections.size());

    std::sort(connections.rbegin(), connections.rend());

    return connections;
}

int clustering_round(ClusterIndex &clusters, KmerIndex &index) {
    std::vector<ClusterConnection> cluster_connections = get_cluster_connections(clusters, index);
    plot_connection_quality(cluster_connections);

    int merge_operations = 0;
    for (ClusterConnection &connection : cluster_connections) {
        if (!connection.is_good) continue;

        if (!(clusters.contains(connection.cluster_x_id) && clusters.contains(connection.cluster_y_id))) continue;

        GenomeReadCluster* cluster_x = clusters[connection.cluster_x_id];
        GenomeReadCluster* cluster_y = clusters[connection.cluster_y_id];

        ClusterID bigger_id, smaller_id;
        if (cluster_x->size() > cluster_y->size()) {
            bigger_id = cluster_x->reference_id;
            smaller_id = cluster_y->reference_id;
        } else {
            bigger_id = cluster_y->reference_id;
            smaller_id = cluster_x->reference_id;
        }

        GenomeReadCluster* bigger = clusters[bigger_id];
        GenomeReadCluster* smaller = clusters[smaller_id];

        bigger->absorb(*smaller);

        for (auto it = begin(smaller->characteristic_kmers); it != end(smaller->characteristic_kmers); it++){
            index[it->first].erase(smaller_id);
            index[it->first][bigger_id] = bigger->characteristic_kmers[it->first].avg_quality(bigger->size());
        }

        clusters.erase(smaller_id);
        merge_operations++;
    }

    return merge_operations;
}


void print_clusters(ClusterIndex &clusters){
    auto sort_by_size = [](GenomeReadCluster* x, GenomeReadCluster* y){
        return x->size() < y->size();
    };

    std::vector<GenomeReadCluster*> cluster_pointers;
    std::transform(
            clusters.begin(),
            clusters.end(),
            std::back_inserter(cluster_pointers),
            [](std::pair<const ClusterID , GenomeReadCluster*> &p) -> GenomeReadCluster* { return p.second; }
            );

    sort(cluster_pointers.rbegin(), cluster_pointers.rend(), sort_by_size);

    int iterate_first = 200;
    for (auto cluster : cluster_pointers){
        std::cout << cluster->consistency() << " ";

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << std::endl;
}


void run_clustering(ClusterIndex &clusters){
    KmerIndex index = get_kmer_index(clusters);
    size_t num_of_components = clusters.size();
    std::cout << fmt::format("{} clusters", num_of_components) << std::endl;
    int merge_operations;
    while ((merge_operations = clustering_round(clusters, index)) > 0){
        std::cout << fmt::format("Performed {} merge operations\n", merge_operations);
        print_clusters(clusters);
    }
}

ClusterIndex get_initial_read_clusters(SequenceRecordIterator &reader, int k, KmerOccurrences &characteristic_kmers) {
    ClusterIndex clusters;

    std::optional<GenomeReadData> read;
    ClusterID cluster_id = 0;
    reader.reset();
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
            InClusterReadData data = {read->header, read->category_id};
            clusters[cluster_id] = new GenomeReadCluster(cluster_id, data, in_read_characteristic);
            cluster_id++;
        }
    }

    return clusters;
}
