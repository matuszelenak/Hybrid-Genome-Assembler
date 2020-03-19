#include <unordered_map>
#include <vector>
#include <set>
#include <iostream>

#include "GenomeReadCluster.h"
#include "SequenceReader.h"

using namespace std;

typedef std::unordered_map<Kmer, std::set<uint64_t>> KmerIndex;
typedef std::vector<uint64_t> ClusterConnection;


bool compare_connections(const ClusterConnection &a, const ClusterConnection &b) { return a[0] < b[0]; }


KmerIndex get_kmer_index(std::vector<GenomeReadCluster> &clusters) {
    KmerIndex index;
    for (auto &cluster : clusters) {
        for (auto kmer : cluster.characteristic_kmers) {
            if (!index.contains(kmer)) {
                index[kmer] = {cluster.reference_id};
            } else {
                index[kmer].insert(cluster.reference_id);
            }
        }
    }
    return index;
}

std::vector<ClusterConnection> get_cluster_connection(std::vector<GenomeReadCluster> &clusters, KmerIndex &index) {
    std::vector<ClusterConnection> connections;

    for (GenomeReadCluster &pivot_cluster : clusters) {
        unordered_map<uint64_t, uint64_t> shared_kmer_counts;
        for (Kmer kmer : pivot_cluster.characteristic_kmers) {
            for (uint64_t merge_candidate : index[kmer]) {
                if (!shared_kmer_counts.contains(merge_candidate)) {
                    shared_kmer_counts[merge_candidate] = 1;
                } else {
                    shared_kmer_counts[merge_candidate] += 1;
                }
            }
        }

        shared_kmer_counts.erase(pivot_cluster.reference_id);

        unordered_map<uint64_t, uint64_t>::iterator iter;
        for (iter = shared_kmer_counts.begin(); iter != shared_kmer_counts.end(); iter++) {
            connections.push_back({iter->second, iter->first, pivot_cluster.reference_id});
        }
    }

    sort(connections.rbegin(), connections.rend(), compare_connections);

    return connections;
}

vector<GenomeReadCluster> clustering_round(vector<GenomeReadCluster> &clusters, KmerIndex &index) {
    unordered_map<uint64_t, GenomeReadCluster> id_to_cluster;
    for (auto &cluster: clusters) { id_to_cluster[cluster.reference_id] = cluster; }

    vector<ClusterConnection> cluster_connections = get_cluster_connection(clusters, index);
    for (ClusterConnection &connection : cluster_connections) {
        if (!(id_to_cluster.contains(connection[1]) && id_to_cluster.contains(connection[2]))) continue;

        GenomeReadCluster cluster_x = id_to_cluster[connection[1]];
        GenomeReadCluster cluster_y = id_to_cluster[connection[2]];

        uint64_t bigger_id, smaller_id;
        if (cluster_x.size() > cluster_y.size()) {
            bigger_id = cluster_x.reference_id;
            smaller_id = cluster_y.reference_id;
        } else {
            bigger_id = cluster_y.reference_id;
            smaller_id = cluster_x.reference_id;
        }

        for (Kmer kmer : id_to_cluster[smaller_id].characteristic_kmers) {
            index[kmer].erase(smaller_id);
            index[kmer].insert(bigger_id);
        }

        id_to_cluster[bigger_id].absorb(id_to_cluster[smaller_id]);
        id_to_cluster.erase(smaller_id);
    }

    unordered_map<uint64_t, GenomeReadCluster>::iterator it;
    vector<GenomeReadCluster> result;
    for (it = id_to_cluster.begin(); it != id_to_cluster.end(); it++){
        result.push_back(it->second);
    }
    return result;
}


void print_clusters(vector<GenomeReadCluster> &clusters){
    auto sort_by_size = [](GenomeReadCluster &x, GenomeReadCluster &y){
        return x.size() < y.size();
    };

    int iterate_first = 200;
    sort(clusters.rend(), clusters.rbegin(), sort_by_size);
    for (auto cluster : clusters){
        cout << cluster.consistency() << " ";

        iterate_first--;
        if (iterate_first == 0) break;
    }
    cout << endl;
}


void run_clustering(vector<GenomeReadCluster> &clusters){
    KmerIndex index = get_kmer_index(clusters);
    size_t num_of_components = clusters.size();
    while (true){
        clusters = clustering_round(clusters, index);
        print_clusters(clusters);

        if (clusters.size() == num_of_components) break;

        num_of_components = clusters.size();
    }
}
