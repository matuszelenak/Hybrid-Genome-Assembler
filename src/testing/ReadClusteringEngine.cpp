#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <fmt/format.h>

#include "ReadClusteringEngine.h"
#include "../common/KmerIterator.h"

namespace algo = boost::algorithm;

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
    construct_read_category_map();

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
