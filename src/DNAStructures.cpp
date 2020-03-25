#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>
#include <iostream>
#include "DNAStructures.h"


namespace algo = boost::algorithm;


GenomeReadCluster::GenomeReadCluster(ClusterID id, InClusterReadData &read_data, std::set<KmerID> &characteristic_kmers) {
    reference_id = id;
    categories.insert(read_data.category_id);
    reads.push_back(read_data);
    this->characteristic_kmer_ids.insert(characteristic_kmers.begin(), characteristic_kmers.end());
}

void GenomeReadCluster::absorb(GenomeReadCluster &cluster) {
    reads.insert(reads.end(), cluster.reads.begin(), cluster.reads.end());
    categories.insert(cluster.categories.begin(), cluster.categories.end());
    characteristic_kmer_ids.insert(cluster.characteristic_kmer_ids.begin(), cluster.characteristic_kmer_ids.end());
}

uint64_t GenomeReadCluster::size() {
    return reads.size();
}

std::string GenomeReadCluster::consistency() {
    std::map<CategoryID , int> category_counts;
    for (const auto& read : reads){
        category_counts.insert(std::map<CategoryID , int>::value_type(read.category_id, 0)).first->second += 1;
    }
    std::vector<std::string>category_count_vector;
    std::transform(
            category_counts.begin(),
            category_counts.end(),
            std::back_inserter(category_count_vector),
            [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
            );
    return algo::join(category_count_vector, "/");
}
