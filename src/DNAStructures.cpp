#include <fmt/format.h>
#include "DNAStructures.h"


GenomeReadCluster::GenomeReadCluster(uint64_t id, InClusterReadData &read_data, std::map<Kmer, InClusterKmerInfo> &characteristic_kmers) {
    reference_id = id;
    _categories = {read_data.category_flag};
    reads.push_back(read_data);
    this->characteristic_kmers = characteristic_kmers;
}

void GenomeReadCluster::absorb(GenomeReadCluster &cluster) {
    reads.insert(reads.end(), cluster.reads.begin(), cluster.reads.end());
    _categories.insert(cluster.categories().begin(), cluster.categories().end());
    for (auto it = begin(cluster.characteristic_kmers); it != end(cluster.characteristic_kmers); it++){
        if (characteristic_kmers.contains(it->first)){
            characteristic_kmers[it->first].sum_of_qualities += it->second.sum_of_qualities;
        } else {
            characteristic_kmers[it->first] = {it->second.sum_of_qualities};
        }
    }
}

std::set<CategoryFlag> GenomeReadCluster::categories() {
    return _categories;
}

uint64_t GenomeReadCluster::size() {
    return reads.size();
}

std::string GenomeReadCluster::consistency() {
    uint64_t first = 0, second = 0;
    for (const auto& read : reads){
        first += read.category_flag;
        second += !read.category_flag;
    }
    return fmt::format("{}/{}", std::max(first, second), this->size());
}
