#include <fmt/format.h>

#include "GenomeReadCluster.h"


GenomeReadCluster::GenomeReadCluster(ClusterID id, std::string &read_header, std::set<KmerID> &characteristic_kmers) {
    reference_id = id;
    read_headers.push_back(read_header);
    this->characteristic_kmer_ids.insert(characteristic_kmers.begin(), characteristic_kmers.end());
}

void GenomeReadCluster::absorb(GenomeReadCluster &cluster) {
    read_headers.insert(read_headers.end(), cluster.read_headers.begin(), cluster.read_headers.end());
    characteristic_kmer_ids.insert(cluster.characteristic_kmer_ids.begin(), cluster.characteristic_kmer_ids.end());
}

uint64_t GenomeReadCluster::size() {
    return read_headers.size();
}
