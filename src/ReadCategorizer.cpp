
#include <optional>
#include <tuple>
#include <algorithm>

#include "ReadCategorizer.h"
#include "SequenceReader.h"
#include "KmerIterator.h"


uint32_t ReadCategorizer::union_metric(std::set< uint32_t> &A, std::set< uint32_t> &B){
    std::set<uint32_t > intersection;
    set_intersection(A.begin(),A.end(),B.begin(),B.end(),
                     std::inserter(intersection,intersection.begin()));
    return (uint32_t)intersection.size();
}


uint64_t ReadCategorizer::get_parent(uint64_t vertex){
    if (parent[vertex] == vertex)
        return vertex;

    parent[vertex] = get_parent(parent[vertex]);
    return parent[vertex];
}


bool ReadCategorizer::unite(uint64_t x, uint64_t y) {
    uint64_t parent_x = get_parent(x);
    uint64_t parent_y = get_parent(y);
    if (parent_x == parent_y)
        return false;

    uint64_t bigger, smaller;
    if (components[parent_x].size() > components[parent_y].size()){
        bigger = parent_x;
        smaller = parent_y;
    }
    else {
        bigger = parent_y;
        smaller = parent_x;
    }

    for (uint64_t vertex : components[smaller]){
        parent[vertex] = bigger;
    }

    components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
    components[smaller].clear();

    kmers_per_component[bigger].insert(kmers_per_component[smaller].begin(), kmers_per_component[smaller].end());
    kmers_per_component[smaller].clear();

    return true;
}


ReadCategorizer::ReadCategorizer(
        const std::string &read_path,
        std::unordered_map<uint64_t,
        uint64_t> &characteristic_kmers,
        int k,
        int final_categories) {
    SequenceReader reader = SequenceReader(read_path);

    // Compress kmers to a sequence of consecutive integers, since we won't care about their semantics anymore
    uint32_t kmer_counter = 0;
    std::unordered_map<uint64_t , uint32_t > kmer_compression_lookup;
    for (auto kmer_data : characteristic_kmers){
        kmer_compression_lookup[kmer_data.first] = kmer_counter;
        ++kmer_counter;
    }

    std::vector<std::set< uint32_t>> kmers_per_read;
    std::optional<GenomeRead> read;
    while ((read = reader.get_next_record()) != std::nullopt){
        KmerIterator it = KmerIterator(*read, k);
        std::optional<uint64_t> kmer_signature;

        std::set< uint32_t> reads_kmers;
        while ((kmer_signature = it.get_next_kmer()) != std::nullopt) {
            if (characteristic_kmers.find(*kmer_signature) != characteristic_kmers.end()){
                reads_kmers.insert(kmer_compression_lookup[*kmer_signature]);
            }
        }
        kmers_per_read.push_back(reads_kmers);
    }

    uint64_t num_of_reads = kmers_per_read.size();

    std::vector< std::tuple<uint64_t, uint64_t, uint64_t > > edges;
    for (uint64_t component_A = 0; component_A < num_of_reads; component_A++){
        for (uint64_t component_B = 0; component_B < num_of_reads; component_B++){
            if (component_A == component_B)
                continue;

            uint64_t edge_weight = union_metric(kmers_per_read[component_A], kmers_per_read[component_B]);
            if (edge_weight == 0)
                continue;

            edges.emplace_back(std::make_tuple(edge_weight, component_A, component_B));
        }
    }
    std::sort(edges.rbegin(), edges.rend());
    kmers_per_component = kmers_per_read;

    for (auto edge : edges){
        if (unite(std::get<1>(edge), std::get<2>(edge))){

        }
    }

    unsigned long size_treshold = 10;

    std::vector<uint64_t > large_component_indices;
    for (uint64_t i = 0; i < num_of_reads; i++){
        if (components[i].size() >= size_treshold ) {
            large_component_indices.push_back(i);
        }
    }

    edges.clear();
    for (uint64_t component_A: large_component_indices){
        for (uint64_t component_B: large_component_indices){
            if (component_A == component_B)
                continue;

            uint64_t edge_weight = union_metric(kmers_per_component[component_A], kmers_per_component[component_B]);
            if (edge_weight == 0)
                continue;

            edges.emplace_back(std::make_tuple(edge_weight, component_A, component_B));
        }
    }

    uint64_t num_of_components = large_component_indices.size();
    std::sort(edges.rbegin(), edges.rend());
    for (auto edge : edges){
        if (unite(std::get<1>(edge), std::get<2>(edge)))
            num_of_components--;

        if (num_of_components == final_categories)
            break;
    }
}