#include <set>
#include <fmt/format.h>

#include "Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H

typedef std::pair<uint32_t, bool> Endpoint;


class GenomeReadCluster {
public:
    GenomeReadCluster(ClusterID id, GenomeReadData &initial_read, std::vector<KmerID> &discriminative_kmer_ids, CategoryID category){
        reference_id = id;
        read_headers.push_back(initial_read.header);
        this->discriminative_kmer_ids = discriminative_kmer_ids;
        categories.insert(category);
        endpoints = {{initial_read.start, true}, {initial_read.end, false}};
    };

    ClusterID reference_id = 0;
    std::vector<std::string> read_headers;
    std::vector<KmerID> discriminative_kmer_ids;
    std::set<CategoryID> categories;
    std::vector<Endpoint> endpoints;

    [[nodiscard]] uint64_t size() const{ return read_headers.size(); };

    std::string repr(){
        std::string cats;
        for (auto cat : categories){
            cats += std::to_string(cat);
            cats += " ";
        }
        return fmt::format("#{} size {} categories {}", reference_id, size(), cats);
    }

    std::string intervals(){
        std::string result;
        for (int i = 0; i < endpoints.size(); i+=2){
            result += (std::to_string(endpoints[i].first) + "-" + std::to_string(endpoints[i + 1].first)) + " ";
        }
        return result;
    }
};


#endif //SRC_GENOMEREADCLUSTER_H
