#include <set>
#include <fmt/format.h>

#include "Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H

typedef std::pair<uint32_t, bool> Endpoint;


class GenomeReadCluster {
public:
    GenomeReadCluster(GenomeReadData &initial_read, std::vector<KmerID> &discriminative_kmer_ids){
        this->id = initial_read.id;
        this->discriminative_kmer_ids = discriminative_kmer_ids;
        this->read_ids = {initial_read.id};
        this->categories = {initial_read.category_id};

        if (initial_read.start != 0 || initial_read.end != 0){
            endpoints = {{initial_read.start, true}, {initial_read.end, false}};
        }
    };

    ClusterID id = 0;
    std::vector<ReadID> read_ids;
    std::vector<KmerID> discriminative_kmer_ids;
    std::set<CategoryID> categories;
    std::vector<Endpoint> endpoints;

    [[nodiscard]] uint64_t size() const{ return read_ids.size(); };

    std::string repr(){
        std::string cats;
        for (auto cat : categories){
            cats += std::to_string(cat);
            cats += " ";
        }
        return fmt::format("#{} size {} categories {}", id, size(), cats);
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
