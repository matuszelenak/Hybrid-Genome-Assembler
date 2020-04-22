#include <set>
#include <fmt/format.h>

#include "Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H


class GenomeReadCluster {
public:
    GenomeReadCluster(ClusterID id, std::string &read_header, std::vector<KmerID> &discriminative_kmer_ids, CategoryID category){
        reference_id = id;
        read_headers.push_back(read_header);
        this->discriminative_kmer_ids = discriminative_kmer_ids;
        categories.insert(category);
    };

    ClusterID reference_id = 0;
    std::vector<std::string> read_headers;
    std::vector<KmerID> discriminative_kmer_ids;
    std::set<CategoryID> categories;

    [[nodiscard]] uint64_t size() const{ return read_headers.size(); };

    std::string repr(){
        std::string cats;
        for (auto cat : categories){
            cats += std::to_string(cat);
            cats += " ";
        }

        return fmt::format("#{} size {} categories {}", reference_id, size(), cats);
    }
};


#endif //SRC_GENOMEREADCLUSTER_H
