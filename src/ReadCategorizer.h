//
// Created by whiskas on 10/02/19.
//

#include <set>
#include <map>
#include <unordered_map>
#include <vector>

#ifndef SRC_READCATEGORIZER_H
#define SRC_READCATEGORIZER_H


class ReadCategorizer {

private:
    std::vector<uint64_t > parent;
    std::vector<std::vector<uint64_t >> components;
    std::vector<std::set<uint32_t >> kmers_per_component;

    uint64_t get_parent(uint64_t vertex);
    bool unite(uint64_t x, uint64_t y);
    uint32_t union_metric(std::set< uint32_t> &A, std::set< uint32_t> &B);
public:
    explicit ReadCategorizer(
            const std::string &read_path,
            std::set< uint64_t> &characteristic_kmers,
            int k,
            int final_categories);
};


#endif //SRC_READCATEGORIZER_H
