//
// Created by whiskas on 10/02/19.
//

#include <set>
#include <map>
#include <unordered_map>
#include <vector>
#include <tuple>

#ifndef SRC_READCATEGORIZER_H
#define SRC_READCATEGORIZER_H


class ReadCategorizer {

private:
    // Variables needed for union-find
    std::vector<uint64_t > parent;
    std::vector<std::vector<uint64_t >> components;
    std::vector<std::set<uint32_t >> kmers_per_component;

    // Variables for evaluation and testing
    std::unordered_map<std::string, int> read_label_to_id;
    std::vector<int> read_labels;
    std::vector<std::tuple<int, int, int > > component_stats;    // Size, consistency %, dominant label

    uint64_t get_parent(uint64_t vertex);
    bool unite(uint64_t x, uint64_t y);
    uint32_t union_metric(std::set< uint32_t> &A, std::set< uint32_t> &B);

    void print_statistics(unsigned long limit);
    std::tuple<unsigned long , int, int > recalculate_component_stats(uint64_t component_id);
public:
    explicit ReadCategorizer(
            const std::string &read_path,
            std::set< uint64_t> &characteristic_kmers,
            int k,
            int final_categories);
};


#endif //SRC_READCATEGORIZER_H
