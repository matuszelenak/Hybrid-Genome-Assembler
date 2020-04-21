#ifndef SRC_BLOOMFILTER_H
#define SRC_BLOOMFILTER_H

#include <iostream>
#include <vector>
#include <fstream>

namespace bloom{

    template<typename F>
    class BloomFilter {
    public:
        std::vector<bool>* data;
        uint64_t size = 1;
        uint64_t index_mask = 1;
        uint16_t hash_count;

        explicit BloomFilter(std::ifstream &f);
        BloomFilter(uint64_t expected_items, double fp_prob);
        ~BloomFilter(){ free(data); };

        static uint64_t get_size(uint64_t expected_items, double fp_prob);
        static int get_hash_count(uint64_t expected_items, uint64_t actual_size);

        uint64_t cardinality();

        void add(F &item);
        bool contains(F &item);
        bool insert(F &item);

        void dump(std::ofstream &f);
    };
}

#endif //SRC_BLOOMFILTER_H
