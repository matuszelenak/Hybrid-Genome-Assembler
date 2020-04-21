#ifndef SRC_COUNTINGBLOOMFILTER_H
#define SRC_COUNTINGBLOOMFILTER_H

#include <iostream>
#include <vector>
#include <fstream>

namespace counting_bloom{

    template<typename F, typename T>
    class CountingBloomFilter {
    public:
        std::vector<T>* data;
        uint64_t size = 1;
        uint64_t index_mask = 1;
        uint16_t hash_count = 1;

        explicit CountingBloomFilter(std::ifstream &f);
        CountingBloomFilter(uint64_t expected_items, double fp_prob);
        ~CountingBloomFilter(){ free(data); };

        static uint64_t get_size(uint64_t expected_items, double fp_prob);
        static int get_hash_count(uint64_t expected_items, uint64_t actual_size);

        void add(F &item);
        T get_count(F &item);

        void dump(std::ofstream &f);
    };
}


#endif //SRC_COUNTINGBLOOMFILTER_H
