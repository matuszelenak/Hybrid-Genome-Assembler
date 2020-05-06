#include <cmath>

#include "MurmurHash3.h"
#include "CountingBloomFilter.h"


using namespace counting_bloom;


template<typename F, typename T>
void CountingBloomFilter<F, T>::add(F &item) {
    uint64_t hash_output[2];
    for (int i = 0; i < hash_count; i++) {
        MurmurHash3_x64_128(&item, sizeof(item), i * 47, &hash_output[0]);
        ++(*data)[hash_output[0] & index_mask];
    }
}

template<typename F, typename T>
T CountingBloomFilter<F, T>::get_count(F &item) {
    T count;
    count = ~(count & 0u);

    uint64_t hash_output[2];
    for (int i = 0; i < hash_count; i++) {
        MurmurHash3_x64_128(&item, sizeof(item), i * 47, &hash_output[0]);
        count = std::min(count, (*data)[hash_output[0] & index_mask]);
    }
    return count;
}

template<typename F, typename T>
uint64_t CountingBloomFilter<F, T>::get_size(uint64_t expected_items, double fp_prob) {
    return (uint64_t) (-(expected_items * log(fp_prob)) / ((log(2) * log(2))));
}

template<typename F, typename T>
int CountingBloomFilter<F, T>::get_hash_count(uint64_t expected_items, uint64_t actual_size) {
    return (int) (((double) actual_size / (double) expected_items) * log(2));
}

template<typename F, typename T>
CountingBloomFilter<F, T>::CountingBloomFilter(uint64_t expected_items, double fp_prob) {
    uint64_t minimal_size = get_size(expected_items, fp_prob);
    while (size < minimal_size) {
        size <<= 1u;
    }
    std::cout << size << std::endl;
    index_mask = size - 1;
    hash_count = get_hash_count(expected_items, size);
    data = new std::vector<T>(size, 0);
}

template<typename F, typename T>
CountingBloomFilter<F, T>::CountingBloomFilter(std::ifstream &f) {
    f.read((char*)&size, sizeof(size));
    f.read((char*)&hash_count, sizeof(hash_count));
    index_mask = size - 1;

    data = new std::vector<T>(size, 0);
    f.read((char*)data->data(), size * sizeof(T));
}

template<typename F, typename T>
void CountingBloomFilter<F, T>::dump(std::ofstream &f) {
    f.write((char*)&size, sizeof(size));
    f.write((char*)&hash_count, sizeof(hash_count));
    f.write((char*)data->data(), size * sizeof(T));
}


template class counting_bloom::CountingBloomFilter<uint64_t , uint16_t>;
