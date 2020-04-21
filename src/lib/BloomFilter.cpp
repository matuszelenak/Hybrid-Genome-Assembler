#include <cmath>

#include "BloomFilter.h"
#include "MurmurHash3.h"


using namespace bloom;

template<typename F>
BloomFilter<F>::BloomFilter(uint64_t expected_items, double fp_prob) {
    uint64_t minimal_size = get_size(expected_items, fp_prob);
    while (size < minimal_size) {
        size <<= 1u;
    }
    index_mask = size - 1;
    hash_count = get_hash_count(expected_items, size);
    data = new std::vector<bool>(size, false);
}

template<typename F>
BloomFilter<F>::BloomFilter(std::ifstream &f) {
    f.read((char*)&size, sizeof(size));
    f.read((char*)&hash_count, sizeof(hash_count));
    index_mask = size - 1;
    data = new std::vector<bool>(size, false);

    int size_in_bytes = ceil((double)size / 8);
    char* bytes = new char[size_in_bytes]();
    f.read(bytes, size_in_bytes);
    for (int i = 0; i < data->size(); i++){
        (*data)[i] = (bool)(bytes[i / 8] & (1 << (i % 8)));
    }
}

template<typename F>
void BloomFilter<F>::dump(std::ofstream &f) {
    f.write((char*)&size, sizeof(size));
    f.write((char*)&hash_count, sizeof(hash_count));

    int size_in_bytes = ceil((double)size / 8);
    char* bytes = new char[size_in_bytes]();
    for (int i = 0; i < data->size(); i++){
        if ((*data)[i]){
            bytes[i / 8] |= (1 << (i % 8));
        }
    }
    f.write(bytes, size_in_bytes);
}

template<typename F>
void BloomFilter<F>::add(F &item) {
    uint64_t hash_output[2];
    for (int i = 0; i < hash_count; i++) {
        MurmurHash3_x64_128(&item, sizeof(item), i * 47, &hash_output[0]);
        (*data)[hash_output[0] & index_mask] = true;
    }
}

template<typename F>
bool BloomFilter<F>::insert(F &item) {
    bool newly_inserted = false;

    uint64_t hash_output[2];
    for (int i = 0; i < hash_count; i++) {
        MurmurHash3_x64_128(&item, sizeof(item), i * 47, &hash_output[0]);
        if (!(*data)[hash_output[0] & index_mask]){
            newly_inserted = true;
            (*data)[hash_output[0] & index_mask] = true;
        }
    }

    return newly_inserted;
}

template<typename F>
bool BloomFilter<F>::contains(F &item) {
    uint64_t hash_output[2];
    bool contains = true;
    for (int i = 0; i < hash_count; i++) {
        MurmurHash3_x64_128(&item, sizeof(item), i * 47, &hash_output[0]);
        if (!(*data)[hash_output[0] & index_mask]){
            return false;
        }
    }
    return contains;
}

template<typename F>
uint64_t BloomFilter<F>::get_size(uint64_t expected_items, double fp_prob) {
    return (uint64_t) (-(expected_items * log(fp_prob)) / ((log(2) * log(2))));
}

template<typename F>
int BloomFilter<F>::get_hash_count(uint64_t expected_items, uint64_t actual_size) {
    return (int) (((double) actual_size / (double) expected_items) * log(2));
}

template<typename F>
uint64_t BloomFilter<F>::cardinality() {
    uint64_t set_bits = std::count(data->begin(), data->end(), true);

    return - (((double)data->size() / (double)hash_count) * log(1 - (double)set_bits / (double)data->size()));
}

template class bloom::BloomFilter<int>;
template class bloom::BloomFilter<uint64_t>;