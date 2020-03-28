#include <cmath>
#include <iostream>
#include <fmt/format.h>

#include "KmerCountingBloomFilter.h"
#include "smhasher/MurmurHash3.h"


uint16_t BIN_SIZE = CACHE_LINE_SIZE / (sizeof(KmerCount) * 8);
uint8_t BITS_FOR_INNER_INDEX = log2(BIN_SIZE);
uint64_t INNER_INDEX_MASK = 0xFFFFFFFFFFFFFFFF >> (uint8_t) (64 - BITS_FOR_INNER_INDEX);


int get_hash_count(uint64_t expected_items, uint64_t actual_size) {
    return (int) (((double) actual_size / (double) expected_items) * log(2));
}

uint64_t get_size(uint64_t expected_items, double fp_prob) {
    return (uint64_t) (-(expected_items * log(fp_prob)) / ((log(2) * log(2))));
}

// This is actually thread safe LMAO
void KmerCountingBloomFilter::add(Kmer kmer) {
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    KmerCount *start_cell = data + (*((BinIndex *) hash_output) % bins) * BIN_SIZE;
    uint8_t *remaining_hash_ptr = (uint8_t *) hash_output + sizeof(BinIndex);
    for (int i = 0; i < hash_count; i++) {
        ++*(start_cell + (*(remaining_hash_ptr + i) & INNER_INDEX_MASK));
    }
    free(hash_output);
}

bool KmerCountingBloomFilter::contains(Kmer kmer) {
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % bins;
    bool contains = true;
    for (int i = 0; i < hash_count; i++) {
        uint8_t inner_index = *((uint8_t *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[bin_index * BIN_SIZE + inner_index] == 0) {
            contains = false;
            break;
        }
    }
    free(hash_output);
    return contains;
}

KmerCount KmerCountingBloomFilter::get_count(Kmer kmer) {
    KmerCount count = KMER_COUNT_MAX;
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % bins;
    for (int i = 0; i < hash_count; i++) {
        uint8_t inner_index = *((uint8_t *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[bin_index * BIN_SIZE + inner_index] < count) {
            count = data[bin_index * BIN_SIZE + inner_index];
        }
    }
    free(hash_output);
    return count;
}

bool KmerCountingBloomFilter::has_kmer_count_in_range(Kmer kmer, int lower, int upper){
    bool result = true;
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % bins;
    for (int i = 0; i < hash_count; i++) {
        uint8_t inner_index = *((uint8_t *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[bin_index * BIN_SIZE + inner_index] < lower || data[bin_index * BIN_SIZE + inner_index] > upper){
            result = false;
            break;
        }
    }
    free(hash_output);
    return result;
}

// TODO figure out a better approximation
uint64_t KmerCountingBloomFilter::cardinality() {
    return get_kmer_count_in_count_range(1, KMER_COUNT_MAX);
}


uint64_t KmerCountingBloomFilter::get_kmer_count_in_count_range(KmerCount lower_bound, KmerCount upper_bound) {
    uint64_t cells = 0;
    for (int i = 0; i < actual_size; i++) {
        cells += (data[i] >= lower_bound && data[i] <= upper_bound);
    }
    return (uint64_t) ((-(double) actual_size / (double) hash_count) * log(1 - ((double) cells / (double) actual_size)));
}

KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, double fp_prob) {
    int bits_for_bin_index = sizeof(BinIndex) * 8;

    actual_size = get_size(expected_items, fp_prob);
    actual_size = (uint64_t) (ceil((double) actual_size / (double) BIN_SIZE)) * BIN_SIZE; // Round up to the multiple of bin size;

    hash_count = std::min(get_hash_count(expected_items, actual_size), (128 - bits_for_bin_index) / 8);

    bins = actual_size / BIN_SIZE;
    std::cout << fmt::format("Kmer bloom filter will take {:.2f}GB of memory\n", (double) actual_size * (double) sizeof(KmerCount) / (double) 1073741824);

    data = new KmerCount[actual_size];
}

KmerCountingBloomFilter::~KmerCountingBloomFilter() {
    free(data);
}