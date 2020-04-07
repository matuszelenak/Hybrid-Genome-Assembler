#include <cmath>
#include <iostream>
#include <fmt/format.h>

#include "KmerCountingBloomFilter.h"


uint16_t INNER_BF_SIZE = CACHE_LINE_SIZE / (sizeof(KmerCount) * 8) * 2;
uint8_t INNER_INDEX_BITS = log2(INNER_BF_SIZE);
uint64_t INNER_INDEX_MASK = 0xFFFFFFFFFFFFFFFF >> (uint8_t) (64 - INNER_INDEX_BITS);


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
    KmerCount *start_cell = data + (*((BinIndex *) hash_output) % inner_bf_count) * INNER_BF_SIZE;
    InnerIndex *remaining_hash_ptr = (InnerIndex *)hash_output + sizeof(BinIndex);
    for (int i = 0; i < hash_count; i++) {
        ++*(start_cell + (*(remaining_hash_ptr + i) & INNER_INDEX_MASK));
    }
    free(hash_output);
}

KmerCount KmerCountingBloomFilter::get_count(Kmer kmer) {
    KmerCount count = KMER_COUNT_MAX;
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % inner_bf_count;
    for (int i = 0; i < hash_count; i++) {
        InnerIndex inner_index = *((InnerIndex *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[bin_index * INNER_BF_SIZE + inner_index] < count) {
            count = data[bin_index * INNER_BF_SIZE + inner_index];
        }
    }
    free(hash_output);
    return count;
}

uint64_t KmerCountingBloomFilter::cardinality() {
    double total = 0;
    for (uint32_t bin_index = 0; bin_index < inner_bf_count; bin_index++){
        uint16_t non_zero_cells = 0;
        for (InnerIndex i = 0; i < INNER_BF_SIZE; i++){
            non_zero_cells += (data[bin_index * INNER_BF_SIZE + i] != 0);
        }
        total += -((double)INNER_BF_SIZE / (double)hash_count) * log(1 - ((double)non_zero_cells) / (double)INNER_BF_SIZE);
    }
    return total;
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items) {
    uint8_t expected_for_nested = ((double)INNER_BF_SIZE / 64) * 8;

    inner_bf_count = ceil((double)expected_items / (double)expected_for_nested);

    actual_size = inner_bf_count * INNER_BF_SIZE;

    hash_count = get_hash_count(expected_for_nested, INNER_BF_SIZE);

    std::cout << fmt::format("Kmer counting bloom filter will take {:.2f}GB of memory\n", (double) actual_size * (double) sizeof(KmerCount) / (double) 1073741824);

    data = new KmerCount[actual_size];

    //single_occurrence_bf = new KmerBloomFilter(expected_items, fp_prob);
}

KmerCountingBloomFilter::~KmerCountingBloomFilter() {
    delete [] data;
    //delete single_occurrence_bf;
}

KmerBloomFilter::KmerBloomFilter(uint64_t expected_items, double fp_prob) {
    actual_size = get_size(expected_items, fp_prob);
    hash_count = get_hash_count(expected_items, actual_size);
    actual_size = (uint64_t)(actual_size / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
    bins = actual_size / CACHE_LINE_SIZE;
    std::cout << fmt::format("Kmer bloom filter will take {:.2f}GB of memory\n", (double) actual_size / 8  / (double) 1073741824);
    data = new uint8_t[actual_size / 8];
}

void KmerBloomFilter::add(Kmer kmer){
    auto hash_output = (uint8_t*)malloc(32);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 1, hash_output + 16);

    uint8_t *start_cell = data + (*((BinIndex *) hash_output) % bins) * 64;
    auto remaining_hash_ptr = (uint16_t *)(hash_output + sizeof(BinIndex));
    for (int i = 0; i < hash_count; i++) {
        *(start_cell + (*(remaining_hash_ptr + i) & 0b111111u)) |= (1u << ((*(remaining_hash_ptr + i) & 0b111000000u) >> 6u));
    }
    free(hash_output);
}

bool KmerBloomFilter::contains(Kmer kmer){
    bool contains = true;

    auto hash_output = (uint8_t*)malloc(32);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 1, hash_output + 16);

    uint8_t *start_cell = data + (*((BinIndex *) hash_output) % bins) * 64;
    auto remaining_hash_ptr = (uint16_t *)(hash_output + sizeof(BinIndex));
    for (int i = 0; i < hash_count; i++) {
        if(!(*(start_cell + (*(remaining_hash_ptr + i) & 0b111111u)) & (1u << ((*(remaining_hash_ptr + i) & 0b111000000u) >> 6u)))){
            contains = false;
            break;
        }
    }
    free(hash_output);
    return contains;
}
