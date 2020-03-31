#include <cmath>
#include <iostream>
#include <fmt/format.h>

#include "KmerCountingBloomFilter.h"


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
    BinIndex bin_index = *((BinIndex *) hash_output) % bins;
    for (int i = 0; i < hash_count; i++) {
        InnerIndex inner_index = *((InnerIndex *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[bin_index * BIN_SIZE + inner_index] < count) {
            count = data[bin_index * BIN_SIZE + inner_index];
        }
    }
    free(hash_output);
    return count;
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, double fp_prob) {
    int bits_for_bin_index = sizeof(BinIndex) * 8;

    actual_size = get_size(expected_items, fp_prob);
    actual_size = (uint64_t) (ceil((double) actual_size / (double) BIN_SIZE)) * BIN_SIZE; // Round up to the multiple of bin size;

    hash_count = std::min(get_hash_count(expected_items, actual_size), (128 - bits_for_bin_index) / 8);

    bins = actual_size / BIN_SIZE;
    std::cout << fmt::format("Kmer bloom filter will take {:.2f}GB of memory\n", (double) actual_size * (double) sizeof(KmerCount) / (double) 1073741824);

    data = new KmerCount[actual_size];

    //single_occurrence_bf = new KmerBloomFilter(expected_items, fp_prob);
}

KmerCountingBloomFilter::~KmerCountingBloomFilter() {
    delete [] data;
    //delete single_occurrence_bf;
}

Histogram KmerCountingBloomFilter::get_histogram(int lower_coverage, int upper_coverage) {
    Histogram hist;
    for (int i = 0; i < actual_size; i++) {
        if (data[i] < lower_coverage || data[i] > upper_coverage) continue;
        hist.insert(Histogram::value_type(data[i], 0)).first->second++;
    }
    for (auto it = begin(hist); it != end(hist);){
        // TODO this formula seems to consistently under-represent, FIX
        it->second = (uint64_t) ((-(double) actual_size / (double) hash_count) * log(1 - ((double) it->second / (double) actual_size)));
        if (it->second == 0){
            it = hist.erase(it);
        } else ++it;
    }
    return hist;
}

//KmerBloomFilter::KmerBloomFilter(uint64_t expected_items, double fp_prob) {
//    actual_size = get_size(expected_items, fp_prob);
//    hash_count = get_hash_count(expected_items, actual_size);
//    actual_size = (int)(actual_size / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
//    bins = actual_size / CACHE_LINE_SIZE;
//    data = new uint8_t[actual_size / 8];
//}
//
//void KmerBloomFilter::add(Kmer kmer){
//    auto hash_output = (uint8_t*)malloc(32);
//    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
//    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 1, hash_output + 16);
//
//    KmerCount *start_cell = data + (*((BinIndex *) hash_output) % bins) * BIN_SIZE;
//    auto remaining_hash_ptr = (uint16_t *)(hash_output + sizeof(BinIndex));
//    for (int i = 0; i < hash_count; i++) {
//        *(start_cell + (*(remaining_hash_ptr + i) & 0b111111u)) |= (1u << ((*(remaining_hash_ptr + i) & 0b111000000u) >> 6u));
//    }
//    free(hash_output);
//}
//
//bool KmerBloomFilter::contains(Kmer kmer){
//    bool contains = true;
//
//    auto hash_output = (uint8_t*)malloc(32);
//    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
//    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 1, hash_output + 16);
//
//    KmerCount *start_cell = data + (*((BinIndex *) hash_output) % bins) * BIN_SIZE;
//    auto remaining_hash_ptr = (uint16_t *)(hash_output + sizeof(BinIndex));
//    for (int i = 0; i < hash_count; i++) {
//        if(!(*(start_cell + (*(remaining_hash_ptr + i) & 0b111111u)) & (1u << ((*(remaining_hash_ptr + i) & 0b111000000u) >> 6u)))){
//            contains = false;
//            break;
//        }
//    }
//    free(hash_output);
//    return contains;
//}
