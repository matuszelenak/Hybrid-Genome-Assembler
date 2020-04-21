#include <cmath>
#include <fstream>
#include <iostream>
#include <fmt/format.h>
#include "KmerCountingBloomFilter.h"

uint16_t CACHE_LINES_PER_BIN = 2;
uint16_t SLOTS_IN_BIN = (CACHE_LINE_SIZE * CACHE_LINES_PER_BIN) / (sizeof(KmerCount) * 8);
uint16_t EXPECTED_FOR_BIN = 5 * CACHE_LINES_PER_BIN;
uint8_t INNER_INDEX_BITS = log2(SLOTS_IN_BIN);
uint64_t INNER_INDEX_MASK = 0xFFFFFFFFFFFFFFFF >> (uint64_t) (64 - INNER_INDEX_BITS);


void KmerCountingBloomFilter::add(Kmer kmer) {
    add(kmer, 0);
}


void KmerCountingBloomFilter::add(Kmer kmer, CategoryID category){
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);

    BinIndex bin_index = (*((BinIndex *) hash_output)) % bins;
    KmerCount *bin_starting_cell_ptr = data + category * actual_size + bin_index * SLOTS_IN_BIN;

    auto *remaining_hash_ptr = (InnerIndex *)((uint8_t *)hash_output + sizeof(BinIndex));
    for (int i = 0; i < hash_count; i++) {
        ++*(bin_starting_cell_ptr + (*(remaining_hash_ptr + i) & INNER_INDEX_MASK));
    }
    free(hash_output);
}

KmerCount KmerCountingBloomFilter::get_count(Kmer kmer) {
    KmerCount total = 0;
    for (CategoryID category_id = 0; category_id < categories; category_id++){
        total += get_count(kmer, category_id);
    }
    return total;
}

KmerCount KmerCountingBloomFilter::get_count(Kmer kmer, CategoryID category) {
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % bins;
    KmerCount *bin_starting_cell_ptr = data + category * actual_size + bin_index * SLOTS_IN_BIN;

    auto *remaining_hash_ptr = (InnerIndex *)((uint8_t *)hash_output + sizeof(BinIndex));

    KmerCount count = KMER_COUNT_MAX;
    for (int i = 0; i < hash_count; i++) {
        count = std::min(count, *(bin_starting_cell_ptr + (*(remaining_hash_ptr + i) & INNER_INDEX_MASK)));
    }
    free(hash_output);
    return count;
}


std::pair<KmerCount, double> KmerCountingBloomFilter::get_count_with_specificity(Kmer kmer){
    KmerCount prevalent_count = 0;
    KmerCount total_count = 0;
    for (int category_id = 0; category_id < categories; category_id++) {
        KmerCount count = get_count(kmer, category_id);
        prevalent_count = std::max(count, prevalent_count);
        total_count += count;
    }
    return {total_count, ((double) prevalent_count / (double) total_count) * 100};
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, int k) : KmerBloomFilter() {
    this->k = k;

    bins = ceil((double)expected_items / (double)EXPECTED_FOR_BIN);

    actual_size = bins * SLOTS_IN_BIN;

    hash_count = get_hash_count(EXPECTED_FOR_BIN, SLOTS_IN_BIN);

    data = new KmerCount[actual_size];
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, int k, int categories) : KmerCountingBloomFilter(expected_items, k){
    this->categories = categories;
    data = (KmerCount*) realloc((void*)data, actual_size * categories * sizeof(KmerCount));
    std::cout << fmt::format("Kmer counting bloom filter will take {:.2f}GB of memory\n", (double) actual_size * categories / (double) 1073741824 * sizeof(KmerCount));
}

KmerCountingBloomFilter::KmerCountingBloomFilter(std::string &path) : KmerBloomFilter(){
    auto in = std::ifstream(path, std::ios::in | std::ios::binary);

    in.seekg (0, std::basic_ifstream<char>::end);
    auto data_length = (int)in.tellg() - (sizeof(k) + sizeof(categories));
    in.seekg (0, std::basic_ifstream<char>::beg);

    std::cout << fmt::format("Kmer counting bloom filter will take {:.2f}GB of memory\n", (double) data_length / (double) 1073741824);
    data = new KmerCount[data_length / sizeof(KmerCount)];

    in.read((char*)&k, sizeof(k));
    in.read((char*)&categories, sizeof(categories));
    in.read((char*)&data[0], data_length);
    in.close();

    actual_size = (data_length / categories) / sizeof(KmerCount);

    bins = actual_size / SLOTS_IN_BIN;

    uint8_t expected_for_nested = ((double)SLOTS_IN_BIN / 64) * 8;
    hash_count = get_hash_count(expected_for_nested, SLOTS_IN_BIN);
}

void KmerCountingBloomFilter::dump(std::string &path) {
    auto out = std::ofstream(path, std::ios::out | std::ios::binary);

    out.write((char*)&k, sizeof(k));
    out.write((char*)&categories, sizeof(categories));
    out.write((char*)&data[0], actual_size * categories * sizeof(KmerCount));
    out.close();
}
