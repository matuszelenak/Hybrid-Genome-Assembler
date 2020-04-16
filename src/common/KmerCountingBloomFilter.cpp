#include <cmath>
#include <iostream>
#include <fmt/format.h>
#include <fstream>

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

void KmerCountingBloomFilter::add(Kmer kmer) {
    add(kmer, 0);
}


void KmerCountingBloomFilter::add(Kmer kmer, CategoryID category){
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    KmerCount *start_cell = data + (*((BinIndex *) hash_output) % inner_bf_count) * INNER_BF_SIZE + category * actual_size;
    InnerIndex *remaining_hash_ptr = (InnerIndex *)hash_output + sizeof(BinIndex);
    for (int i = 0; i < hash_count; i++) {
        ++*(start_cell + (*(remaining_hash_ptr + i) & INNER_INDEX_MASK));
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
    KmerCount count = KMER_COUNT_MAX;
    void *hash_output = malloc(16);
    MurmurHash3_x64_128(&kmer, sizeof(Kmer), 0, hash_output);
    BinIndex bin_index = *((BinIndex *) hash_output) % inner_bf_count;
    for (int i = 0; i < hash_count; i++) {
        InnerIndex inner_index = *((InnerIndex *) hash_output + i + sizeof(BinIndex)) & INNER_INDEX_MASK;
        if (data[category * actual_size + bin_index * INNER_BF_SIZE + inner_index] < count) {
            count = data[category * actual_size + bin_index * INNER_BF_SIZE + inner_index];
        }
    }
    free(hash_output);
    return count;
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, int k) {
    this->k = k;
    uint8_t expected_for_nested = ((double)INNER_BF_SIZE / 64) * 8;

    inner_bf_count = ceil((double)expected_items / (double)expected_for_nested);

    actual_size = inner_bf_count * INNER_BF_SIZE;

    hash_count = get_hash_count(expected_for_nested, INNER_BF_SIZE);

    data = new KmerCount[actual_size];
}


KmerCountingBloomFilter::KmerCountingBloomFilter(uint64_t expected_items, int k, int categories) : KmerCountingBloomFilter(expected_items, k){
    this->categories = categories;
    delete [] data;
    data = new KmerCount[actual_size * categories];
}

KmerCountingBloomFilter::KmerCountingBloomFilter(std::string &path){
    auto in = std::ifstream(path, std::ios::in | std::ios::binary);

    in.seekg (0, std::basic_ifstream<char>::end);
    auto data_length = (int)in.tellg() - (sizeof(k) + sizeof(categories));
    in.seekg (0, std::basic_ifstream<char>::beg);

    data = new KmerCount[data_length / sizeof(KmerCount)];

    in.read((char*)&k, sizeof(k));
    in.read((char*)&categories, sizeof(categories));
    in.read((char*)&data[0], data_length);
    in.close();

    actual_size = (data_length / categories) / sizeof(KmerCount);

    inner_bf_count = actual_size / INNER_BF_SIZE;

    uint8_t expected_for_nested = ((double)INNER_BF_SIZE / 64) * 8;
    hash_count = get_hash_count(expected_for_nested, INNER_BF_SIZE);
}

void KmerCountingBloomFilter::dump(std::string &path) {
    auto out = std::ofstream(path, std::ios::out | std::ios::binary);

    out.write((char*)&k, sizeof(k));
    out.write((char*)&categories, sizeof(categories));
    out.write((char*)&data[0], actual_size * categories * sizeof(KmerCount));
    out.close();
}

KmerCountingBloomFilter::~KmerCountingBloomFilter() {
    delete [] data;
}



KmerBloomFilter::KmerBloomFilter(uint64_t expected_items) {
    bins = ceil((double)expected_items / (double)54); // We can fit ~54 items into 512 bits of memory with fp = 0.01
    actual_size = bins * CACHE_LINE_SIZE;
    hash_count = get_hash_count(54, CACHE_LINE_SIZE);

    std::cout << fmt::format("Kmer bloom filter will take {:.2f}GB of memory\n", (double) actual_size / 8  / (double) 1073741824);
    data = new uint8_t[actual_size / 8];
}

KmerBloomFilter::KmerBloomFilter(std::string &path){
    auto in = std::ifstream(path, std::ios::in | std::ios::binary);

    in.seekg (0, std::basic_ifstream<char>::end);
    auto data_length = in.tellg();
    in.seekg (0, std::basic_ifstream<char>::beg);

    actual_size = data_length * 8;
    bins = actual_size / CACHE_LINE_SIZE;
    hash_count = get_hash_count(54, CACHE_LINE_SIZE);
    data = new uint8_t[data_length];

    in.read((char*)&data[0], data_length);
    in.close();
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

void KmerBloomFilter::dump(std::string &path) {
    auto out = std::ofstream(path, std::ios::out | std::ios::binary);
    out.write((char*)&data[0], actual_size / 8);
    out.close();
}