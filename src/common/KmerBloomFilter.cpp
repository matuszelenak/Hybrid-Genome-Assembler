#include <fstream>
#include <cmath>
#include <iostream>
#include <fmt/format.h>
#include "KmerBloomFilter.h"


KmerBloomFilter::KmerBloomFilter(uint64_t expected_items, int k): k(k) {
    bins = ceil((double)expected_items / (double)54); // We can fit ~54 items into 512 bits of memory with fp = 0.01
    actual_size = (uint64_t)bins * (uint64_t)CACHE_LINE_SIZE;
    hash_count = get_hash_count(54, CACHE_LINE_SIZE);

    std::cout << fmt::format("Kmer bloom kmers will take {:.2f}GB of memory\n", (double) actual_size / 8  / (double) 1073741824);
    data = new uint8_t[actual_size / 8]();
}

KmerBloomFilter::KmerBloomFilter(std::string &path): k(0){
    auto in = std::ifstream(path, std::ios::in | std::ios::binary);

    in.seekg (0, std::basic_ifstream<char>::end);
    auto data_length = (int)in.tellg() - sizeof(k);
    in.seekg (0, std::basic_ifstream<char>::beg);

    actual_size = data_length * 8;
    bins = actual_size / CACHE_LINE_SIZE;
    hash_count = get_hash_count(54, CACHE_LINE_SIZE);
    data = new uint8_t[data_length];

    in.read((char*)&k, sizeof(k));
    in.read((char*)&data[0], data_length);
    in.close();
}


int KmerBloomFilter::get_hash_count(uint64_t _expected_items, uint64_t _actual_size) {
    return (int) (((double) _actual_size / (double) _expected_items) * log(2));
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
    out.write((char*)&k, sizeof(k));
    out.write((char*)&data[0], actual_size / 8);
    out.close();
}