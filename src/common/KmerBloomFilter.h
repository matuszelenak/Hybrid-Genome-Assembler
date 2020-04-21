#ifndef SRC_KMERBLOOMFILTER_H
#define SRC_KMERBLOOMFILTER_H

#include "../lib/MurmurHash3.h"
#include "Types.h"

typedef uint32_t BinIndex;

#define CACHE_LINE_SIZE 512

class KmerBloomFilter {
protected:
    uint32_t bins;
    uint8_t *data;
    uint64_t actual_size;
    int hash_count;

public:
    KmerBloomFilter() = default;

    explicit KmerBloomFilter(std::string &path);

    explicit KmerBloomFilter(uint64_t expected_items, int k);

    void add(Kmer kmer);

    bool contains(Kmer kmer);

    void dump(std::string &path);

    static int get_hash_count(uint64_t _expected_items, uint64_t _actual_size);

    int k;
};


#endif //SRC_KMERBLOOMFILTER_H
