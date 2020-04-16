#include <map>
#include "../lib/MurmurHash3.h"
#include "Types.h"

#ifndef SRC_COUNTINGBLOOM_H
#define SRC_COUNTINGBLOOM_H


typedef uint16_t KmerCount;
typedef uint32_t BinIndex;
typedef uint8_t InnerIndex;
typedef std::map<KmerCount, uint64_t> Histogram;

#define KMER_COUNT_MAX UINT16_MAX
#define CACHE_LINE_SIZE 512
#define MAX_BIN_COUNT UINT32_MAX


class KmerBloomFilter {
protected:
    uint32_t bins;
    uint8_t *data;
    uint64_t actual_size;
    int hash_count;
public:
    explicit KmerBloomFilter(std::string &path);

    explicit KmerBloomFilter(uint64_t expected_items);

    ~KmerBloomFilter() { delete[] data; };

    void add(Kmer kmer);

    bool contains(Kmer kmer);

    void dump(std::string &path);
};


class KmerCountingBloomFilter {
public:
    uint32_t inner_bf_count;
    CategoryID categories = 1;
    explicit KmerCountingBloomFilter(std::string &path);

    KmerCountingBloomFilter(uint64_t expected_items, int k);

    KmerCountingBloomFilter(uint64_t expected_items, int k, int categories);

    ~KmerCountingBloomFilter();

    KmerCount *data;
    uint64_t actual_size;
    int hash_count;
    int k;

    void add(Kmer kmer, CategoryID category);
    void add(Kmer kmer);

    KmerCount get_count(Kmer kmer, CategoryID category);
    KmerCount get_count(Kmer kmer);

    void dump(std::string &path);
};

int get_hash_count(uint64_t expected_items, uint64_t actual_size);


#endif //SRC_COUNTINGBLOOM_H
