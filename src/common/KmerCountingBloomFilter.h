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
    explicit KmerBloomFilter(uint64_t expected_items, double fp_prob);

    ~KmerBloomFilter() { delete[] data; };

    void add(Kmer kmer);

    bool contains(Kmer kmer);
};


class KmerCountingBloomFilter {
private:
    uint32_t inner_bf_count;
    //KmerBloomFilter* single_occurrence_bf;

public:
    explicit KmerCountingBloomFilter(uint64_t expected_items);

    ~KmerCountingBloomFilter();

    KmerCount *data;
    uint64_t actual_size;
    int hash_count;

    void add(Kmer kmer);

    KmerCount get_count(Kmer kmer);

    uint64_t cardinality();
};

int get_hash_count(uint64_t expected_items, uint64_t actual_size);


#endif //SRC_COUNTINGBLOOM_H
