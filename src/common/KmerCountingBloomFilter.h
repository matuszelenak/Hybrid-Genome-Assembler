#include "MurmurHash3.h"
#include "Types.h"

#ifndef SRC_COUNTINGBLOOM_H
#define SRC_COUNTINGBLOOM_H


#define KMER_COUNT_MAX UINT8_MAX
#define CACHE_LINE_SIZE 512
#define MAX_BIN_COUNT UINT32_MAX


typedef uint8_t KmerCount;
typedef uint32_t BinIndex;


class KmerCountingBloomFilter {
private:
    uint32_t bins;
public:
    explicit KmerCountingBloomFilter(uint64_t expected_items, double fp_prob);

    ~KmerCountingBloomFilter();

    KmerCount *data;
    uint64_t actual_size;
    int hash_count;

    void add(Kmer kmer);

    bool contains(Kmer kmer);

    KmerCount get_count(Kmer kmer);

    bool has_kmer_count_in_range(Kmer kmer, int lower, int upper);

    uint64_t cardinality();

    uint64_t get_kmer_count_in_count_range(KmerCount lower_bound, KmerCount upper_bound);
};


#endif //SRC_COUNTINGBLOOM_H
