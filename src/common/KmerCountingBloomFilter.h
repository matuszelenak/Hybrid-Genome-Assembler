#include <map>
#include "MurmurHash3.h"
#include "Types.h"

#ifndef SRC_COUNTINGBLOOM_H
#define SRC_COUNTINGBLOOM_H


#define KMER_COUNT_MAX UINT8_MAX
#define CACHE_LINE_SIZE 512
#define MAX_BIN_COUNT UINT32_MAX


typedef uint8_t KmerCount;
typedef uint32_t BinIndex;
typedef uint8_t InnerIndex;
typedef std::map<KmerCount, uint64_t> Histogram;


//class KmerBloomFilter {
//protected:
//    uint32_t bins;
//    uint8_t *data;
//    uint64_t actual_size;
//    int hash_count;
//public:
//    explicit KmerBloomFilter(uint64_t expected_items, double fp_prob);
//    ~KmerBloomFilter(){delete [] data;};
//
//    void add(Kmer kmer);
//    bool contains(Kmer kmer);
//};


class KmerCountingBloomFilter {
private:
    uint32_t bins;
    //KmerBloomFilter* single_occurrence_bf;

public:
    explicit KmerCountingBloomFilter(uint64_t expected_items, double fp_prob);

    ~KmerCountingBloomFilter();

    KmerCount *data;
    uint64_t actual_size;
    int hash_count;

    void add(Kmer kmer);

    KmerCount get_count(Kmer kmer);

    Histogram get_histogram(int lower_coverage, int upper_coverage);
};


#endif //SRC_COUNTINGBLOOM_H
