#ifndef SRC_COUNTINGBLOOM_H
#define SRC_COUNTINGBLOOM_H

#include "KmerBloomFilter.h"


typedef uint16_t KmerCount;
typedef uint8_t InnerIndex;
typedef std::map<KmerCount, uint64_t> Histogram;

#define KMER_COUNT_MAX UINT16_MAX


class KmerCountingBloomFilter : public KmerBloomFilter {
public:
    KmerCount* data;

    CategoryID categories = 1;
    explicit KmerCountingBloomFilter(std::string &path);

    KmerCountingBloomFilter(uint64_t expected_items, int k);

    KmerCountingBloomFilter(uint64_t expected_items, int k, int categories);

    void add(Kmer kmer, CategoryID category);
    void add(Kmer kmer);

    KmerCount get_count(Kmer kmer, CategoryID category);
    KmerCount get_count(Kmer kmer);
    std::pair<KmerCount, double> get_count_with_specificity(Kmer kmer);

    void dump(std::string &path);
};


#endif //SRC_COUNTINGBLOOM_H
