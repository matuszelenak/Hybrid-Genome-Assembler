#ifndef SRC_KMEROCCURRENCECOUNTER_H
#define SRC_KMEROCCURRENCECOUNTER_H


#include "../common/SequenceRecordIterator.h"
#include "../lib/CountingBloomFilter.h"
#include "../lib/BloomFilter.h"

typedef std::vector<counting_bloom::CountingBloomFilter<Kmer, uint16_t>* > OccurrenceFilters;

class KmerOccurrenceCounter {
    SequenceRecordIterator* reader = nullptr;
    OccurrenceFilters filters;

public:
    KmerOccurrenceCounter() = default;
    uint32_t expected_num_of_kmers = 0;
    int k = 0;

    KmerSpecificity specificity;
    Histogram histogram;

    explicit KmerOccurrenceCounter(SequenceRecordIterator &read_iterator);
    KmerOccurrenceCounter(SequenceRecordIterator &read_iterator, int k);

    void compute_filters();

    KmerSpecificity get_specificity(std::vector<double> &thresholds);

    Histogram get_histogram();

    void export_kmers_in_range(int lower_bound, int upper_bound, std::string &path);

    bloom::BloomFilter<Kmer> count_non_singletons();
};


#endif //SRC_KMEROCCURRENCECOUNTER_H
