#ifndef SRC_KMEROCCURRENCECOUNTER_H
#define SRC_KMEROCCURRENCECOUNTER_H


#include "SequenceRecordIterator.h"
#include "../lib/CountingBloomFilter.h"
#include "../lib/BloomFilter.h"

typedef std::vector<counting_bloom::CountingBloomFilter<Kmer, uint16_t>* > OccurrenceFilters;
typedef std::map<uint16_t, uint64_t> Histogram;

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

    void compute_filters_thread();
    void compute_filters();

    void specificity_thread(bloom::BloomFilter<Kmer> &processed, std::mutex &mut);
    KmerSpecificity get_specificity(std::vector<double> &thresholds);

    void get_histogram_thread(uint64_t *counts, bloom::BloomFilter<Kmer> &processed, std::mutex &mut);
    Histogram get_histogram();

    void export_kmers_thread(bloom::BloomFilter<Kmer> &exported_kmers, bloom::BloomFilter<Kmer> &processed, int lower_bound, int upper_bound);
    void export_kmers_in_range(int lower_bound, int upper_bound, std::string &path);
};


#endif //SRC_KMEROCCURRENCECOUNTER_H
