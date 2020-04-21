#include <vector>

#include "SequenceRecordIterator.h"
#include "KmerCountingBloomFilter.h"

#include "../lib/CountingBloomFilter.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H


typedef std::vector<counting_bloom::CountingBloomFilter<Kmer, uint16_t>* > OccurrenceFilters;

uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k);

std::pair<int, uint64_t> get_unique_k_length(SequenceRecordIterator &read_iterator);

Histogram kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, uint32_t kmer_occurrence_histogram);

void visualize_kmer_positions_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int lower_coverage, int upper_coverage, std::ofstream &output);

void visualize_kmer_positions(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, int lower_coverage, int upper_coverage);

#endif //SRC_KMERANALYSIS_H
