#include <vector>

#include "../common/SequenceRecordIterator.h"
#include "../lib/CountingBloomFilter.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H


typedef std::vector<counting_bloom::CountingBloomFilter<Kmer, uint16_t>* > OccurrenceFilters;

uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k);

std::pair<int, uint64_t> get_unique_k_length(SequenceRecordIterator &read_iterator);

#endif //SRC_KMERANALYSIS_H
