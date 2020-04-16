#include <vector>

#include "SequenceRecordIterator.h"
#include "KmerCountingBloomFilter.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H

int get_coverage(std::vector<ReadFileMetaData> &read_meta_data, int genome_size);

int get_genome_size(std::vector<ReadFileMetaData> &read_meta_data, int coverage);

std::vector<unsigned int> get_k_sizes(int max_genome_size);

uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k);

std::pair<int, uint64_t> get_unique_k_length(SequenceRecordIterator &read_iterator);

Histogram kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, uint32_t kmer_occurrence_histogram);

KmerCountingBloomFilter kmer_occurrence_filter(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers);

KmerSpecificity get_kmer_specificity(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, uint32_t expected_num_of_kmers);

void visualize_kmer_positions(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, int lower_coverage, int upper_coverage);

#endif //SRC_KMERANALYSIS_H
