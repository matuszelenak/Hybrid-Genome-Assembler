#include <vector>

#include "SequenceRecordIterator.h"
#include "KmerCountingBloomFilter.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H

int get_coverage(std::vector<ReadFileMetaData> &read_meta_data, int genome_size);

int get_genome_size(std::vector<ReadFileMetaData> &read_meta_data, int coverage);

std::vector<uint8_t> get_k_sizes(int max_genome_size);

uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k);

Histogram kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k);

KmerCountingBloomFilter kmer_occurrence_filter(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers);

#endif //SRC_KMERANALYSIS_H
