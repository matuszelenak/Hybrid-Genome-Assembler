#include "../common/Types.h"
#include "../common/SequenceRecordIterator.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_PRODUCTIONKMERANALYSIS_H
#define SRC_PRODUCTIONKMERANALYSIS_H

void kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, uint32_t expected_num_of_kmers);

KmerCountingBloomFilter kmer_occurrences(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers);

#endif //SRC_PRODUCTIONKMERANALYSIS_H
