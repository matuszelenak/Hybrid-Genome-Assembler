#include "../common/Types.h"
#include "../common/SequenceRecordIterator.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_PRODUCTIONKMERANALYSIS_H
#define SRC_PRODUCTIONKMERANALYSIS_H

KmerCountingBloomFilter kmer_occurrences(SequenceRecordIterator &read_iterator, int k, int genome_size, int coverage);

#endif //SRC_PRODUCTIONKMERANALYSIS_H
