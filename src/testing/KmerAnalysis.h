#include <set>
#include <map>
#include <tsl/robin_map.h>

#include "../common/SequenceRecordIterator.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_TESTKMERANALYSIS_H
#define SRC_TESTKMERANALYSIS_H

KmerSpecificity get_kmer_specificity(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k);

#endif //SRC_TESTKMERANALYSIS_H
