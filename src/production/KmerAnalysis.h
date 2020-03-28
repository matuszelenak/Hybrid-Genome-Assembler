#include "../common/Types.h"
#include "../common/SequenceRecordIterator.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H



KmerOccurrences kmer_occurrences_bf(SequenceRecordIterator &read_iterator, int k, int genome_size, int coverage);


#endif //SRC_KMERANALYSIS_H
