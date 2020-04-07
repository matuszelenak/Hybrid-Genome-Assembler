#include <vector>

#include "SequenceRecordIterator.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H

int get_coverage(std::vector<ReadFileMetaData> &read_meta_data, int genome_size);

int get_genome_size(std::vector<ReadFileMetaData> &read_meta_data, int coverage);

uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k);

#endif //SRC_KMERANALYSIS_H
