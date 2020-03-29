#include <vector>

#include "SequenceRecordIterator.h"

#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H

int get_coverage(std::vector<ReadFileMetaData> &read_meta_data, int genome_size);

int get_genome_size(std::vector<ReadFileMetaData> &read_meta_data, int coverage);

uint32_t get_num_of_expected_kmers(uint k, uint genome_size, uint coverage, uint read_length, double error_rate);

#endif //SRC_KMERANALYSIS_H
