#include <set>
#include <map>
#include <tsl/robin_map.h>

#include "../common/SequenceRecordIterator.h"
#include "../common/KmerCountingBloomFilter.h"

#ifndef SRC_TESTKMERANALYSIS_H
#define SRC_TESTKMERANALYSIS_H

typedef int UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;


void plot_kmer_specificity(std::vector<ReadFileMetaData> &meta, std::map<int, KmerSpecificity> &specificities, int max_coverage);

std::vector<uint8_t> get_k_sizes(int max_genome_size);

KmerSpecificity get_kmer_specificity(SequenceRecordIterator &read_iterator, std::vector<KmerCountingBloomFilter*> &filters, int k, uint32_t expected_num_of_kmers);

std::vector<KmerCountingBloomFilter*> kmer_occurrences(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers);


#endif //SRC_TESTKMERANALYSIS_H
