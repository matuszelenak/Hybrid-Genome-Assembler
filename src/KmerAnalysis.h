#include <map>

#include "DNAStructures.h"
#include "ReadDataLoader.h"


#ifndef SRC_KMERANALYSIS_H
#define SRC_KMERANALYSIS_H

typedef int UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;

void plot_kmer_specificity(std::map<int, KmerSpecificity> &specificities, int max_coverage);
KmerSpecificity get_kmer_specificity(KmerOccurrences &occurrences);
KmerOccurrences filter_characteristic_kmers(KmerOccurrences &occurrences, int coverage_lower_bound, int coverage_upper_bound);
KmerOccurrences kmer_occurrences(ReadDataLoader &reader, int k);

#endif //SRC_KMERANALYSIS_H
