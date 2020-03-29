#include <set>
#include <map>
#include <tsl/robin_map.h>

#include "../common/SequenceRecordIterator.h"

#ifndef SRC_TESTKMERANALYSIS_H
#define SRC_TESTKMERANALYSIS_H

typedef int UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;

struct KmerInfo {
    uint16_t in_first_category = 0;
    uint16_t in_second_category = 0;
    uint32_t sum_of_qualities = 0;

    [[nodiscard]] uint16_t total_occurrences() const {
        return this->in_first_category + this->in_second_category;
    };
    Quality avg_quality(){
        return (Quality)(this->sum_of_qualities / this->total_occurrences());
    };
};

typedef tsl::robin_map<Kmer, KmerInfo> KmerOccurrences;

void plot_kmer_specificity(std::vector<ReadFileMetaData> &meta, std::map<int, KmerSpecificity> &specificities, int max_coverage);

KmerSpecificity get_kmer_specificity(KmerOccurrences &occurrences);

KmerOccurrences filter_characteristic_kmers(KmerOccurrences &occurrences, int coverage_lower_bound, int coverage_upper_bound);

KmerOccurrences kmer_occurrences(SequenceRecordIterator &read_iterator, int k, int max_genome_size, int coverage);

std::vector<int> get_k_sizes(int max_genome_size);

#endif //SRC_TESTKMERANALYSIS_H
