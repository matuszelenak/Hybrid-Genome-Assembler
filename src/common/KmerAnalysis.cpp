#include <cmath>

#include "KmerAnalysis.h"

int _get_genome_size_or_coverage(std::vector<ReadFileMetaData> &read_meta_data, int known) {
    uint64_t unknown = 0;
    for (auto &meta : read_meta_data) {
        unknown = std::max(unknown, (uint64_t) (meta.total_bases / known));
    }
    return unknown;
}

// For the case when only the coverage is known
int get_genome_size(std::vector<ReadFileMetaData> &read_meta_data, int coverage) {
    return _get_genome_size_or_coverage(read_meta_data, coverage);
}

// For the case when only the genome size is known
int get_coverage(std::vector<ReadFileMetaData> &read_meta_data, int genome_size) {
    return _get_genome_size_or_coverage(read_meta_data, genome_size);
}

uint32_t get_num_of_expected_kmers(uint k, uint genome_size, uint coverage, uint read_length, double error_rate) {
    return genome_size * coverage * ((double) (read_length - k + 1) / (double) read_length) * (1 - pow(1 - error_rate, k));
}
