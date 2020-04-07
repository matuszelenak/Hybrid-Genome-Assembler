#include <cmath>

#include "KmerAnalysis.h"
#include "../lib/HyperLogLog.hpp"
#include "KmerIterator.h"

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


void approximate_kmer_count_thread(SequenceRecordIterator &read_iterator, int k, hll::HyperLogLog &hyper){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            hyper.add((void*)&it.current_kmer, sizeof(Kmer));
        }
    }
}


uint64_t get_approximate_kmer_count(SequenceRecordIterator &read_iterator, int k){
    auto hyper = hll::HyperLogLog(10);
    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    read_iterator.reset();
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(approximate_kmer_count_thread, std::ref(read_iterator), k, std::ref(hyper));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return hyper.estimate();
}