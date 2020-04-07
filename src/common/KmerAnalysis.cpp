#include <cmath>

#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "KmerCountingBloomFilter.h"
#include "Utils.h"

#include "../lib/HyperLogLog.hpp"

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


std::vector<uint8_t> get_k_sizes(int max_genome_size) {
    std::vector<uint8_t> k_sizes;
    uint8_t k_guess = (uint8_t) ceil(log(max_genome_size) / log(4)) | 1u;
    for (int k_value = k_guess; k_value <= 19; k_value += 2) {
        k_sizes.push_back(k_value);
    }
    return k_sizes;
}


void kmer_occurrence_histogram_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, KmerBloomFilter &processed, int k, uint64_t* counts){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            // Who needs thread safety anyways? it's not like two identical kmers can be processed simultaneously...right?
            if (!processed.contains(it.current_kmer)){
                processed.add(it.current_kmer);
                ++*(counts + filter.get_count(it.current_kmer));
            }
        }
    }
}


Histogram kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k){
    KmerBloomFilter processed = KmerBloomFilter(filter._expected_items, 0.001);

    read_iterator.reset();
    auto* counts = new uint64_t[KMER_COUNT_MAX + 1]();
    auto runner = ThreadRunner(kmer_occurrence_histogram_thread, std::ref(read_iterator), std::ref(filter), std::ref(processed), k, counts);

    Histogram hist;
    for (int i = 5; i <= KMER_COUNT_MAX; i++){
        if (counts[i] < 50) continue;
        hist[i] = counts[i];
    }
    delete [] counts;
    return hist;
}


void kmer_occurrences_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k) {
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            filter.add(it.current_kmer, read->category_id);
        }
    }
}


KmerCountingBloomFilter kmer_occurrence_filter(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers) {
    KmerCountingBloomFilter filter(expected_num_of_kmers, read_iterator.categories);
    read_iterator.reset();
    auto runner = ThreadRunner(kmer_occurrences_thread, std::ref(read_iterator), std::ref(filter), k);
    return filter;
}