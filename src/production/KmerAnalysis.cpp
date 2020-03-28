#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include "KmerAnalysis.h"
#include "../common/KmerIterator.h"
#include "../common/KmerCountingBloomFilter.h"


namespace algo = boost::algorithm;
namespace pt = boost::posix_time;


std::mutex mut;


void kmer_occurrences_first_pass(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            bf.add(it.current_kmer);
        }
    }
}


void kmer_occurrences_second_pass(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, int lower_bound, int upper_bound, KmerOccurrences &occ){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            if (bf.has_kmer_count_in_range(it.current_kmer, lower_bound, upper_bound)){
                mut.lock();
                occ.insert(KmerOccurrences::value_type(it.current_kmer, {bf.get_count(it.current_kmer), 0, 0}));
                mut.unlock();
            }
        }
    }
}


KmerOccurrences kmer_occurrences_bf(SequenceRecordIterator &read_iterator, int k, int genome_size, int coverage){
    uint32_t expected_kmers = get_num_of_expected_kmers(k, genome_size, coverage, read_iterator.average_read_length(), 0.01);
    auto bf = KmerCountingBloomFilter(expected_kmers, 0.001);

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    // First pass to fill the bloom filter
    read_iterator.reset();
    pt::ptime ts_start = pt::microsec_clock::local_time();
    for (auto &i : t) i = std::thread(kmer_occurrences_first_pass, std::ref(read_iterator), std::ref(bf), k);
    for (auto &i : t) i.join();
    std::cout << fmt::format("First pass done in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    ts_start = pt::microsec_clock::local_time();
    KmerOccurrences occurrences;
    int lower = 1, upper = 255;
    uint64_t expected_kmers_in_range = bf.get_kmer_count_in_count_range(lower, upper);
    std::cout << fmt::format("Expecting {} kmers in range <{}, {}>\n", expected_kmers_in_range, lower, upper);
    occurrences.rehash(expected_kmers_in_range);
    std::cout << fmt::format("Occurence map constructed in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    read_iterator.reset();
    ts_start = pt::microsec_clock::local_time();
    for (auto &i : t) i = std::thread(kmer_occurrences_second_pass, std::ref(read_iterator), std::ref(bf), k, lower, upper, std::ref(occurrences));
    for (auto &i : t) i.join();
    std::cout << fmt::format("Second pass done in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    return occurrences;
}