#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/thread.hpp>
#include <iostream>

#include "KmerAnalysis.h"

#include "../common/KmerIterator.h"
#include "../common/KmerCountingBloomFilter.h"
#include "../common/KmerAnalysis.h"
#include "../common/Utils.h"


namespace algo = boost::algorithm;
namespace pt = boost::posix_time;


std::mutex mut;


void kmer_occurrence_histogram_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, KmerBloomFilter &processed, int k, uint64_t* counts){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            KmerCount c = bf.get_count(it.current_kmer);
            if (c < 3) continue;

            // Who needs thread safety anyways? it's not like two identical kmers can be processed simultaneously...right?
            if (!processed.contains(it.current_kmer)){
                processed.add(it.current_kmer);
                ++*(counts + c);
            }
        }
    }
}


void kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, uint32_t expected_num_of_kmers){
    Histogram hist;
    KmerBloomFilter processed = KmerBloomFilter(expected_num_of_kmers, 0.01);

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    read_iterator.reset();
    pt::ptime ts_start = pt::microsec_clock::local_time();
    auto* counts = new uint64_t[KMER_COUNT_MAX + 1];
    for (int i = 0; i <= KMER_COUNT_MAX; i++) *(counts + i) = 0;

    for (auto &i : t) i = std::thread(kmer_occurrence_histogram_thread, std::ref(read_iterator), std::ref(bf), std::ref(processed), k, counts);
    for (auto &i : t) i.join();
    std::cout << fmt::format("Histogram calculated in in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    uint64_t total = 0;
    for (int i = 1; i <= KMER_COUNT_MAX; i++){
        if (counts[i] < 50) continue;
        hist[i] = counts[i];
        total += counts[i];
    }
    delete [] counts;


    std::vector<std::string> occurrence_strings;
    std::transform(hist.begin(), hist.end(), std::back_inserter(occurrence_strings), [](Histogram::value_type occ) -> std::string {
        return fmt::format("{}: {}", occ.first, occ.second);
    });
    std::string plot_input = fmt::format("{{{}}}", algo::join(occurrence_strings, ", "));
    run_command_with_input("python common/python/plot.py", plot_input);
}


void kmer_occurrences_first_pass(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            bf.add(it.current_kmer);
        }
    }
}


KmerCountingBloomFilter kmer_occurrences(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers){
    auto bf = KmerCountingBloomFilter(expected_num_of_kmers);

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    // First pass to fill the bloom filter
    read_iterator.reset();
    pt::ptime ts_start = pt::microsec_clock::local_time();
    for (auto &i : t) i = std::thread(kmer_occurrences_first_pass, std::ref(read_iterator), std::ref(bf), k);
    for (auto &i : t) i.join();
    std::cout << fmt::format("First pass done in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    return bf;
}