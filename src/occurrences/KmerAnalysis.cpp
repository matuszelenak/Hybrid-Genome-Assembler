#include <cmath>

#include "KmerAnalysis.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"

#include "../lib/BloomFilter.h"
#include "../lib/HyperLogLog.hpp"


namespace cbf = counting_bloom;
namespace bf = bloom;


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

    read_iterator.rewind();
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(approximate_kmer_count_thread, std::ref(read_iterator), k, std::ref(hyper));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return hyper.estimate();
}


std::pair<int, uint64_t > get_unique_k_length(SequenceRecordIterator &read_iterator){
    uint64_t previous_count = get_approximate_kmer_count(read_iterator, 11);
    int k = 13;
    std::cout << fmt::format("k=11 : ~{} kmers", previous_count) << std::endl;
    while (k < 33){
        uint64_t count = get_approximate_kmer_count(read_iterator, k);
        std::cout << fmt::format("k={} : ~{} kmers", k, count) << std::endl;
        if (count < previous_count || ((double)std::min(count, previous_count) / (double)std::max(count, previous_count)) > 0.9){
            return {k, count};
        }
        k += 2;
        previous_count = count;
    }
    return {k, previous_count};
}
