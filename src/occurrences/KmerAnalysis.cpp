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
    int k = 11;
    long long previous_count = get_approximate_kmer_count(read_iterator, k);
    std::cout << fmt::format("k=11 : ~{} kmers", previous_count) << std::endl;
    while (k < 33){
        long long count = get_approximate_kmer_count(read_iterator, k + 2);
        std::cout << fmt::format("k={} : ~{} kmers", k + 2, count) << std::endl;
        if (((double)abs(count - previous_count) / ((double)(count + previous_count) / 2.0)) < 0.1){
            return {k, previous_count};
        }

        k += 2;
        previous_count = count;
    }
    return {k, previous_count};
}
