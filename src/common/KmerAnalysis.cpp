#include <cmath>
#include <numeric>

#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "KmerCountingBloomFilter.h"
#include "Utils.h"

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


Histogram kmer_occurrence_histogram(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, uint32_t expected_num_of_kmers){
    KmerBloomFilter processed = KmerBloomFilter(expected_num_of_kmers, k);

    read_iterator.rewind();
    auto* counts = new uint64_t[KMER_COUNT_MAX + 1]();
    run_in_threads(kmer_occurrence_histogram_thread, std::ref(read_iterator), std::ref(filter), std::ref(processed), k, counts);

    Histogram hist;
    for (int i = 5; i <= KMER_COUNT_MAX; i++){
        if (counts[i] < 50) continue;
        hist[i] = counts[i];
    }
    delete [] counts;
    return hist;
}


std::mutex visualization_mut;

void visualize_kmer_positions_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int lower_coverage, int upper_coverage, std::ofstream &output){
    std::optional<GenomeReadData> read;
    std::vector<KmerCount> kmer_counts(filter.categories, 0);

    std::map<UpperSpecificity, char> specificity_color_codes = {
            {70, 'K'},
            {85, 'R'},
            {90, 'O'},
            {95, 'Y'},
            {99, 'B'},
            {101, 'G'}
    };

    KmerQuality threshold_quality = {(Quality)';', (Quality)'C'};

    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, filter.k);

        char specificity_colors[read->sequence.length() + 1];
        uint32_t position_in_read = 0;

        while (it.next_kmer()) {
            KmerCount count = filter.get_count(it.current_kmer);
            if (lower_coverage <= count && count <= upper_coverage){
                for (int category_id = 0; category_id < filter.categories; category_id++) {
                    kmer_counts[category_id] = filter.get_count(it.current_kmer, category_id);
                }
                uint32_t total_count = std::accumulate(kmer_counts.begin(), kmer_counts.end(), 0);
                KmerCount prevalent_count = *std::max_element(kmer_counts.begin(), kmer_counts.end());

                specificity_colors[position_in_read] = (*specificity_color_codes.upper_bound(((double) prevalent_count / (double) total_count) * 100)).second;
            } else specificity_colors[position_in_read] = 'U';

            position_in_read++;
        }
        specificity_colors[position_in_read] = 0;

        visualization_mut.lock();
        output << "\t'" << specificity_colors << "',\n";
        visualization_mut.unlock();
    }
}


void visualize_kmer_positions(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, int lower_coverage, int upper_coverage) {
    std::ofstream output;
    output.open("reads.js");
    output << "reads = [\n";
    read_iterator.rewind();
    run_in_threads(visualize_kmer_positions_thread, std::ref(read_iterator), std::ref(filter), lower_coverage, upper_coverage, std::ref(output));

    output << "]";
    output.close();
}