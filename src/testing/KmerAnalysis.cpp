#include <iostream>
#include <fmt/format.h>
#include <numeric>

#include "../common/KmerAnalysis.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"
#include "KmerAnalysis.h"

namespace pt = boost::posix_time;

std::mutex mut;
std::mutex visualization_mut;
std::set<UpperSpecificity> upper_specificity_bounds = {70, 85, 90, 95, 99, 101};
std::map<UpperSpecificity, char> specificity_color_codes = {
        {70, 'K'},
        {85, 'R'},
        {90, 'O'},
        {95, 'Y'},
        {99, 'B'},
        {101, 'G'}
};

void kmer_specificity_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, KmerBloomFilter &processed, KmerSpecificity &specificity){
    std::optional<GenomeReadData> read;
    std::vector<KmerCount> kmer_counts(filter.categories, 0);

    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            if (!processed.contains(it.current_kmer)) {
                processed.add(it.current_kmer);

                for (int category_id = 0; category_id < filter.categories; category_id++) {
                    kmer_counts[category_id] = filter.get_count(it.current_kmer, category_id);
                }
                uint32_t total_count = std::accumulate(kmer_counts.begin(), kmer_counts.end(), 0);
                KmerCount prevalent_count = *std::max_element(kmer_counts.begin(), kmer_counts.end());
                UpperSpecificity kmer_specificity = *upper_specificity_bounds.upper_bound(((double) prevalent_count / (double) total_count) * 100);

                mut.lock();
                specificity[kmer_specificity].insert(std::pair<NumOfOccurrences, UniqueKmerCount>(total_count, 0)).first->second += 1;
                mut.unlock();
            }
        }
    }
}


KmerSpecificity get_kmer_specificity(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k) {
    KmerSpecificity specificity = {
            {70,  {}},
            {85,  {}},
            {90,  {}},
            {95,  {}},
            {101, {}}
    };
    KmerBloomFilter processed = KmerBloomFilter(filter._expected_items, 0.001);
    read_iterator.reset();
    auto runner = ThreadRunner(kmer_specificity_thread, std::ref(read_iterator), std::ref(filter), k, std::ref(processed), std::ref(specificity));

    return specificity;
}


void visualize_kmer_positions_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &filter, int k, int lower_coverage, int upper_coverage, std::ofstream &output){
    std::optional<GenomeReadData> read;
    std::vector<KmerCount> kmer_counts(filter.categories, 0);

    KmerQuality threshold_quality = {(Quality)';', (Quality)'C'};

    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        auto q_it = KmerQualityIterator(read->qualities, k);

        char specificity_colors[read->sequence.length() + 1];
        uint32_t position_in_read = 0;

        while (it.next_kmer() && q_it.next_quality()) {
            KmerCount count = filter.get_count(it.current_kmer);
            if (lower_coverage <= count && count <= upper_coverage && q_it.current_quality > threshold_quality ){
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
    read_iterator.reset();
    auto runner = ThreadRunner(visualize_kmer_positions_thread, std::ref(read_iterator), std::ref(filter), k, lower_coverage, upper_coverage, std::ref(output));

    output << "]";
    output.close();
}