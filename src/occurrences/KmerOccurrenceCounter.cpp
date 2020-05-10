#include "KmerOccurrenceCounter.h"
#include "KmerAnalysis.h"
#include "../common/Utils.h"
#include "../common/KmerIterator.h"


KmerOccurrenceCounter::KmerOccurrenceCounter(SequenceRecordIterator &read_iterator) {
    reader = &read_iterator;
    auto k_and_count = get_unique_k_length(read_iterator);
    k = k_and_count.first;
    expected_num_of_kmers = k_and_count.second;
    compute_filters();
}

KmerOccurrenceCounter::KmerOccurrenceCounter(SequenceRecordIterator &read_iterator, int k) : k(k) {
    reader = &read_iterator;
    expected_num_of_kmers = get_approximate_kmer_count(read_iterator, k);
    compute_filters();
}

bloom::BloomFilter<Kmer> KmerOccurrenceCounter::count_non_singletons(){
    auto first_occurrence = bloom::BloomFilter<Kmer>(expected_num_of_kmers, 0.01);
    auto second_occurrence = bloom::BloomFilter<Kmer>(expected_num_of_kmers, 0.01);

    auto count_non_singletons_thread = [this, &first_occurrence, &second_occurrence](){
        std::optional<GenomeReadData> read;
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            while (it.next_kmer()) {
                if (!first_occurrence.insert(it.current_kmer)){
                    second_occurrence.insert(it.current_kmer);
                }
            }
        }
    };

    reader->rewind();
    run_in_threads(count_non_singletons_thread);

    return second_occurrence;
}

void KmerOccurrenceCounter::compute_filters() {

    std::cout << "Counting kmers with more than one occurrence\n";
    auto more_than_once = count_non_singletons();
    auto more_than_once_count = more_than_once.cardinality();

    std::cout << "Computing filters\n";
    filters.clear();
    for (int temp = 0; temp < reader->categories; temp++) {
        filters.push_back(new counting_bloom::CountingBloomFilter<Kmer, uint16_t>(more_than_once_count, 0.05));
    }
    reader->rewind();

    auto compute_filters_thread = [this, &more_than_once](){
        std::optional<GenomeReadData> read;
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            while (it.next_kmer()) {
                if (more_than_once.contains(it.current_kmer)){
                    filters[read->category_id]->add(it.current_kmer);
                }
            }
        }
    };
    run_in_threads(compute_filters_thread);
}

KmerSpecificity KmerOccurrenceCounter::get_specificity(std::vector<double> &thresholds) {

    std::cout << "Calculating specificity\n";
    specificity.clear();
    for (auto threshold : thresholds) {
        specificity.insert({threshold, {}});
    }

    bloom::BloomFilter<Kmer> processed(expected_num_of_kmers, 0.01);
    std::mutex mut;
    reader->rewind();

    auto get_specificity_thread = [this, &processed, &mut](){
        std::set<double> threshold_set;
        for (auto it = begin(specificity); it != end(specificity); ++it) {
            threshold_set.insert(it->first);
        }
        std::optional<GenomeReadData> read;
        std::vector<uint16_t> kmer_counts(reader->categories, 0);
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            while (it.next_kmer()) {
                if (processed.insert(it.current_kmer)) {
                    uint16_t prevalent_count = 0;
                    uint16_t total_count = 0;
                    for (auto filter : filters){
                        auto count = filter->get_count(it.current_kmer);
                        prevalent_count = std::max(count, prevalent_count);
                        total_count += count;
                    }
                    UpperSpecificity kmer_specificity = *threshold_set.upper_bound(((double) prevalent_count / (double) total_count) * 100);
                    mut.lock();
                    specificity[kmer_specificity].insert(std::pair<NumOfOccurrences, UniqueKmerCount>(total_count, 0)).first->second += 1;
                    mut.unlock();
                }
            }
        }
    };
    run_in_threads(get_specificity_thread);
    return specificity;
}

Histogram KmerOccurrenceCounter::get_histogram() {
    std::cout << "Calulating histogram\n";

    auto *counts = new uint64_t[UINT16_MAX + 1]();

    bloom::BloomFilter<Kmer> processed(expected_num_of_kmers, 0.01);
    std::mutex mut;

    reader->rewind();
    auto get_histogram_thread = [this, counts, &processed, &mut](){
        std::optional<GenomeReadData> read;
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            while (it.next_kmer()) {
                if (processed.insert(it.current_kmer)) {
                    uint16_t total_count = 0;
                    for (auto filter : filters){
                        total_count += filter->get_count(it.current_kmer);
                    }
                    ++*(counts + total_count);
                }
            }
        }
    };
    run_in_threads(get_histogram_thread);

    histogram.clear();
    for (int i = 1; i <= UINT16_MAX; i++) {
        if (counts[i] > 100) histogram[i] = counts[i];
    }
    delete[] counts;

    return histogram;
}

void KmerOccurrenceCounter::export_kmers_in_range(int lower_bound, int upper_bound, std::string &path) {
    if (histogram.empty()) get_histogram();

    if (path.empty()){
        path = fmt::format("./data/kmers/{}_{}-{}__kmers.bin", reader->meta.filename, lower_bound, upper_bound);
    }

    uint64_t exported_kmer_count = 0;
    for (auto it = histogram.lower_bound(lower_bound); it != histogram.upper_bound(upper_bound + 1); it++) {
        exported_kmer_count += it->second;
    }
    std::cout << "Exporting ~" << exported_kmer_count << " kmers\n";

    bloom::BloomFilter<Kmer> exported_kmers(exported_kmer_count, 0.01);
    bloom::BloomFilter<Kmer> processed(expected_num_of_kmers, 0.01);
    reader->rewind();

    auto export_kmers_thread = [this, &exported_kmers, &processed, lower_bound, upper_bound](){
        std::optional<GenomeReadData> read;
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            while (it.next_kmer()) {
                if (!processed.insert(it.current_kmer)) continue;

                uint16_t total_count = 0;
                for (int category_id = 0; category_id < reader->categories; category_id++) {
                    total_count += filters[category_id]->get_count(it.current_kmer);
                }
                if (lower_bound <= total_count && total_count <= upper_bound) {
                    exported_kmers.add(it.current_kmer);
                }
            }
        }
    };
    run_in_threads(export_kmers_thread);

    auto out = std::ofstream(path, std::ios::out | std::ios::binary);
    out.write((char*)&k, sizeof(k));
    exported_kmers.dump(out);
    out.close();
}
