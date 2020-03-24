#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <cmath>

#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "Utils.h"

namespace algo = boost::algorithm;


std::mutex mut;

size_t PARTIAL_MERGE_AFTER = 100000;


void plot_kmer_specificity(std::vector<ReadFileMetaData> &meta, std::map<int, KmerSpecificity> &specificities, int max_coverage) {
    // Lord forgive me for what I am about to do
    std::vector<std::string> k_spec_strings;
    for (const auto &k_specs : specificities) {
        std::vector<std::string> bound_strings;
        for (const auto &bound_specs : k_specs.second) {
            std::vector<std::string> counts;
            for (const auto &coverage_counts : bound_specs.second) {
                counts.push_back(fmt::format("({}, {})", coverage_counts.first, coverage_counts.second));
            }
            bound_strings.push_back(fmt::format("({}, [{}])", bound_specs.first, algo::join(counts, ", ")));
        }
        k_spec_strings.push_back(fmt::format("({}, [{}])", k_specs.first, algo::join(bound_strings, ", ")));
    }
    std::vector<std::string> filenames;
    std::transform(meta.begin(), meta.end(), std::back_inserter(filenames), [](ReadFileMetaData &d) -> std::string { return d.filename; });
    std::string plot_input = fmt::format("{}\n{} {}\n{}", algo::join(filenames, "/"), specificities.size(), max_coverage, algo::join(k_spec_strings, "\n"));
    std::cout << run_command_with_input("python python_scripts/plot_histogram.py", plot_input) << std::endl;
}


KmerSpecificity get_kmer_specificity(KmerOccurrences &occurrences) {
    KmerSpecificity specificity;

    std::set<UpperSpecificity> upper_specificity_bounds = {70, 85, 90, 95, 99, 101};
    for (auto s : upper_specificity_bounds) specificity[s] = {};

    int progress = 1;
    for (auto iter = begin(occurrences); iter != end(occurrences); iter++, progress++) {
        uint16_t total_count = iter->second.total_occurrences();
        uint64_t prevalent = std::max(iter->second.in_first_category, iter->second.in_second_category);
        UpperSpecificity kmer_specificity = *upper_specificity_bounds.upper_bound(((double) prevalent / (double) total_count) * 100);
        if (!specificity[kmer_specificity].contains(total_count)) {
            specificity[kmer_specificity][total_count] = 0;
        }
        specificity[kmer_specificity][total_count] += 1;
        if (progress % 100 == 0 || progress == occurrences.size())
            show_progress(progress, occurrences.size(), "Computing specificities");
    }
    return specificity;
}


KmerOccurrences filter_characteristic_kmers(KmerOccurrences &occurrences, int coverage_lower_bound, int coverage_upper_bound) {
    for (auto it = begin(occurrences); it != end(occurrences);) {
        if (it->second.total_occurrences() < coverage_lower_bound || it->second.total_occurrences() > coverage_upper_bound) {
            it = occurrences.erase(it);
        } else ++it;
    }
    return occurrences;
}


void _merge_partial_occurrences(KmerOccurrences &merged_into, KmerOccurrences &to_merge) {
    for (auto it = begin(to_merge); it != end(to_merge); it++) {
        KmerOccurrences::iterator iter = merged_into.insert( KmerOccurrences::value_type(it->first, {0, 0, 0}) ).first;
        iter.value().in_first_category += it->second.in_first_category;
        iter.value().in_second_category += it->second.in_second_category;
        iter.value().sum_of_qualities += it->second.sum_of_qualities;
    }
}


void kmer_occurrences_thread(SequenceRecordIterator &read_iterator, KmerOccurrences &total_occurrences, int k, int thread_id, int &merger_id, int num_of_threads) {
    KmerOccurrences partial_occurrences;
    partial_occurrences.rehash(PARTIAL_MERGE_AFTER * 2);
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k);
        std::optional<std::pair<Kmer, KmerQuality>> kmer_info;

        while ((kmer_info = it.get_next_kmer()) != std::nullopt) {
            KmerOccurrences::iterator iter = partial_occurrences.insert( KmerOccurrences::value_type(kmer_info->first, {0, 0, (uint32_t) kmer_info->second.avg_quality}) ).first;
            iter.value().in_first_category += (int) read->category_id;
            iter.value().in_second_category += (int) !read->category_id;
            iter.value().sum_of_qualities += kmer_info->second.avg_quality;
        }

        if (partial_occurrences.size() > PARTIAL_MERGE_AFTER && thread_id == merger_id) {
            mut.lock();
            _merge_partial_occurrences(total_occurrences, partial_occurrences);
            //std::cout << fmt::format("Thread {} merging {} items\n", thread_id, partial_occurrences.size());
            merger_id = (merger_id + 1) % num_of_threads;
            mut.unlock();
            partial_occurrences.clear();
        }
    }
    mut.lock();
    _merge_partial_occurrences(total_occurrences, partial_occurrences);
    mut.unlock();
}


KmerOccurrences kmer_occurrences(SequenceRecordIterator &read_iterator, int k, int max_genome_size) {
    KmerOccurrences occurrences;
    // TODO finetune this formula, turns out it has quite an effect on speed
    int expected_kmers = (int)round(max_genome_size * ((double)k / 10) * 2);
    std::cout << fmt::format("Expecting {} kmers\n", expected_kmers);
    occurrences.rehash((KmerOccurrences::size_type)expected_kmers);

    unsigned int num_threads = std::thread::hardware_concurrency();

    boost::posix_time::ptime mst1 = boost::posix_time::microsec_clock::local_time();
    read_iterator.reset();
    int approved_merger_id = 0;

    std::thread t[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(kmer_occurrences_thread, std::ref(read_iterator), std::ref(occurrences), k, i, std::ref(approved_merger_id), num_threads);
    }

    for (int i = 0; i < num_threads; ++i) {
        t[i].join();
    }

    boost::posix_time::ptime mst2 = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration msdiff = mst2 - mst1;
    std::cout << fmt::format("Occurrences computed in {} ms. Counted {} of {} possible kmers\n", msdiff.total_milliseconds(), occurrences.size(), pow(4, k));

    return occurrences;
}


std::vector<int> get_k_sizes(int max_genome_size) {
    std::vector<int> k_sizes;
    int k_guess = (int) ceil(log(max_genome_size) / log(4));
    for (int k_value = k_guess; k_value < std::min(k_guess + 4, 15); k_value++) {
        k_sizes.push_back(k_value);
    }
    return k_sizes;
}

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
