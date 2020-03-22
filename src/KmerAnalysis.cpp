#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <fmt/format.h>

#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "Utils.h"

namespace algo = boost::algorithm;


std::mutex mut;


void plot_kmer_specificity(std::map<int, KmerSpecificity> &specificities, int max_coverage) {
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
    std::string plot_input = fmt::format("{} {}\n{}", specificities.size(), max_coverage, algo::join(k_spec_strings, "\n"));
    std::cout << run_command_with_input("python python_scripts/plot_histogram.py", plot_input) << std::endl;
}


KmerSpecificity get_kmer_specificity(KmerOccurrences &occurrences) {
    KmerSpecificity specificity;

    std::set<UpperSpecificity> upper_specificity_bounds = {70, 85, 90, 95, 99, 101};
    for (auto s : upper_specificity_bounds) specificity[s] = {};

    KmerOccurrences::iterator iter;
    int progress;
    for (iter = occurrences.begin(), progress = 1; iter != occurrences.end(); iter++, progress++) {
        uint64_t total_count = iter->second.total_occurrences();
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


void kmer_occurrences_thread(SequenceRecordIterator &read_iterator, KmerOccurrences &total_occurrences, int k, int &finish_order) {
    KmerOccurrences partial_occurrences;
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k);
        std::optional<std::pair<Kmer, KmerQuality>> kmer_info;

        while ((kmer_info = it.get_next_kmer()) != std::nullopt) {
            if (!partial_occurrences.contains(kmer_info->first)) {
                partial_occurrences[kmer_info->first] = {0, 0, (uint64_t) kmer_info->second.avg_quality};
            }
            partial_occurrences[kmer_info->first].in_first_category += (int) read->category_flag;
            partial_occurrences[kmer_info->first].in_second_category += (int) !read->category_flag;
            partial_occurrences[kmer_info->first].sum_of_qualities += kmer_info->second.avg_quality;
        }
    }

    // Merge into the total occurrences
    mut.lock();
    std::cout << fmt::format("Thread {} merging...\n", finish_order);
    for (auto it = begin(partial_occurrences); it != end(partial_occurrences); it++) {
        if (!total_occurrences.contains(it->first)) {
            total_occurrences[it->first] = {0, 0, 0};
        }
        total_occurrences[it->first].in_first_category += it->second.in_first_category;
        total_occurrences[it->first].in_second_category += it->second.in_second_category;
        total_occurrences[it->first].sum_of_qualities += it->second.sum_of_qualities;
    }
    finish_order++;
    mut.unlock();
}


KmerOccurrences kmer_occurrences(SequenceRecordIterator &read_iterator, int k) {
    KmerOccurrences occurrences;

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    int finish_order = 0;
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(kmer_occurrences_thread, std::ref(read_iterator), std::ref(occurrences), k, std::ref(finish_order));
    }

    for (int i = 0; i < num_threads; ++i) {
        t[i].join();
    }

    return occurrences;
}
