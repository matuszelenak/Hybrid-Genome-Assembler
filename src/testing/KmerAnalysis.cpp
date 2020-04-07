#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <numeric>

#include "../common/KmerAnalysis.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"
#include "KmerAnalysis.h"

namespace algo = boost::algorithm;
namespace pt = boost::posix_time;


std::mutex mut;


void plot_kmer_specificity(std::vector<ReadFileMetaData> &meta, std::map<int, KmerSpecificity> &specificities, int max_coverage) {
    // Lord forgive me for what I am about to do
    std::vector<std::string> k_spec_strings;
    for (const auto &k_specs : specificities) {
        std::vector<std::string> bound_strings;
        for (const auto &bound_specs : k_specs.second) {
            std::vector<std::string> counts;
            for (const auto &coverage_counts : bound_specs.second) {
                if (coverage_counts.first < (double)max_coverage / 6.0) continue;
                counts.push_back(fmt::format("({}, {})", coverage_counts.first, coverage_counts.second));
            }
            bound_strings.push_back(fmt::format("({}, [{}])", bound_specs.first, algo::join(counts, ", ")));
        }
        k_spec_strings.push_back(fmt::format("({}, [{}])", k_specs.first, algo::join(bound_strings, ", ")));
    }
    std::vector<std::string> filenames;
    std::transform(meta.begin(), meta.end(), std::back_inserter(filenames), [](ReadFileMetaData &d) -> std::string { return d.filename; });
    std::string plot_input = fmt::format("{}\n{} {}\n{}", algo::join(filenames, "/"), specificities.size(), max_coverage, algo::join(k_spec_strings, "\n"));
    std::cout << run_command_with_input("python3 common/python/plot_histogram.py", plot_input) << std::endl;
}

std::vector<uint8_t> get_k_sizes(int max_genome_size) {
    std::vector<uint8_t> k_sizes;
    uint8_t k_guess = (uint8_t) ceil(log(max_genome_size) / log(4)) | 1u;
    for (int k_value = k_guess; k_value <= 19; k_value+=2) {
        k_sizes.push_back(k_value);
    }
    return k_sizes;
}


KmerSpecificity get_kmer_specificity(SequenceRecordIterator &read_iterator, std::vector<KmerCountingBloomFilter*> &filters, int k, uint32_t expected_num_of_kmers){
    KmerSpecificity specificity;

    std::set<UpperSpecificity> upper_specificity_bounds = {70, 85, 90, 95, 99, 101};
    for (auto s : upper_specificity_bounds) specificity[s] = {};

    KmerBloomFilter processed = KmerBloomFilter(expected_num_of_kmers * 2, 0.01);

    read_iterator.reset();
    std::optional<GenomeReadData> read;
    std::vector<KmerCount> kmer_counts(filters.size(), 0);

    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            if (!processed.contains(it.current_kmer)){
                processed.add(it.current_kmer);

                for (int filter_id = 0; filter_id < filters.size(); filter_id++){
                    kmer_counts[filter_id] = filters[filter_id]->get_count(it.current_kmer);
                }
                uint32_t total_count = std::accumulate(kmer_counts.begin(), kmer_counts.end(), 0);
                if (total_count < 3) continue;

                KmerCount prevalent_count = *std::max_element(kmer_counts.begin(), kmer_counts.end());
                UpperSpecificity kmer_specificity = *upper_specificity_bounds.upper_bound(((double) prevalent_count / (double) total_count) * 100);

                specificity[kmer_specificity].insert(std::pair<NumOfOccurrences, UniqueKmerCount>(total_count, 0)).first->second += 1;
            }
        }
    }
    return specificity;
}



void kmer_occurrences_first_pass(SequenceRecordIterator &read_iterator, std::vector<KmerCountingBloomFilter*> &filters, int k){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            filters[read->category_id]->add(it.current_kmer);
        }
    }
}


std::vector<KmerCountingBloomFilter*> kmer_occurrences(SequenceRecordIterator &read_iterator, int k, uint32_t expected_num_of_kmers){
    std::vector<KmerCountingBloomFilter*> filters;
    for (int category = 0; category < read_iterator.file_meta.size(); category++){
        filters.push_back(new KmerCountingBloomFilter(expected_num_of_kmers));
    }

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    // First pass to fill the bloom filter
    read_iterator.reset();
    pt::ptime ts_start = pt::microsec_clock::local_time();
    for (auto &i : t) i = std::thread(kmer_occurrences_first_pass, std::ref(read_iterator), std::ref(filters), k);
    for (auto &i : t) i.join();
    std::cout << fmt::format("First pass done in {}ms\n", (pt::microsec_clock::local_time() - ts_start).total_milliseconds());

    return filters;
}
