#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <tsl/robin_set.h>

#include "KmerCountingBloomFilter.h"
#include "SequenceRecordIterator.h"
#include "KmerIterator.h"
#include "KmerAnalysis.h"


namespace po = boost::program_options;


std::mutex cntr_mutex;


void occurrence_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k) {
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k, true);
        while (it.next_kmer()) {
            bf.add(it.current_kmer);
        }
    }
}


void true_cardinality_thread(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, tsl::robin_set<Kmer> &counter) {
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k, true);
        while (it.next_kmer()) {
            if (bf.contains(it.current_kmer)) {
                cntr_mutex.lock();
                counter.insert(it.current_kmer);
                cntr_mutex.unlock();
            }
        }
    }
}

uint32_t true_cardinality(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k) {
    tsl::robin_set<Kmer> counter;
    std::thread t[8];
    read_iterator.reset();
    for (auto &i : t) {
        i = std::thread(true_cardinality_thread, std::ref(read_iterator), std::ref(bf), k, std::ref(counter));
    }

    for (auto &i : t) {
        i.join();
    }
    return counter.size();
}


uint32_t get_true_number_of_kmers(SequenceRecordIterator &read_iterator, int k) {
    std::set<Kmer> s;
    std::optional<GenomeReadData> read;
    read_iterator.reset();
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k, true);
        while (it.next_kmer()) {
            s.insert(it.current_kmer);
        }
    }
    return s.size();
}


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("genome,g", po::value<int>(), "Estimated size of the larger genome")
            ("coverage,c", po::value<int>(), "Estimated coverage (for one read file)")
            ("k-size,k", po::value<int>(), "Size of kmer to analyze & select");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (read_paths.empty()) {
        throw std::invalid_argument("You need to specify paths to read files");
    }

    int k = vm["k-size"].as<int>();

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths);

    int max_genome_size = 0, max_coverage = 0;
    if (vm.count("coverage") || vm.count("genome")) {
        if (vm.count("coverage")) {
            max_coverage = vm["coverage"].as<int>();

            max_genome_size = get_genome_size(read_iterator.meta, max_coverage);
            std::cout << fmt::format("Determined genome size {}\n", max_genome_size);
        }
        if (vm.count("genome")) {
            max_genome_size = vm["genome"].as<int>();
            if (max_coverage == 0) {
                max_coverage = get_coverage(read_iterator.meta, max_genome_size);
            }
            std::cout << fmt::format("Determined coverage {}\n", max_coverage);
        }
    } else {
        throw std::invalid_argument("Please specify either the estimated genome size or coverage");
    }

    uint32_t estimated_kmer_count = get_num_of_expected_kmers(k, max_genome_size, max_coverage, 150, 0.01);
    std::cout << fmt::format("Expecting {} kmers\n", estimated_kmer_count);
    auto bf = KmerCountingBloomFilter(estimated_kmer_count, 0.001);
    std::cout << bf.actual_size << std::endl;
    std::cout << bf.hash_count << std::endl;


    boost::posix_time::ptime mst1 = boost::posix_time::microsec_clock::local_time();
    std::thread t[8];
    for (auto &i : t) {
        i = std::thread(occurrence_thread, std::ref(read_iterator), std::ref(bf), k);
    }

    for (auto &i : t) {
        i.join();
    }
    boost::posix_time::ptime mst2 = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration msdiff = mst2 - mst1;
    std::cout << fmt::format("Occurrences computed in {} ms. Counted {} of {} possible kmers\n", msdiff.total_milliseconds(), bf.cardinality(), pow(4, k));
    std::cout << fmt::format("BF actual size {} cardinality {}\n", bf.actual_size, bf.cardinality());
    std::cout << fmt::format("BF true cardinality is {}\n", true_cardinality(read_iterator, bf, k));
}