#include <iostream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/thread.hpp>
#include <fmt/format.h>
#include <tsl/robin_set.h>

#include "../common/KmerCountingBloomFilter.h"
#include "../common/SequenceRecordIterator.h"
#include "../common/KmerIterator.h"
#include "../common/KmerAnalysis.h"
#include "../lib/HyperLogLog.hpp"


namespace po = boost::program_options;


void approximate_kmer_count_thread(SequenceRecordIterator &read_iterator, int k, hll::HyperLogLog &hyper){
    std::optional<GenomeReadData> read;
    while ((read = read_iterator.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);
        while (it.next_kmer()) {
            hyper.add((void*)&it.current_kmer, sizeof(Kmer));
        }
    }
}



uint32_t get_true_number_of_kmers(SequenceRecordIterator &read_iterator, int k) {
    auto hyper = hll::HyperLogLog(10);
    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(approximate_kmer_count_thread, std::ref(read_iterator), k, std::ref(hyper));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return hyper.estimate();
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

//    int max_genome_size = 0, max_coverage = 0;
//    if (vm.count("coverage") || vm.count("genome")) {
//        if (vm.count("coverage")) {
//            max_coverage = vm["coverage"].as<int>();
//
//            max_genome_size = get_genome_size(read_iterator.meta, max_coverage);
//            std::cout << fmt::format("Determined genome size {}\n", max_genome_size);
//        }
//        if (vm.count("genome")) {
//            max_genome_size = vm["genome"].as<int>();
//            if (max_coverage == 0) {
//                max_coverage = get_coverage(read_iterator.meta, max_genome_size);
//            }
//            std::cout << fmt::format("Determined coverage {}\n", max_coverage);
//        }
//    } else {
//        throw std::invalid_argument("Please specify either the estimated genome size or coverage");
//    }
//
//    uint32_t estimated_kmer_count = get_num_of_expected_kmers(k, max_genome_size, max_coverage, 150, 0.01);
//    std::cout << fmt::format("Expecting {} kmers\n", estimated_kmer_count);
//    auto bf = KmerCountingBloomFilter(estimated_kmer_count, 0.001);
//    std::cout << bf.actual_size << std::endl;
//    std::cout << bf.hash_count << std::endl;
//
//
//    boost::posix_time::ptime mst1 = boost::posix_time::microsec_clock::local_time();
//    std::thread t[8];
//    for (auto &i : t) {
//        i = std::thread(occurrence_thread, std::ref(read_iterator), std::ref(bf), k);
//    }
//
//    for (auto &i : t) {
//        i.join();
//    }
//    boost::posix_time::ptime mst2 = boost::posix_time::microsec_clock::local_time();
//    boost::posix_time::time_duration msdiff = mst2 - mst1;
//    std::cout << fmt::format("Occurrences computed in {} ms. Counted {} of {} possible kmers\n", msdiff.total_milliseconds(), bf.cardinality(), pow(4, k));
//    std::cout << fmt::format("BF actual size {} cardinality {}\n", bf.actual_size, bf.cardinality());
    std::cout << fmt::format("Number of kmers is {}\n", get_true_number_of_kmers(read_iterator, k));
}