#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <fmt/format.h>

#include "../common/SequenceRecordIterator.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"
#include "../common/KmerAnalysis.h"

#include "ReadClusteringEngine.h"
#include "KmerAnalysis.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("genome,g", po::value<int >(), "Estimated size of the larger genome")
            ("coverage,c", po::value<int >(), "Estimated coverage (for one read file)")
            ("k-size,k", po::value<int> (), "Size of kmer to analyze & select")
            ("cov-lower,l", po::value<int> (), "Lower bound for coverage of characteristic kmers")
            ("cov-upper,u", po::value<int> (), "Upper bound for coverage of characteristic kmers");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    if (read_paths.empty()) {
        throw std::invalid_argument("You need to specify paths to read files");
    }

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths);

    int max_genome_size = 0, max_coverage = 0;
    if (vm.count("coverage") || vm.count("genome")){
        if (vm.count("coverage")){
            max_coverage = vm["coverage"].as<int>();

            max_genome_size = get_genome_size(read_iterator.file_meta, max_coverage);
            std::cout << fmt::format("Determined genome size {}\n", max_genome_size);
        }
        if (vm.count("genome")){
            max_genome_size = vm["genome"].as<int>();
            if (max_coverage == 0){
                max_coverage = get_coverage(read_iterator.file_meta, max_genome_size);
            }
            std::cout << fmt::format("Determined coverage {}\n", max_coverage);
        }
    } else {
        throw std::invalid_argument("Please specify either the estimated genome size or coverage");
    }

    int selected_k = 0, cov_lower = 0, cov_upper = 0;
    if (vm.count("k-size")) selected_k = vm["k-size"].as<int>();
    if (vm.count("cov-lower")) cov_lower = vm["cov-lower"].as<int>();
    if (vm.count("cov-upper")) cov_upper = vm["cov-upper"].as<int>();

    auto expected_num_of_kmers = get_approximate_kmer_count(read_iterator, selected_k);

    std::cout << fmt::format("Expecting {} kmers\n", expected_num_of_kmers);
    KmerCountingBloomFilter bf = kmer_occurrences(read_iterator, selected_k, expected_num_of_kmers);
    kmer_occurrence_histogram(read_iterator, bf, selected_k, expected_num_of_kmers);

    auto engine = ReadClusteringEngine(read_iterator, bf, selected_k, cov_lower, cov_upper);
    engine.run_clustering();
    return 0;
}
