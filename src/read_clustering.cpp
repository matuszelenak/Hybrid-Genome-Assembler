#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <fmt/format.h>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"
#include "common/KmerAnalysis.h"
#include "common/Plotting.h"
#include "common/BaseReadClusteringEngine.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("occurrences,o", po::value<std::string> (), "Path to binary file with kmer occurrences CBF")
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
    if (!vm.count("occurrences")) throw std::invalid_argument("You need to specify paths to read files");

    int cov_lower = 0, cov_upper = 0;
    if (vm.count("cov-lower")) cov_lower = vm["cov-lower"].as<int>();
    if (vm.count("cov-upper")) cov_upper = vm["cov-upper"].as<int>();

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths, true);
    std::string bf_path = vm["occurrences"].as<std::string>();
    auto filter = KmerCountingBloomFilter(bf_path);

    uint32_t expected_kmer_count = get_approximate_kmer_count(read_iterator, filter.k);
    KmerSpecificity spec = get_kmer_specificity(read_iterator, filter, expected_kmer_count);
    std::map<int, KmerSpecificity> spec_map = {{filter.k, spec}};
    plot_kmer_specificity(spec_map, 100);
    if (cov_lower == 0){
        std::cout << "Enter the lower bound for occurrences of discriminative kmers" << std::endl;
        std::cin >> cov_lower;
    }
    if (cov_upper == 0){
        std::cout << "Enter the upper bound for occurrences of discriminative kmers" << std::endl;
        std::cin >> cov_upper;
    }
    if (cov_lower > cov_upper) throw std::logic_error("Range must be non-empty");

    auto engine = BaseReadClusteringEngine(read_iterator, filter, filter.k, cov_lower, cov_upper);
    engine.run_clustering();
    return 0;
}
