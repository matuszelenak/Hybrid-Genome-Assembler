#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <iostream>

#include "common/SequenceRecordIterator.h"
#include "common/KmerAnalysis.h"
#include "common/KmerCountingBloomFilter.h"
#include "common/KmerIterator.h"

namespace po = boost::program_options;


int main(int argc, char *argv[]){
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("references", -1);
    std::vector<std::string> reference_paths;

    desc.add_options()
            ("help,h", "Help screen")
            ("references", po::value<std::vector<std::string>>(&reference_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("bloom,b", po::value<std::string> (), "Path to the CBF")
            ("output,o", po::value<std::string> (), "Output path ")
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

    std::string kmers_path = vm["bloom"].as<std::string>();
    auto kmers = KmerCountingBloomFilter(kmers_path);

    int cov_lower = 0, cov_upper = 0;
    if (vm.count("cov-lower")) cov_lower = vm["cov-lower"].as<int>();
    if (vm.count("cov-upper")) cov_upper = vm["cov-upper"].as<int>();

    std::ofstream output;
    output.open("positions_in_reference.js");
    output << "positions = [\n";
    for (auto reference_path : reference_paths){
        auto reader = SequenceRecordIterator(reference_path);
        visualize_kmer_positions_thread(reader, kmers, cov_lower, cov_upper, output);
    }
    output << "]";
    output.close();
}