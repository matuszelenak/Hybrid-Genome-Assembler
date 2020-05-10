#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <iostream>
#include "occurrences/JellyfishOccurrenceReader.h"
#include "common/Plotting.h"

namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("k-size,k", po::value<int>(), "Size of kmer to analyze & select")
            ("output,o", po::value<std::string> (), "Output path for the counting bloom filter");

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

    std::string output_path = "mers.txt";
    if (vm.count("output")) {
        output_path = vm["output"].as<std::string>();
    }

    auto counter = JellyfishOccurrenceReader(read_paths, vm["k-size"].as<int>());

    std::set<double> spec_thresholds = {70, 85, 90, 95, 99, 101};
    KmerSpecificity spec = counter.get_specificity(spec_thresholds);
    std::map<int, KmerSpecificity> spec_map = {{vm["k-size"].as<int>(), spec}};
    plot_kmer_specificity(spec_map, 200);

    int lower, upper;
    std::cout << "Enter lower and upper bounds for exported kmers\n";
    std::cin >> lower >> upper;

    counter.export_kmers(lower, upper, output_path);
    return 0;
}