#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <iostream>
#include <fmt/format.h>
#include "occurrences/JellyfishOccurrenceReader.h"
#include "common/Plotting.h"
#include "occurrences/KmerAnalysis.h"

namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;
    std::string output_path;
    int k = 11;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("k-size,k", po::value<int>(&k), "Size of kmer to analyze & select")
            ("output,o", po::value<std::string> (&output_path), "Output path for the counting bloom filter");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    if (read_paths.empty()) throw std::invalid_argument("You need to specify paths to read files");
    if (!vm.contains("k-size")){
        SequenceRecordIterator reader(read_paths, true);
        auto k_count_pair = get_unique_k_length(reader);
        k = k_count_pair.first;
    };

    auto counter = JellyfishOccurrenceReader(read_paths, k);

    std::set<double> spec_thresholds = {70, 85, 90, 95, 99, 100, 100.01};
    KmerSpecificity spec = counter.get_specificity(spec_thresholds);
    std::map<int, KmerSpecificity> spec_map = {{k, spec}};
    plot_kmer_specificity(spec_map, 200);

    int lower, upper;
    double percent;
    std::cout << "Enter lower and upper bounds for exported kmers as well as percentage\n";
    std::cin >> lower >> upper >> percent;

    if (output_path.empty()) output_path = fmt::format("{}-mers_{}_{}_{}%.txt", k, lower, upper, percent * 100);
    counter.export_kmers(lower, upper, percent, output_path);
    return 0;
}