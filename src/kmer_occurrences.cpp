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
#include "common/KmerOccurrenceCounter.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("k-size,k", po::value<int> (), "Size of kmer to analyze & select")
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

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths, true);
    auto occurrence_counter = KmerOccurrenceCounter();

    int selected_k;
    if (vm.count("k-size")){
        occurrence_counter = KmerOccurrenceCounter(read_iterator, vm["k-size"].as<int>());
    } else {
        occurrence_counter = KmerOccurrenceCounter(read_iterator);
    }

    std::string output_path;
    if (vm.count("output")) {
        output_path = vm["output"].as<std::string>();
    } else {
        output_path = read_iterator.meta.filename + "__kmers.bin";
    }

    if (read_iterator.categories == 1){
        auto hist = occurrence_counter.get_histogram();
        plot_histogram(hist);
    } else {
        std::vector<double> spec_thresholds = {70, 85, 90, 95, 99, 101};
        KmerSpecificity spec = occurrence_counter.get_specificity(spec_thresholds);
        std::map<int, KmerSpecificity> spec_map = {{selected_k, spec}};
        plot_kmer_specificity(spec_map, 100);
    }

    int lower, upper;
    std::cout << "Enter lower and upper bounds for exported kmers\n";
    std::cin >> lower >> upper;

    occurrence_counter.export_kmers_in_range(lower, upper, output_path);
    return 0;
}
