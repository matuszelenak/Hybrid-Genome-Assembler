#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"
#include "common/KmerAnalysis.h"
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

    int selected_k;
    uint64_t expected_num_of_kmers;
    if (vm.count("k-size")){
        selected_k = vm["k-size"].as<int>();
        expected_num_of_kmers = get_approximate_kmer_count(read_iterator, selected_k);
    } else {
        auto k_and_count = get_unique_k_length(read_iterator);
        selected_k = k_and_count.first;
        expected_num_of_kmers = k_and_count.second;
    }

    std::cout << fmt::format("Expecting {} kmers\n", expected_num_of_kmers);
    KmerCountingBloomFilter bf = kmer_occurrence_filter(read_iterator, selected_k, expected_num_of_kmers);

    Histogram hist = kmer_occurrence_histogram(read_iterator, bf, selected_k, expected_num_of_kmers);
    plot_histogram(hist);

    std::string output_path;
    if (vm.count("output")) {
        output_path = vm["output"].as<std::string>();
    } else {
        std::vector<std::string> filenames;
        std::transform(read_iterator.file_meta.begin(), read_iterator.file_meta.end(), std::back_inserter(filenames), [](ReadFileMetaData &meta) -> std::string { return meta.filename; });
        output_path = boost::algorithm::join(filenames, "__") + "__filter.bin";
    }

    bf.dump(output_path);
    return 0;
}
