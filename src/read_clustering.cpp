#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"
#include "common/KmerAnalysis.h"
#include "common/Plotting.h"
#include "common/BaseReadClusteringEngine.h"
#include "lib/BloomFilter.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("kmers,k", po::value<std::string>(), "Path to binary file with discriminative kmers");

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
    if (!vm.count("kmers")) throw std::invalid_argument("You need to specify path to kmers");

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths, true);

    auto in = std::ifstream(vm["kmers"].as<std::string>(), std::ios::in | std::ios::binary);
    int k;
    in.read((char *) &k, sizeof(k));
    auto kmers = bloom::BloomFilter<Kmer>(in);

    auto engine = BaseReadClusteringEngine(read_iterator, k, kmers);
    engine.run_clustering();
    return 0;
}
