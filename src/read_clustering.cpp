#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"
#include "common/KmerAnalysis.h"
#include "common/ReadClusteringEngine.h"
#include "lib/BloomFilter.h"


namespace po = boost::program_options;


std::pair<tsl::robin_set<Kmer>, int> load_kmers(std::string path){
    std::ifstream in;
    in.open(path);
    std::string kmer;

    int k;
    tsl::robin_set<Kmer> s;
    while (std::getline(in, kmer)){
        k = kmer.length();
        KmerIterator k_it(kmer, k);
        k_it.next_kmer();
        s.insert(k_it.current_kmer);
//        std::cout << k_it.number_to_sequence(k_it.current_kmer) << std::endl;
//        std::cin >> kmer;
    }
    return {s, k};
}

int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("bloom,b", po::value<std::string>(), "Path to binary file with discriminative kmers")
            ("kmers,k", po::value<std::string>(), "Path to text file with kmers")
            ("platform,p", po::value<std::string>(), "Specify origin platform of the reads [pacbio|nanopore]");

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

    Platform platform;
    if (vm.count("platform")){
        auto platform_str = vm["platform"].as<std::string>();
        if (platform_str == "pacbio") platform = PacBio;
        else if (platform_str == "nanopore") platform = Nanopore;
        else throw std::invalid_argument("Invalid platform specified");
    } else platform = PacBio;

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths, true);
    for (auto meta: read_iterator.file_meta){
        std::cout << meta.repr();
    }

    if (vm.count("bloom")){
        auto in = std::ifstream(vm["bloom"].as<std::string>(), std::ios::in | std::ios::binary);
        int k;
        in.read((char *) &k, sizeof(k));
        auto kmers = bloom::BloomFilter<Kmer>(in);

        auto engine = ReadClusteringEngine(read_iterator, k, kmers, platform);
        engine.run_clustering();
    } else if (vm.count("kmers")){
        auto set_k_pair = load_kmers(vm["kmers"].as<std::string>());

        std::cout << fmt::format("Loaded {} kmers of length {}\n", set_k_pair.first.size(), set_k_pair.second);
        auto engine = ReadClusteringEngine(read_iterator, set_k_pair.second, set_k_pair.first, platform);
        engine.run_clustering();
    } else {
        throw std::invalid_argument("You need to specify path to kmers");
    }
    return 0;
}
