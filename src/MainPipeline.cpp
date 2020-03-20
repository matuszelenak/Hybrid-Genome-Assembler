#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <jsoncpp/json/json.h>
#include <fmt/format.h>

#include "ReadDataLoader.h"
#include "DNAStructures.h"
#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "ReadClustering.h"
#include "Utils.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_path", 1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_path", po::value<std::string>(), "Path to file with reads (FASTA or FASTQ)");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    std::string read_path, meta_path, out_path;

    if (vm.count("read_path")) {
        read_path = vm["read_path"].as<std::string>();
        meta_path = read_path + "_meta.json";
    } else throw std::invalid_argument("You need to specify a read file");

    Json::Value meta;
    std::ifstream file(meta_path);
    file >> meta;
    std::cout << meta << std::endl;

    int k_guess = (int)ceil(log(meta["genome_size"].asUInt64()) / log(4)) + 3;

    std::map<int, KmerSpecificity> specificities = {};
    std::map<int, KmerOccurrences> per_k_occurrences;
    for (int k_length = k_guess; k_length < k_guess + 2; k_length++) {
        std::cout << k_length << std::endl;

        ReadDataLoader loader = ReadDataLoader(read_path);
        KmerOccurrences occurrences = kmer_occurrences(loader, k_length);
        specificities[k_length] = get_kmer_specificity(occurrences);
        per_k_occurrences[k_length] = occurrences;
    }

    plot_kmer_specificity(specificities, meta["coverage"].asInt() * 3);
    int selected_k, cov_low, cov_high;
    std::cin >> selected_k >> cov_low >> cov_high;
    std::cout << fmt::format("{} total kmers\n", per_k_occurrences[selected_k].size());
    KmerOccurrences characteristic_kmers = filter_characteristic_kmers(per_k_occurrences[selected_k], cov_low, cov_high);

    std::cout << fmt::format("{} characteristic kmers\n", characteristic_kmers.size());

    ReadDataLoader loader = ReadDataLoader(read_path);
    std::vector<GenomeReadCluster> initial_clusters = get_initial_read_clusters(loader, selected_k, characteristic_kmers);
    run_clustering(initial_clusters);
    return 0;
}
