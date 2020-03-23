#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <jsoncpp/json/json.h>
#include <fmt/format.h>

#include "SequenceRecordIterator.h"
#include "DNAStructures.h"
#include "KmerAnalysis.h"
#include "KmerIterator.h"
#include "ReadClustering.h"
#include "Utils.h"


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

            max_genome_size = get_genome_size(read_iterator, max_coverage);
            std::cout << fmt::format("Determined genome size {}\n", max_genome_size);
        }
        if (vm.count("genome")){
            max_genome_size = vm["genome"].as<int>();
            if (max_coverage == 0){
                max_coverage = get_coverage(read_iterator, max_genome_size);
            }
            std::cout << fmt::format("Determined coverage {}\n", max_coverage);
        }
    } else {
        for (auto &path : read_paths){
            Json::Value meta;
            std::ifstream file(path.substr(0, path.find_last_of('.')) + "_meta.json");
            file >> meta;
            std::cout << meta << std::endl;

            max_genome_size = std::max(max_genome_size, meta["genome_size"].asInt());
            max_coverage = std::max(max_coverage, meta["coverage"].asInt());
        }
    }

    std::vector<int>k_sizes;
    if (vm.count("k-size")){
        k_sizes = {vm["k-size"].as<int>()};
    } else {
        k_sizes = get_k_sizes(max_genome_size);
    }

    std::map<int, KmerSpecificity> per_k_specificities = {};
    for (auto k_length : k_sizes) {
        std::cout << fmt::format("\n ### Running analysis for k-mer size {} ###\n", k_length);
        KmerOccurrences occ = kmer_occurrences(read_iterator, k_length);
        per_k_specificities[k_length] = get_kmer_specificity(occ);
    }

    int selected_k = 0, cov_lower = 0, cov_upper = 0;
    if (vm.count("k-size")) selected_k = vm["k-size"].as<int>();
    if (vm.count("cov-lower")) cov_lower = vm["cov-lower"].as<int>();
    if (vm.count("cov-upper")) cov_upper = vm["cov-upper"].as<int>();

    if (!(selected_k && cov_lower && cov_upper)){
        plot_kmer_specificity(per_k_specificities, max_coverage * 3);

        if (!selected_k){
            std::cout << "Enter the size of k-mer\n";
            std::cin >> selected_k;
        }

        if (!(cov_lower && cov_upper)){
            std::cout << "Enter the lower and upper bound for coverage of characteristic kmers\n";
            std::cin >> cov_lower >> cov_upper;
        }
    }

    KmerOccurrences selected_occurrences = kmer_occurrences(read_iterator, selected_k);
    std::cout << fmt::format("{} total kmers\n", selected_occurrences.size());
    KmerOccurrences characteristic_kmers = filter_characteristic_kmers(selected_occurrences, cov_lower, cov_upper);

    std::cout << fmt::format("{} characteristic kmers\n", characteristic_kmers.size());

    std::vector<GenomeReadCluster> initial_clusters = get_initial_read_clusters(read_iterator, selected_k, characteristic_kmers);
    run_clustering(initial_clusters);
    return 0;
}
