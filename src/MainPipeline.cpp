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
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    for (auto &s : read_paths){
        std::cout << s << std::endl;
    }

    int max_genome_size, max_coverage;

    if (read_paths.empty()) {
        throw std::invalid_argument("You need to specify a read file");
    } else {
        std::vector<int> genome_sizes, coverages;
        for (auto &path : read_paths){
            Json::Value meta;
            std::ifstream file(path.substr(0, path.find_last_of('.')) + "_meta.json");
            std::cout << fmt::format("Loading meta from {}\n", path.substr(0, path.find_last_of('.')) + "_meta.json");
            file >> meta;
            std::cout << meta << std::endl;

            genome_sizes.push_back(meta["genome_size"].asInt());
            coverages.push_back(meta["coverage"].asInt());
        }

        max_genome_size = *std::max_element(genome_sizes.begin(), genome_sizes.end());
        max_coverage = *std::max_element(coverages.begin(), coverages.end());
    }

    int k_guess = (int)ceil(log(max_genome_size) / log(4)) + 1;

    std::map<int, KmerSpecificity> specificities = {};
    std::map<int, KmerOccurrences> per_k_occurrences;
    for (int k_length = k_guess; k_length < k_guess + 4; k_length++) {
        std::cout << k_length << std::endl;

        SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths);
        KmerOccurrences occurrences = kmer_occurrences(read_iterator, k_length);
        specificities[k_length] = get_kmer_specificity(occurrences);
        per_k_occurrences[k_length] = occurrences;
    }

    plot_kmer_specificity(specificities, max_coverage * 2);
    int selected_k, cov_low, cov_high;
    std::cin >> selected_k >> cov_low >> cov_high;
    std::cout << fmt::format("{} total kmers\n", per_k_occurrences[selected_k].size());
    KmerOccurrences characteristic_kmers = filter_characteristic_kmers(per_k_occurrences[selected_k], cov_low, cov_high);

    std::cout << fmt::format("{} characteristic kmers\n", characteristic_kmers.size());

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths);
    std::vector<GenomeReadCluster> initial_clusters = get_initial_read_clusters(read_iterator, selected_k, characteristic_kmers);
    run_clustering(initial_clusters);
    return 0;
}
