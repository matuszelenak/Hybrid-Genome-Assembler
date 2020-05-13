#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"

#include "clustering/ReadClusteringEngine.h"


namespace po = boost::program_options;


std::pair<std::unordered_set<Kmer>, int> load_text_file_kmers(std::string &path){
    std::ifstream in;
    in.open(path);
    std::string kmer;

    int k;
    in.seekg(0);
    std::unordered_set<Kmer> s;
    while (std::getline(in, kmer)){
        k = kmer.length();
        KmerIterator k_it(kmer, k);
        k_it.next_kmer();
        s.insert(k_it.current_kmer);
    }
    return {s, k};
}

int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;
    std::string kmer_path, output_folder_path;
    ReadClusteringConfig config;
    bool debug;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("kmers,k", po::value<std::string>(&kmer_path), "Path to text file with kmers")
            ("output,o", po::value<std::string>(&output_folder_path), "Path to folder with exported clusters")
            ("core_size", po::value<int>(&config.core_component_min_size), "Minimum size for a core component")
            ("core_conn", po::value<double>(&config.core_forming_connections), "Minimum score for a core component forming connection")
            ("tail_amplification", po::value<ConnectionScore>(&config.tail_amplification_min_score), "Minimal score for tail amplifying connections")
            ("core_enrichment", po::value<ConnectionScore>(&config.enrichment_connections_min_score), "Minimal score for connections enriching core components")
            ("debug,d", po::bool_switch(&debug), "Debug flag. Treat read files as separate haplotype reads.");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }
    if (read_paths.empty()) throw std::invalid_argument("You need to specify paths to read files");
    if (kmer_path.empty()) throw std::invalid_argument("You need to specify path to kmers");
    auto kmers_k_pair = load_text_file_kmers(kmer_path);

    SequenceRecordIterator read_iterator = SequenceRecordIterator(read_paths, debug);
    for (auto meta: read_iterator.file_meta){
        std::cout << meta.repr();
    }

    if (output_folder_path.empty()) output_folder_path = fmt::format("./{}_clusters/", read_iterator.meta.filename);

    auto engine = ReadClusteringEngine(read_iterator, config);
    auto cluster_ids = engine.run_clustering(kmers_k_pair.first, kmers_k_pair.second);
    engine.export_components(cluster_ids, output_folder_path);

    return 0;
}
