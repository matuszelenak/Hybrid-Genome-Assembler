#include <iostream>
#include <map>
#include <string>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "common/SequenceRecordIterator.h"
#include "common/KmerIterator.h"
#include "common/Utils.h"
#include "common/KmerAnalysis.h"
#include "common/ReadClusteringEngine.h"
#include "lib/BloomFilter.h"


namespace po = boost::program_options;


std::pair<tsl::robin_set<Kmer>, int> load_text_file_kmers(std::string &path){
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
    }
    return {s, k};
}

std::pair<bloom::BloomFilter<Kmer>*, int> load_bloom_filter_kmers(std::string &path){
    auto in = std::ifstream(path, std::ios::in | std::ios::binary);
    int k;
    in.read((char *) &k, sizeof(k));
    return {new bloom::BloomFilter<Kmer>(in), k};
}

int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_paths", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_paths", po::value<std::vector<std::string>>(&read_paths)->multitoken(), "Path to file with reads (FASTA or FASTQ)")
            ("kmers,k", po::value<std::string>(), "Path to text file with kmers")
            ("format,f", po::value<std::string>(), "Specify format of the kmer file [plain|bloom]")
            ("output,o", po::value<std::string>(), "Path to folder with exported clusters");

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
    for (auto meta: read_iterator.file_meta){
        std::cout << meta.repr();
    }

    auto engine = ReadClusteringEngine(read_iterator);

    if (vm.count("kmers")){
        std::string kmer_path = vm["kmers"].as<std::string>();

        if (vm.count("format")){
            std::string format = vm["format"].as<std::string>();
            if (format == "bloom"){
                auto set_k_pair = load_bloom_filter_kmers(kmer_path);
                engine.set_kmers(set_k_pair.first, set_k_pair.second);
            } else if (format == "plain"){
                auto set_k_pair = load_text_file_kmers(kmer_path);
                engine.set_kmers(set_k_pair.first, set_k_pair.second);
            } else {
                throw std::invalid_argument("Invalid kmer file format");
            }
        } else {
            throw std::invalid_argument("Specify the kmer file format");
        }

    } else {
        throw std::invalid_argument("You need to specify path to kmers");
    }

    auto cluster_ids = engine.run_clustering();

    std::experimental::filesystem::path output_folder_path = fmt::format("./{}_clusters/", read_iterator.meta.filename);
    if (vm.count("output")){
        output_folder_path = vm["output"].as<std::string>();
    }
    engine.export_clusters(cluster_ids, output_folder_path);

    return 0;
}
