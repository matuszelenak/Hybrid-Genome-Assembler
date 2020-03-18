#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <jsoncpp/json/json.h>

#include "structs.h"
#include "SequenceReader.h"
#include "KmerIterator.h"


namespace po = boost::program_options;

using namespace std;


std::ostream &operator<<(std::ostream &os, KmerQuality const &quality) {
    return os << "Min:" << quality.min_quality << "; Avg:" << quality.avg_quality;
}


Json::Value occurrences_to_json(CategoryKmerCounts occurrences) {
    Json::Value root;
    CategoryKmerCounts::iterator iter;
    for (iter = occurrences.begin(); iter != occurrences.end(); iter++) {
        KmerCounts &counts = iter->second;
        KmerCounts::iterator inner;

        Json::Value kmer_counts;
        for (inner = counts.begin(); inner != counts.end(); inner++){
            kmer_counts[to_string(inner->first)] = Json::UInt64(inner->second);
        }
        root[iter->first] = kmer_counts;
    }
    return root;
}


CategoryKmerCounts kmer_occurrences(SequenceReader &reader, int k) {
    CategoryKmerCounts occurrences;

    optional<GenomeRead> read;
    while ((read = reader.get_next_record()) != nullopt) {
        KmerIterator it = KmerIterator(*read, k);
        std::optional<std::pair<Kmer, KmerQuality>> kmer_info;

        if (occurrences.find(read->category) == occurrences.end()){
            occurrences[read->category] = {};
        }

        while ((kmer_info = it.get_next_kmer()) != nullopt) {
            if (occurrences[read->category].find(kmer_info->first) == occurrences[read->category].end()){
                occurrences[read->category][kmer_info->first] = 1;
            } else {
                occurrences[read->category][kmer_info->first] += 1;
            }
        }
    }
    return occurrences;
}


int main(int argc, char *argv[]) {
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_path", 1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_path", po::value<string>(), "Path to file with reads (FASTA or FASTQ)");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc;
        return 0;
    }

    std::string read_path, meta_path, out_path;

    if (vm.count("read_path")) {
        read_path = vm["read_path"].as<string>();
        meta_path = read_path + "_meta.json";
    } else throw invalid_argument("You need to specify a read file");

    Json::Value meta;
    std::ifstream file(meta_path);
    file >> meta;
    cout << meta << endl;

    int k_guess = ceil(log(meta["genome_size"].asUInt64()) / log(4)) + 3;
    cout << k_guess << endl;
    SequenceReader reader = SequenceReader(read_path);
    CategoryKmerCounts occurrences = kmer_occurrences(reader, k_guess);

    std::ofstream file_id;
    file_id.open(read_path + "_occurrences.json");
    Json::FastWriter writer;
    file_id << writer.write(occurrences_to_json(occurrences));
    file_id.close();
    return 0;
}