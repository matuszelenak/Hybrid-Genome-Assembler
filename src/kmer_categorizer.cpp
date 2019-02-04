#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "SequenceReader.h"
#include "KmerIterator.h"

namespace po = boost::program_options;

using namespace std;

unordered_map<uint64_t, uint64_t> kmer_occurences(SequenceReader &reader, int k){
    unordered_map<uint64_t, uint64_t> result;

    while (true) {
        optional<string> read = reader.get_next_record();
        if (read == nullopt) {
            break;
        }

        KmerIterator it = KmerIterator(*read, k);
        while (true) {
            optional<uint64_t> kmer_signature = it.get_next_kmer();
            if (kmer_signature == nullopt) break;

            if (result.find(*kmer_signature) == result.end()) {
                result[*kmer_signature] = 0;
            }
            result[*kmer_signature] += 1;
        }
    }

    return result;
}

int main(int argc, char* argv[]) {
    po::options_description desc{"Options"};
    desc.add_options()
            ("help,h", "Help screen")
            ("k", po::value<int >()->default_value(11), "Length of k-mer")
            ("cov,c", po::value<int >()->default_value(30), "Coverage of the read file")
            ("in,i", po::value<string >(), "Path to file with reads");

    po::variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    int k = 11;
    string in_path;

    if (vm.count("help")){
        cout << desc;
    }
    if (vm.count("k")){
        k = vm["k"].as<int>();
    }

    if (vm.count("in")){
        in_path = vm["in"].as<string>();
    }

    int c = 30;
    if (vm.count("cov")){
        c = vm["cov"].as<int>();
    }

    SequenceReader reader = SequenceReader(in_path);
    unordered_map<uint64_t, uint64_t> occurences = kmer_occurences(reader, k);

    for (std::pair<uint64_t, uint64_t> element : occurences)
    {
        if (element.second > 4 * c || element.second == 1) continue;
        cout << element.second << ' ';
    }
    return 0;
}