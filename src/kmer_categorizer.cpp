#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

namespace po = boost::program_options;

using namespace std;

unordered_map<char, uint8_t>BASE_TO_NUM = {
        {'A', 0b00},
        {'C', 0b01},
        {'G', 0b10},
        {'T', 0b11},
};

unordered_map<char, uint8_t> COMPLEMENT = {
        {'A', 0b11},
        {'C', 0b10},
        {'G', 0b01},
        {'T', 0b00}
};

uint8_t BITS_PER_BASE = 2;

uint64_t ALL_SET = 0xFFFFFFFFFFFFFFFF;

struct GenomeRead{
    uint32_t id;
    string sequence;
};

vector<GenomeRead> load_reads(const string &path){
    vector<GenomeRead> result;
    ifstream infile;
    infile.open(path);
    string line;
    uint32_t read_id = 0;
    while(infile >> line){
        infile >> line;
        GenomeRead read = {read_id, line.substr(0, line.size() - 1)};
        result.push_back(read);
    }
    return result;
}

pair<uint64_t, uint64_t> dna_string_to_number(string kmer, int k){
    uint64_t forward_direction = 0;
    uint64_t reverse_direction = 0;
    auto complement_shift_by = (((k - 1) * BITS_PER_BASE));

    uint64_t clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    for (char base : kmer) {
        forward_direction <<= BITS_PER_BASE;
        forward_direction |= BASE_TO_NUM[base];
        forward_direction &= clearing_mask;

        reverse_direction >>= BITS_PER_BASE;
        reverse_direction |= ((uint64_t)COMPLEMENT[base] << complement_shift_by);
    }

    return make_pair(forward_direction, reverse_direction);
}

unordered_map<uint64_t, uint64_t> kmer_occurences(vector<GenomeRead> &mixed_reads, int k){
    unordered_map<uint64_t, uint64_t> kmer_occurences;

    uint64_t clearing_mask = ALL_SET >> (sizeof(uint64_t) * 8 - (k * BITS_PER_BASE));
    auto complement_shift_by = (((k - 1) * BITS_PER_BASE));

    for (auto &read : mixed_reads) {
        pair<uint64_t, uint64_t > initial = dna_string_to_number(read.sequence.substr(0, k - 1), k);
        uint64_t forward_strand_kmer = initial.first;
        uint64_t complementary_strand_kmer = initial.second;

        for (int j = k - 1; j < read.sequence.size(); j++){
            forward_strand_kmer <<= BITS_PER_BASE;
            forward_strand_kmer |= BASE_TO_NUM[read.sequence[j]];
            forward_strand_kmer &= clearing_mask;

            complementary_strand_kmer >>= BITS_PER_BASE;
            complementary_strand_kmer |= ((uint64_t)COMPLEMENT[read.sequence[j]] << complement_shift_by);

            uint64_t kmer_signature = min(forward_strand_kmer, complementary_strand_kmer);

            if (kmer_occurences.find(kmer_signature) == kmer_occurences.end()){
                kmer_occurences[kmer_signature] = 0;
            }
            kmer_occurences[kmer_signature] += 1;
        }
    }

    return kmer_occurences;
}

int main(int argc, char* argv[]) {
    po::options_description desc{"Options"};
    desc.add_options()
            ("help,h", "Help screen")
            ("k", po::value<int >()->default_value(11), "Length of k-mer")
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

    vector<GenomeRead> reads = load_reads(in_path);
    unordered_map<uint64_t, uint64_t> occurences = kmer_occurences(reads, k);

    for (std::pair<uint64_t, uint64_t> element : occurences)
    {
        cout << element.second << ' ';
    }
    return 0;
}