#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <experimental/filesystem>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>

#include "SequenceReader.h"
#include "KmerIterator.h"

namespace po = boost::program_options;

using namespace std;

unordered_map<uint64_t , uint64_t > load_kmer_counts(const string &path){
    ifstream infile;
    infile.open(path);
    uint64_t kmer, count;
    unordered_map<uint64_t , uint64_t > result;
    while (infile >> kmer >> count){
        result[kmer] = count;
    }
    return result;
}

void characteristic_kmer_positions(const string &path, unordered_map<uint64_t, uint64_t> &kmers, int k, const string &out_path){
    SequenceReader reader = SequenceReader(path);
    ofstream out_file;
    out_file.open(out_path);
    while (true){
        optional<string> read = reader.get_next_record();
        if (read == nullopt) {
            break;
        }

        KmerIterator it = KmerIterator(*read, k);
        unsigned long position = 0;
        while (true) {
            optional<uint64_t> kmer_signature = it.get_next_kmer();
            if (kmer_signature == nullopt) break;

            if (kmers.find(*kmer_signature) != kmers.end()){
                out_file << position << " ";
            }
            position ++;
        }
        out_file << endl;
    }
    out_file.close();
}

void kmer_counts(SequenceReader &reader, int k, const string &out_path){
    unordered_map<uint64_t, uint64_t> counts;

    ofstream out_file;
    out_file.open(out_path);
    while (true) {
        optional<string> read = reader.get_next_record();
        if (read == nullopt) {
            break;
        }

        KmerIterator it = KmerIterator(*read, k);
        while (true) {
            optional<uint64_t> kmer_signature = it.get_next_kmer();
            if (kmer_signature == nullopt) break;

            if (counts.find(*kmer_signature) == counts.end()) {
                counts[*kmer_signature] = 0;
            }
            counts[*kmer_signature] += 1;
        };
    }
    for (std::pair<uint64_t, uint64_t> element : counts)
    {
        out_file << element.first << ' ' << element.second << endl;
    }
    out_file.close();
}

int main(int argc, char* argv[]) {
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("command", 1).add("read_path", 1);
    desc.add_options()
            ("help,h", "Help screen")
            ("k", po::value<int >()->default_value(11), "Length of k-mer")
            ("cov,c", po::value<int >()->default_value(30), "Coverage of the read file")
            ("command", po::value<string >(), "Command")
            ("read_path", po::value<string >(), "Path to file with reads (FASTA or FASTQ)")
            ("counts", po::value<string >(), "Path to file with computed kmer counts")
            ("out,o", po::value<string >(), "Path to output file")
            ("low,l", po::value<int >(), "Lower bound for coverage of characteristic kmers")
            ("upper,u", po::value<int >(), "Upper bound for coverage of characteristic kmers");


    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")){
        cout << desc;
        return 0;
    }

    std::experimental::filesystem::path read_path, out_path;
    string command;
    int k = 11;
    int lower = 1;
    int upper = 1000;

    if (vm.count("command")){
        command = vm["command"].as<string>();
    } else throw invalid_argument("You need to specify a command");

    if (vm.count("read_path")){
        read_path = vm["read_path"].as<string>();
    } else throw invalid_argument("You need to specify a read file");

    std::experimental::filesystem::path path(read_path);

    if (vm.count("k")){
        k = vm["k"].as<int>();
    }

    if (vm.count("out")){
        out_path = vm["out"].as<string>();
    }

    if (command == "count"){
        int c = 30;
        if (vm.count("cov")){
            c = vm["cov"].as<int>();
        }

        if (out_path.empty()){
            out_path = path.replace_extension("counts");
        }

        SequenceReader reader = SequenceReader(read_path);
        kmer_counts(reader, k, out_path);
    }
    else if (command == "analyze"){

        if (vm.count("lower")){
            lower = vm["lower"].as<int>();
        }

        if (vm.count("upper")){
            upper = vm["upper"].as<int>();
        }

        string counts;
        if (vm.count("counts")){
            counts = vm["counts"].as<string>();
        } else throw invalid_argument("You need to specify file with computed counts");

        if (out_path.empty()){
            out_path = path.replace_extension("analysis");
        }

        unordered_map<uint64_t, uint64_t> kmer_counts = load_kmer_counts(counts);
        for(auto it = begin(kmer_counts); it != end(kmer_counts);)
        {
            if (it->second < lower || it-> second > upper)
            {
                it = kmer_counts.erase(it);
            }
            else
                ++it;
        }
        characteristic_kmer_positions(read_path, kmer_counts, k, out_path);
    }
//

    return 0;
}