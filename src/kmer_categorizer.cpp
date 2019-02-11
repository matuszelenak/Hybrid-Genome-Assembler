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
#include "KmerAnalysisWriter.h"
#include "ReadCategorizer.h"

namespace po = boost::program_options;

using namespace std;


set<uint64_t > load_kmer_counts(const string &path, int lower, int upper){
    set<uint64_t > result;
    SequenceReader reader = SequenceReader(path);
    optional<GenomeRead> read;
    while ((read = reader.get_next_record()) != nullopt){
        pair<uint64_t , uint64_t > nums = KmerIterator::sequence_to_number((*read).sequence);
        int kmer_count = stoi((*read).header);
        if (kmer_count >= lower && kmer_count <= upper){
            result.insert(min(nums.first, nums.second));
        }
    }
    return result;
}

void get_kmer_positions(SequenceReader &reader, KmerAnalysisWriter &writer, std::set<uint64_t>characteristic_kmers, int k){
    optional<GenomeRead> read;
    while ((read = reader.get_next_record()) != nullopt){
        vector<pair<uint64_t, unsigned long> > positions;
        unsigned long position = 0;

        KmerIterator it = KmerIterator(*read, k);
        optional<uint64_t> kmer_signature;

        while ((kmer_signature = it.get_next_kmer()) != nullopt) {
            if (characteristic_kmers.find(*kmer_signature) != characteristic_kmers.end()){
                positions.emplace_back(make_pair(*kmer_signature, position));
            }
            position ++;
        }
        writer.write_sequence_data((*read).header, positions);
    }
}

void characteristic_kmer_positions(const string &path, set<uint64_t> &kmers, const string &out_path, int k){
    SequenceReader reader = SequenceReader(path);
    KmerAnalysisWriter writer = KmerAnalysisWriter(out_path);

    int num_threads = thread::hardware_concurrency();

    thread t[num_threads];

    for (int i = 0; i < num_threads; ++i) {
        t[i] = thread(get_kmer_positions, ref(reader), ref(writer), ref(kmers), k);
    }

    for (int i = 0; i < num_threads; ++i) {
        t[i].join();
    }
}

void categorize_reads(const string &path, set<uint64_t> &kmers, int k, int num_categories){
    ReadCategorizer categorizer = ReadCategorizer(path, kmers, k, num_categories);
}

void kmer_counts(SequenceReader &reader, int k, const string &out_path){
    unordered_map<uint64_t, uint64_t> counts;

    ofstream out_file;
    out_file.open(out_path);
    optional<GenomeRead> read;
    while ((read = reader.get_next_record()) != nullopt){

        KmerIterator it = KmerIterator(*read, k);
        optional<uint64_t> kmer_signature;

        while ((kmer_signature = it.get_next_kmer()) != nullopt) {
            if (counts.find(*kmer_signature) == counts.end()) {
                counts[*kmer_signature] = 0;
            }
            counts[*kmer_signature] += 1;
        };
    }
    for (std::pair<uint64_t, uint64_t> element : counts)
    {
        out_file << ">" << element.second << endl;
        out_file << KmerIterator::number_to_sequence(element.first, k) << endl;
    }
    out_file.close();
}

int main(int argc, char* argv[]) {
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("command", 1).add("read_path", 1);
    desc.add_options()
            ("help,h", "Help screen")
            ("kmer,k", po::value<int >()->default_value(11), "Length of k-mer")
            ("categories,c", po::value<int >()->default_value(2), "Number of desired read categories")
            ("command", po::value<string >(), "Command")
            ("read_path", po::value<string >(), "Path to file with reads (FASTA or FASTQ)")
            ("counts", po::value<string >(), "Path to file with computed kmer counts")
            ("out,o", po::value<string >(), "Path to output file")
            ("lower,l", po::value<int >(), "Lower bound for coverage of characteristic kmers")
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

    if (vm.count("kmer")){
        k = vm["kmer"].as<int>();
    }

    if (vm.count("out")){
        out_path = vm["out"].as<string>();
    }

    if (command == "count"){
        if (out_path.empty()){
            out_path = path.replace_extension("counts");
        }

        SequenceReader reader = SequenceReader(read_path);
        kmer_counts(reader, k, out_path);
    }
    else if (command == "analyze" || command == "categorize"){

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

        set<uint64_t > characteristic_kmers = load_kmer_counts(counts, lower, upper);
        if (command == "analyze")
            characteristic_kmer_positions(read_path, characteristic_kmers, out_path, k);
        if (command == "categorize"){
            int cats;
            if (vm.count("categories")){
                cats = vm["categories"].as<int>();
            }
            categorize_reads(read_path, characteristic_kmers, k, cats);
        }
    }
    return 0;
}