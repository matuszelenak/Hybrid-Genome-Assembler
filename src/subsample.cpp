#include <iostream>
#include <map>
#include <string>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <fmt/format.h>
#include <random>

#include "common/SequenceRecordIterator.h"


namespace po = boost::program_options;


int main(int argc, char *argv[]) {
    std::vector<std::string> read_paths;

    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("input", -1);
    desc.add_options()
            ("help,h", "Help screen")
            ("input", po::value<std::string> (), "Output path for the counting bloom filter");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    std::string s = vm["input"].as<std::string>();
    auto reader = SequenceRecordIterator(s);

    std::cout << "Enter the approximate number of reads you want\n";
    int wanted_reads;
    std::cin >> wanted_reads;

    double percentage = (double) wanted_reads / (double) reader.meta.records;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::ofstream subsampled_file;
    std::string cluster_file_path = fmt::format("{}_subsampled", s);

    subsampled_file.open(cluster_file_path);

    std::optional<GenomeReadData> read;
    while ((read = reader.get_next_record()) != std::nullopt){
        if (dis(rng) < percentage){
            subsampled_file << read->fastX_string() << std::endl;
        }
    }
    subsampled_file.close();
}