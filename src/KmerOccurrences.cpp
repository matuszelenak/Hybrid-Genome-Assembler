#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/algorithm/string/join.hpp>
#include <jsoncpp/json/json.h>
#include <fmt/format.h>

#include "structs.h"
#include "SequenceReader.h"
#include "KmerIterator.h"
#include "Utils.h"


namespace po = boost::program_options;
namespace algo = boost::algorithm;

typedef int UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;


std::ostream &operator<<(std::ostream &os, KmerQuality const &quality) {
    return os << "Min:" << quality.min_quality << "; Avg:" << quality.avg_quality;
}

void plot_kmer_specificity(std::map<int, KmerSpecificity> &specificities) {
    // Lord forgive me for what I am about to do
    std::vector<std::string> k_spec_strings;
    for (const auto &k_specs : specificities) {
        std::vector<std::string> bound_strings;
        for (const auto &bound_specs : k_specs.second) {
            std::vector<std::string> counts;
            for (const auto &coverage_counts : bound_specs.second) {
                counts.push_back(fmt::format("({}, {})", coverage_counts.first, coverage_counts.second));
            }
            bound_strings.push_back(fmt::format("({}, [{}])", bound_specs.first, algo::join(counts, ", ")));
        }
        k_spec_strings.push_back(fmt::format("({}, [{}])", k_specs.first, algo::join(bound_strings, ", ")));
    }
    std::string plot_input = fmt::format("{}\n{}", specificities.size(), algo::join(k_spec_strings, "\n"));
    std::cout << run_command_with_input("python python_scripts/plot_histogram.py", plot_input) << std::endl;
}

KmerSpecificity get_kmer_specificity(KmerOccurrences &occurrences) {
    KmerSpecificity specificity;

    std::set<UpperSpecificity> upper_specificity_bounds = {70, 85, 90, 95, 99, 101};
    for (auto s : upper_specificity_bounds) {
        specificity[s] = {};
    }

    KmerOccurrences::iterator iter;
    for (iter = occurrences.begin(); iter != occurrences.end(); iter++) {
        uint64_t total_count = iter->second.in_first_category + iter->second.in_second_category;
        uint64_t prevalent = std::max(iter->second.in_first_category, iter->second.in_second_category);
        UpperSpecificity kmer_specificity = *upper_specificity_bounds.upper_bound(((double) prevalent / (double) total_count) * 100);
        if (!specificity[kmer_specificity].contains(total_count)) {
            specificity[kmer_specificity][total_count] = 0;
        }
        specificity[kmer_specificity][total_count] += 1;
    }
    return specificity;
}

KmerOccurrences kmer_occurrences(SequenceReader &reader, int k) {
    KmerOccurrences occurrences;

    std::optional<GenomeRead> read;
    std::map<Category, bool> categories;
    while ((read = reader.get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(*read, k);
        std::optional<std::pair<Kmer, KmerQuality>> kmer_info;

        if (!categories.contains(read->category)) {
            categories[read->category] = categories.empty();
            if (categories.size() > 2) {
                throw std::logic_error("More than two categories found in reads");
            }
        }

        while ((kmer_info = it.get_next_kmer()) != std::nullopt) {
            if (!occurrences.contains(kmer_info->first)) {
                occurrences[kmer_info->first] = {0, 0, (uint64_t) kmer_info->second.avg_quality};
            }
            occurrences[kmer_info->first].in_first_category += (int) categories[read->category];
            occurrences[kmer_info->first].in_second_category += (int) !categories[read->category];
        }
    }

    KmerOccurrences::iterator it;
    for (it = occurrences.begin(); it != occurrences.end(); it++) {
        it->second.avg_quality = it->second.avg_quality / (it->second.in_first_category + it->second.in_second_category);
    }
    return occurrences;
}


int main(int argc, char *argv[]) {
    po::options_description desc{"Options"};
    po::positional_options_description p;
    p.add("read_path", 1);
    desc.add_options()
            ("help,h", "Help screen")
            ("read_path", po::value<std::string>(), "Path to file with reads (FASTA or FASTQ)");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc;
        return 0;
    }

    std::string read_path, meta_path, out_path;

    if (vm.count("read_path")) {
        read_path = vm["read_path"].as<std::string>();
        meta_path = read_path + "_meta.json";
    } else throw std::invalid_argument("You need to specify a read file");

    Json::Value meta;
    std::ifstream file(meta_path);
    file >> meta;
    std::cout << meta << std::endl;

    int k_guess = ceil(log(meta["genome_size"].asUInt64()) / log(4)) + 3;

    std::map<int, KmerSpecificity> specificities = {};
    std::map<int, KmerOccurrences> per_k_occurrences;
    for (int k_length = k_guess; k_length < k_guess + 2; k_length++) {
        std::cout << k_length << std::endl;

        SequenceReader reader = SequenceReader(read_path);
        KmerOccurrences occurrences = kmer_occurrences(reader, k_length);
        specificities[k_length] = get_kmer_specificity(occurrences);
        per_k_occurrences[k_length] = occurrences;
    }

    plot_kmer_specificity(specificities);
    int selected_k, cov_low, cov_high;
    //std::cin >> selected_k >> cov_low >> cov_high;

    return 0;
}