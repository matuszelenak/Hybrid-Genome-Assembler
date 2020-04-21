#include <vector>
#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>

#include "Plotting.h"
#include "Utils.h"

namespace algo = boost::algorithm;

void plot_histogram(Histogram &hist){
    std::vector<std::string> occurrence_strings;
    std::transform(hist.begin(), hist.end(), std::back_inserter(occurrence_strings), [](Histogram::value_type occ) -> std::string {
        return fmt::format("{}: {}", occ.first, occ.second);
    });
    std::string plot_input = fmt::format("{{{}}}", algo::join(occurrence_strings, ", "));
    run_command_with_input("python common/python/plot.py", plot_input);
}


void plot_kmer_specificity(std::map<int, KmerSpecificity> &specificities, int max_coverage) {
    // Lord forgive me for what I am about to do
    std::vector<std::string> k_spec_strings;
    for (const auto &k_specs : specificities) {
        std::vector<std::string> bound_strings;
        for (const auto &bound_specs : k_specs.second) {
            std::vector<std::string> counts;
            for (const auto &coverage_counts : bound_specs.second) {
                if (coverage_counts.second < 50) continue;
                //if (coverage_counts.first < (double)max_coverage / 6.0) continue;
                counts.push_back(fmt::format("({}, {})", coverage_counts.first, coverage_counts.second));
            }
            bound_strings.push_back(fmt::format("({}, [{}])", bound_specs.first, algo::join(counts, ", ")));
        }
        k_spec_strings.push_back(fmt::format("({}, [{}])", k_specs.first, algo::join(bound_strings, ", ")));
    }
    std::string plot_input = fmt::format("{} {}\n{}", specificities.size(), max_coverage, algo::join(k_spec_strings, "\n"));
    std::cout << run_command_with_input("python3 common/python/plot_histogram.py", plot_input) << std::endl;
}