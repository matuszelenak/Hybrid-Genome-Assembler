#include <fstream>
#include <filesystem>
#include <numeric>
#include <random>
#include "JellyfishOccurrenceReader.h"
#include "../common/Utils.h"


CountPair parse_line(std::string &s){
    auto space_pos = s.find(' ');
    auto kmer = s.substr(0, space_pos);
    auto count = std::stoi(s.substr(space_pos, std::string::npos));
    return {kmer, count};
}

JellyfishOccurrenceReader::JellyfishOccurrenceReader(std::vector<std::string> &paths, int k) {
    this->k = k;

    for (auto read_path : paths){
        std::string sorted_path = fmt::format("{}_{}-mers_sorted", read_path, k);
        if (!std::filesystem::exists(sorted_path)){
            std::string jellyfish_cmd = fmt::format("./occurrences/run_jellyfish.sh {} {} {}", read_path, k, sorted_path);
            run_command_with_input(jellyfish_cmd.c_str(), "");
        }
        sorted_paths.push_back(sorted_path);
    }

    auto q_compare = [](std::pair<CountPair, int>&x, std::pair<CountPair, int>&y){
        return x.first.first > y.first.first;
    };

    sorted_kmers_queue = std::priority_queue<
            std::pair<CountPair, int>,
            std::vector<std::pair<CountPair, int>>,
            std::function<bool(std::pair<CountPair, int>&, std::pair<CountPair, int>&)>>(q_compare);

    reset_reader();
}

void JellyfishOccurrenceReader::reset_reader(){
    while (!sorted_kmers_queue.empty()) sorted_kmers_queue.pop();

    std::string line;
    sorted_files.resize(sorted_paths.size());
    for (int i = 0; i < sorted_files.size(); i++){
        if (sorted_files[i] != nullptr){
            sorted_files[i]->close();
            free(sorted_files[i]);
        }
        sorted_files[i] = new std::ifstream;
        sorted_files[i]->open(sorted_paths[i]);

        getline(*sorted_files[i], line);
        auto parsed = parse_line(line);
        sorted_kmers_queue.push({parsed, i});
    }

    current_kmer_counts.resize(sorted_paths.size(), 0);
    current_kmer = sorted_kmers_queue.top().first.first;
    exhausted = false;
}

bool JellyfishOccurrenceReader::get_next_kmer(std::string &kmer, std::vector<int> &counts){
    if (exhausted) return false;

    std::fill(counts.begin(), counts.end(), 0);
    while (!sorted_kmers_queue.empty() && sorted_kmers_queue.top().first.first == current_kmer){
        auto top = sorted_kmers_queue.top();
        counts[top.second] = top.first.second;
        sorted_kmers_queue.pop();

        std::string line;
        if (getline(*sorted_files[top.second], line)){
            auto parsed = parse_line(line);
            sorted_kmers_queue.push({parsed, top.second});
        }
    }
    kmer = current_kmer;

    if (sorted_kmers_queue.empty()){
        exhausted = true;
    } else {
        current_kmer = sorted_kmers_queue.top().first.first;
    }
    return true;
}

KmerSpecificity JellyfishOccurrenceReader::get_specificity(std::set<double> &thresholds){
    KmerSpecificity result;
    for (auto threshold : thresholds) {
        result.insert({threshold, {}});
    }

    std::string kmer;
    std::vector<int> counts(sorted_paths.size(), 0);
    while (get_next_kmer(kmer, counts)){
        int prevalent_count = 0;
        int total_count = 0;
        for (auto count : counts){
            prevalent_count = std::max(count, prevalent_count);
            total_count += count;
        }
        UpperSpecificity kmer_specificity = *thresholds.upper_bound(((double) prevalent_count / (double) total_count) * 100);
        result[kmer_specificity].insert(std::pair<NumOfOccurrences, UniqueKmerCount>(total_count, 0)).first->second += 1;
    }

    return result;
}

void JellyfishOccurrenceReader::export_kmers(int lower, int upper, double percent, std::string &path){
    reset_reader();
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::ofstream out;
    out.open(path);
    std::string kmer;
    std::vector<int> counts(sorted_paths.size(), 0);
    while (get_next_kmer(kmer, counts)) {
        int total_count = std::accumulate(counts.begin(), counts.end(), 0);
        if (lower <= total_count && total_count <= upper && dis(rng) < percent){
            out << fmt::format("{}\n", kmer);
        }
        if (std::count_if(counts.begin(), counts.end(), [](int c){ return c > 0; }) == 1){
            //out << fmt::format("{}\n", kmer);
        }
    }
    out.close();
}