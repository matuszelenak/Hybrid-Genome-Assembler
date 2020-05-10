#include <fstream>
#include <filesystem>
#include <numeric>
#include "JellyfishOccurrenceReader.h"
#include "../common/Utils.h"

bool queue_compare(std::pair<std::string, std::ifstream*> &x, std::pair<std::string, std::ifstream*> &y){
    return x.first > y.first;
}

CountPair parse_line(std::string &s){
    auto space_pos = s.find(' ');
    auto kmer = s.substr(0, space_pos);
    auto count = std::stoi(s.substr(space_pos, std::string::npos));
    return {kmer, count};
}

std::string merge_sorted_files(std::vector<std::string> &file_paths, std::string output_path){
    std::string line;

    std::priority_queue<
            std::pair<std::string, std::ifstream*>,
            std::vector<std::pair<std::string, std::ifstream*>>,
            std::function<bool(std::pair<std::string, std::ifstream*>&, std::pair<std::string, std::ifstream*>&)>
    > q(queue_compare);

    std::vector<std::ifstream*> files(file_paths.size());
    for (int i = 0; i < file_paths.size(); i++){
        files[i] = new std::ifstream;
        files[i]->open(file_paths[i]);
        getline(*files[i], line);
        q.push({line, files[i]});
    }

    std::ofstream output;
    output.open(output_path);
    while (!q.empty()){
        auto val_file_pair = q.top();
        output << val_file_pair.first << std::endl;
        q.pop();

        if (getline(*val_file_pair.second, line)){
            q.push({line, val_file_pair.second});
        }
    }
    output.close();
    for (auto file_ptr : files){
        file_ptr->close();
        free(file_ptr);
    }
    for (const auto& path : file_paths){
        std::filesystem::remove(path);
    }
    return output_path;
}

JellyfishOccurrenceReader::JellyfishOccurrenceReader(std::vector<std::string> &paths, int k) {
    this->k = k;

    for (auto read_path : paths){
        std::string sorted_path = read_path + "_kmers_sorted";
        if (!std::filesystem::exists(sorted_path)){
            auto partial_sorted_paths = create_sorted_files(read_path);
            merge_sorted_files(partial_sorted_paths, sorted_path);
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

std::vector<std::string> JellyfishOccurrenceReader::create_sorted_files(std::string &read_file_path) {
    std::string jellyfish_cmd = fmt::format("./occurrences/run_jellyfish.sh {} {}", read_file_path, k);
    run_command_with_input(jellyfish_cmd.c_str(), "");

    std::ifstream dumped_kmers;
    dumped_kmers.open("mers.txt");

    std::string kmer_line;
    std::vector<std::string> sorted_files_paths;

    std::vector<std::string> kmer_block;
    int block_counter = 0;

    auto sort_and_dump = [&kmer_block, &block_counter, &read_file_path](){
        std::sort(kmer_block.begin(), kmer_block.end());

        std::ofstream block_file;
        std::string block_file_path = fmt::format("{}_kmers_{}", read_file_path, block_counter);
        block_file.open(block_file_path);
        for (const auto& kmer_count_pair : kmer_block){
            block_file << kmer_count_pair << std::endl;
        }
        block_file.close();
        kmer_block.clear();
        return block_file_path;
    };

    while (getline(dumped_kmers, kmer_line)){
        kmer_block.push_back(kmer_line);

        if (kmer_block.size() == 1000000){
            sorted_files_paths.push_back(sort_and_dump());
            block_counter++;
        }
    }
    if (!kmer_block.empty()) sorted_files_paths.push_back(sort_and_dump());

    std::filesystem::remove("mers.txt");
    return sorted_files_paths;
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

void JellyfishOccurrenceReader::export_kmers(int lower, int upper, std::string &path){
    reset_reader();

    std::ofstream out;
    out.open(path);
    std::string kmer;
    std::vector<int> counts(sorted_paths.size(), 0);
    while (get_next_kmer(kmer, counts)) {
        int total_count = std::accumulate(counts.begin(), counts.end(), 0);
        if (lower <= total_count && total_count <= upper){
            out << fmt::format("{}\n", kmer);
        }
    }
    out.close();
}