#ifndef SRC_JELLYFISHOCCURRENCEREADER_H
#define SRC_JELLYFISHOCCURRENCEREADER_H

#include <vector>
#include <string>
#include <queue>
#include <boost/function.hpp>
#include <map>
#include <set>

typedef std::pair<std::string, int> CountPair;
typedef double UpperSpecificity;
typedef int NumOfOccurrences;
typedef int UniqueKmerCount;
typedef std::map<uint16_t, uint64_t> Histogram;
typedef std::map<UpperSpecificity, std::map<NumOfOccurrences, UniqueKmerCount>> KmerSpecificity;

class JellyfishOccurrenceReader {
public:
    int k;
    std::vector<std::string> sorted_paths;
    std::vector<std::ifstream*> sorted_files;
    bool exhausted = false;

    std::string current_kmer;
    std::vector<int> current_kmer_counts;
    std::priority_queue<
            std::pair<CountPair, int>,
            std::vector<std::pair<CountPair, int>>,
            std::function<bool(std::pair<CountPair, int>&, std::pair<CountPair, int>&)>> sorted_kmers_queue;

    JellyfishOccurrenceReader(std::vector<std::string> &paths, int k);

    bool get_next_kmer(std::string &kmer, std::vector<int> &counts);

    void reset_reader();

    KmerSpecificity get_specificity(std::set<double> &thresholds);

    void export_kmers(int lower, int upper, double percent, std::string &path);
};


#endif //SRC_JELLYFISHOCCURRENCEREADER_H
