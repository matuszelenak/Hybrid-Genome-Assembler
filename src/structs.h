#include <string>
#include <unordered_map>
#include <vector>

#ifndef SRC_STRUCTS_H
#define SRC_STRUCTS_H

typedef uint64_t Kmer;
typedef std::string Category;
typedef int Quality;
typedef std::unordered_map<Kmer, uint64_t> KmerCounts;
typedef std::unordered_map<Category , KmerCounts> CategoryKmerCounts;

struct GenomeRead {
    std::string header;
    std::string sequence;
    std::vector<Quality > qualities;
    std::string category;
};

struct KmerInfo {
    int in_first_category = 0;
    int in_second_category = 0;
    uint64_t avg_quality = 0;
};

typedef std::unordered_map<Kmer, KmerInfo> KmerOccurrences;

struct KmerQuality {
    Quality min_quality;
    Quality avg_quality;
};

#endif //SRC_STRUCTS_H
