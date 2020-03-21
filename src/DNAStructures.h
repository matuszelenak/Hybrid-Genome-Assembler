#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <map>

#ifndef SRC_DNASTRUCTURES_H
#define SRC_DNASTRUCTURES_H


typedef uint64_t Kmer;
typedef bool CategoryFlag;
typedef int Quality;

struct GenomeReadData {
    std::string header;
    std::string sequence;
    std::vector<Quality > qualities;
    CategoryFlag category_flag;
};

struct InClusterReadData {
    std::string header;
    CategoryFlag category_flag;
};

struct KmerInfo {
    int in_first_category = 0;
    int in_second_category = 0;
    uint64_t sum_of_qualities = 0;

    int total_occurrences(){
        return this->in_first_category + this->in_second_category;
    };
    Quality avg_quality(){
        return (Quality)(this->sum_of_qualities / this->total_occurrences());
    };
};

typedef std::unordered_map<Kmer, KmerInfo> KmerOccurrences;

struct InClusterKmerInfo {
    uint64_t sum_of_qualities = 0;

    Quality avg_quality(int total_count){
        return (Quality)(this->sum_of_qualities / total_count);
    };
};

struct KmerQuality {
    Quality min_quality;
    Quality avg_quality;
};


class GenomeReadCluster {
private:
    std::set<CategoryFlag> _categories;

public:
    explicit GenomeReadCluster(uint64_t id, InClusterReadData &read_data, std::map<Kmer, InClusterKmerInfo> &characteristic_kmers);

    uint64_t reference_id;
    std::vector<InClusterReadData> reads;
    std::map<Kmer, InClusterKmerInfo> characteristic_kmers;

    void absorb(GenomeReadCluster &cluster);
    std::set<CategoryFlag> categories();
    std::string consistency();
    uint64_t size();
};


#endif //SRC_DNASTRUCTURES_H
