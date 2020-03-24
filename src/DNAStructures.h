#include <vector>
#include <set>
#include <string>
#include <map>
#include <tsl/robin_map.h>

#ifndef SRC_DNASTRUCTURES_H
#define SRC_DNASTRUCTURES_H


typedef uint64_t Kmer;
typedef int CategoryID;
typedef int Quality;
typedef uint32_t ClusterID;

struct GenomeReadData {
    std::string header;
    std::string sequence;
    std::vector<Quality > qualities;
    CategoryID category_id;
};

struct InClusterReadData {
    std::string header;
    CategoryID category_id;
};

struct KmerInfo {
    uint16_t in_first_category = 0;
    uint16_t in_second_category = 0;
    uint32_t sum_of_qualities = 0;

    [[nodiscard]] uint16_t total_occurrences() const {
        return this->in_first_category + this->in_second_category;
    };
    Quality avg_quality(){
        return (Quality)(this->sum_of_qualities / this->total_occurrences());
    };
};

//typedef tsl::robin_map<Kmer, KmerInfo, std::hash<Kmer>, std::equal_to<>, std::allocator<std::pair<Kmer, KmerInfo>>, true> KmerOccurrences;
typedef tsl::robin_map<Kmer, KmerInfo> KmerOccurrences;

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
public:
    explicit GenomeReadCluster(ClusterID id, InClusterReadData &read_data, std::map<Kmer, InClusterKmerInfo> &characteristic_kmers);

    ClusterID reference_id;
    std::vector<InClusterReadData> reads;
    std::map<Kmer, InClusterKmerInfo> characteristic_kmers;

    void absorb(GenomeReadCluster &cluster);
    std::set<CategoryID> categories;
    std::string consistency();
    uint64_t size();
};


#endif //SRC_DNASTRUCTURES_H
