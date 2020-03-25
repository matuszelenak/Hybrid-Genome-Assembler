#include <vector>
#include <set>
#include <string>
#include <map>
#include <tsl/robin_map.h>

#ifndef SRC_DNASTRUCTURES_H
#define SRC_DNASTRUCTURES_H


typedef uint64_t Kmer;
typedef uint32_t KmerID;
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

typedef tsl::robin_map<Kmer, KmerInfo> KmerOccurrences;

struct InClusterKmerInfo {
    uint32_t sum_of_qualities = 0;
    uint32_t total_count = 0;

    Quality avg_quality(){
        return (Quality)(this->sum_of_qualities / this->total_count);
    };
};

struct KmerQuality {
    Quality min_quality;
    Quality avg_quality;
};


class GenomeReadCluster {
public:
    explicit GenomeReadCluster(ClusterID id, InClusterReadData &read_data, std::set<KmerID> &characteristic_kmers);

    ClusterID reference_id;
    std::vector<InClusterReadData> reads;
    std::set<KmerID> characteristic_kmer_ids;

    void absorb(GenomeReadCluster &cluster);
    std::set<CategoryID> categories;
    std::string consistency();
    uint64_t size();
};


#endif //SRC_DNASTRUCTURES_H
