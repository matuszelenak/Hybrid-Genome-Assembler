#include <set>
#include <fmt/format.h>
#include <boost/algorithm/string/join.hpp>

#include "../common/Utils.h"
#include "../common/Types.h"

#ifndef SRC_GENOMEREADCLUSTER_H
#define SRC_GENOMEREADCLUSTER_H

typedef std::pair<uint32_t, bool> Endpoint;
typedef std::pair<uint32_t, uint32_t> Interval;


struct ReadMetaData {
    ReadID id;
    CategoryID category_id;
    uint32_t start;
    uint32_t end;
};


class GenomeReadCluster {
public:
    GenomeReadCluster(ReadMetaData &initial_read_meta, std::vector<KmerID> &discriminative_kmer_ids){
        this->id = initial_read_meta.id;
        this->discriminative_kmer_ids = discriminative_kmer_ids;
        this->contained_reads = {initial_read_meta};
        this->categories = {initial_read_meta.category_id};
    };

    ClusterID id = 0;
    std::vector<ReadMetaData> contained_reads;
    std::vector<KmerID> discriminative_kmer_ids;
    std::set<CategoryID> categories;

    [[nodiscard]] uint64_t size() const{ return contained_reads.size(); };

    std::string consistency(){
        std::map<CategoryID, int> category_counts = {{0, 0}, {1, 0}, {2, 0}};
        for (auto meta : contained_reads) {
            category_counts.insert({meta.category_id, 0}).first->second += 1;
        }
        std::vector<std::string> category_count_vector;
        std::transform(
                category_counts.begin(),
                category_counts.end(),
                std::back_inserter(category_count_vector),
                [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
        );
        return boost::algorithm::join(category_count_vector, "/");
    }

    std::map<CategoryID, std::vector<Interval>> category_intervals(){
        std::map<CategoryID, std::vector<Endpoint>> category_endpoints;
        for (auto meta: contained_reads){
            category_endpoints.insert({meta.category_id, {}}).first->second.emplace_back(meta.start, true);
            category_endpoints.insert({meta.category_id, {}}).first->second.emplace_back(meta.end, false);
        }

        std::map<CategoryID, std::vector<Interval>> result;
        for (auto endpoints : category_endpoints){
            std::sort(endpoints.second.begin(), endpoints.second.end());

            result.insert({endpoints.first, {}});
            uint32_t current_interval_start;
            int opened_intervals = 0;
            for (auto endpoint : endpoints.second) {
                if (endpoint.second) {
                    opened_intervals++;
                    if (opened_intervals == 1) {
                        current_interval_start = endpoint.first;
                    }
                } else {
                    opened_intervals--;
                    if (opened_intervals == 0) {
                        result[endpoints.first].push_back({current_interval_start, endpoint.first});
                    }
                }
            }
        }
        return result;
    }

    std::string to_string(){
        std::vector<std::string> category_intervals_strings;
        for (auto intervals : category_intervals()){
            std::vector<std::string> interval_strings;
            std::transform(
                    intervals.second.begin(),
                    intervals.second.end(),
                    std::back_inserter(interval_strings),
                    [](Interval &i) -> std::string { return fmt::format("({},{})", i.first, i.second); }
            );
            category_intervals_strings.push_back(boost::algorithm::join(interval_strings, ";"));
        }

        return fmt::format("#{} : {} [{}]", this->id, this->consistency(), boost::algorithm::join(category_intervals_strings, " / "));
    }
};

#endif //SRC_GENOMEREADCLUSTER_H
