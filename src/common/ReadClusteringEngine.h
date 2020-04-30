#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <boost/functional/hash.hpp>
#include <experimental/filesystem>
#include <map>
#include <queue>

#include "SequenceRecordIterator.h"
#include "GenomeReadCluster.h"
#include "Utils.h"
#include "../lib/BloomFilter.h"

#ifndef SRC_READCLUSTERINGENGINE_H
#define SRC_READCLUSTERINGENGINE_H

typedef uint32_t KmerID;
typedef tsl::robin_map<Kmer, KmerID> KmerIndex;

typedef uint64_t ConnectionScore;
typedef std::map<ClusterID, GenomeReadCluster*> ClusterIndex;

typedef std::vector<std::vector<ClusterID>> KmerClusterIndex;
typedef tsl::robin_map<KmerID, std::vector<ClusterID>> IndexRemovalMap;
typedef std::vector<ClusterID > IDComponent;

typedef std::vector<std::pair<ClusterID, ClusterID> > EdgeList;
typedef tsl::robin_map<ClusterID, tsl::robin_set<ClusterID> > NeighborMap;


struct ClusterConnection{
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    ConnectionScore score;
    bool is_good;

    bool operator < (const ClusterConnection& conn) const
    {
        return (this->score < conn.score);
    }
    bool operator > (const ClusterConnection& conn) const
    {
        return (this->score > conn.score);
    }
};

struct InterClusterAlignment {
    ClusterID cluster_x_id;
    ClusterID cluster_y_id;
    uint32_t length;
    double identity;

    bool operator < (const InterClusterAlignment& aln) const
    {
        return ((double)this->length * this->identity < (double)aln.length * aln.identity);
    }
    bool operator > (const InterClusterAlignment& aln) const
    {
        return ((double)this->length * this->identity > (double)aln.length * aln.identity);
    }
};


class ReadClusteringEngine {
protected:
    SequenceRecordIterator* reader;
    bloom::BloomFilter<Kmer>* kmers = nullptr;
    tsl::robin_set<Kmer> kmers_set;
    int k;
    Platform platform;

    ClusterIndex cluster_index;
    KmerClusterIndex kmer_cluster_index;
    std::vector<Kmer> kmer_id_to_kmer;
    std::vector<ReadID> ambiguous_reads;

    // TODO remove
    tsl::robin_map<ReadID, std::pair<uint32_t, uint32_t> > read_intervals;

    void construct_indices_thread(KmerIndex &kmer_index);
    int construct_indices();

    void get_connections_thread(ConnectionScore min_score, ConcurrentQueue<ClusterID> &cluster_id_queue, std::vector<ClusterConnection> &accumulator);
    std::vector<ClusterConnection> get_all_connections(ConnectionScore min_score);
    std::vector<ClusterConnection> get_connections(std::vector<ClusterID> &cluster_ids, ConnectionScore min_score);

    void kmer_cluster_index_update(ConcurrentQueue<IndexRemovalMap::value_type> &removal_list_queue);
    void merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal);
    void merge_clusters(std::vector<IDComponent> &components);

    std::vector<std::pair<IDComponent, EdgeList>> union_find(std::vector<ClusterConnection> &connections, std::set<ClusterID> &restricted, int min_component_size);

    std::map<ClusterID, std::vector<std::string>> assemble_clusters(std::vector<ClusterID> &cluster_ids);

    std::vector<InterClusterAlignment> get_alignments(std::map<ClusterID, std::vector<std::string>> &assembly);

    std::vector<ClusterID> filter_clusters(const std::function<bool(GenomeReadCluster*)>& func){
        std::vector<ClusterID> result;
        for (auto cluster_pair : cluster_index){
            if (func(cluster_pair.second)){
                result.push_back(cluster_pair.first);
            }
        }
        return result;
    };

    std::pair<std::vector<ClusterID>, std::vector<ClusterID>> get_spanning_tree_boundary_reads(EdgeList &edges);
public:
    void run_clustering();
    std::map<ClusterID, std::string> export_clusters(std::vector<ClusterID> &cluster_ids, std::experimental::filesystem::path &directory_path);

    ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers, Platform platform);
    ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, tsl::robin_set<Kmer> &kmers, Platform platform);
    ~ReadClusteringEngine();

    void print_clusters(int first_n);
    void print_clusters(std::vector<ClusterID> &ids);

    std::vector<Interval> get_read_coverage_intervals(std::vector<ReadID> &read_ids);

    std::pair<std::vector<KmerID>, std::vector<KmerID>> get_chain_tail_kmers(std::vector<ClusterID> &chain_component, EdgeList &edges);

    std::vector<ClusterConnection> get_chain_connections(std::vector<std::pair<IDComponent, EdgeList>> &components_with_edges);
};

void plot_connection_quality(std::vector<ClusterConnection> &connections);


#endif //SRC_READCLUSTERINGENGINE_H
