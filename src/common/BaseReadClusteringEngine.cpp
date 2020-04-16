#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <numeric>
#include <experimental/filesystem>
#include <boost/algorithm/string/join.hpp>

#include "KmerIterator.h"
#include "BaseReadClusteringEngine.h"
#include "Utils.h"

namespace fs = std::experimental::filesystem;
namespace algo = boost::algorithm;

std::mutex connections_mut;
std::mutex merge_mut;
std::mutex component_erase_mut;
std::mutex index_merge;

uint64_t MIN_SCORE = 5;


void plot_connection_quality(std::vector<ClusterConnection> &connections){
    std::map<bool, std::map<ConnectionScore , int>>score_to_matching;
    for (auto &conn : connections){
        score_to_matching[conn.is_good].insert(std::map<ConnectionScore , int>::value_type(conn.score, 0)).first->second += 1;
    }
    std::vector<std::string>good_connections, bad_connections;
    for (auto it = begin(score_to_matching[false]); it != end(score_to_matching[false]); it++){
        bad_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    for (auto it = begin(score_to_matching[true]); it != end(score_to_matching[true]); it++){
        good_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    std::string hist_input = fmt::format("{{{}}}\n{{{}}}", algo::join(good_connections, ", "), algo::join(bad_connections, ", "));
    run_command_with_input("python3 common/python/plot_connections_histogram.py", hist_input);
}


std::string BaseReadClusteringEngine::cluster_consistency(GenomeReadCluster *cluster) {
    std::map<CategoryID, int> category_counts;
    for (const auto &read_header : cluster->read_headers) {
        category_counts.insert(std::map<CategoryID, int>::value_type(read_category_map[read_header], 0)).first->second += 1;
    }
    std::vector<std::string> category_count_vector;
    std::transform(
            category_counts.begin(),
            category_counts.end(),
            std::back_inserter(category_count_vector),
            [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
    );
    return algo::join(category_count_vector, "/");
}


void BaseReadClusteringEngine::print_clusters(int first_n) {
    std::vector<GenomeReadCluster *> cluster_pointers;
    std::transform(
            cluster_index.begin(),
            cluster_index.end(),
            std::back_inserter(cluster_pointers),
            [](std::pair<const ClusterID, GenomeReadCluster *> &p) -> GenomeReadCluster * { return p.second; }
    );

    sort(cluster_pointers.rbegin(), cluster_pointers.rend(), [](GenomeReadCluster *x, GenomeReadCluster *y) -> bool { return x->size() < y->size(); });

    int iterate_first = first_n;
    if (first_n == -1) iterate_first = cluster_pointers.size();
    for (auto cluster : cluster_pointers) {
        std::cout << cluster_consistency(cluster) << " ";

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << std::endl;
}

void BaseReadClusteringEngine::construct_read_category_map() {
    this->reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = this->reader->get_next_record()) != std::nullopt) {
        read_category_map.insert(tsl::robin_map<std::string, CategoryID>::value_type(read->header, read->category_id));
    }
}



BaseReadClusteringEngine::BaseReadClusteringEngine(SequenceRecordIterator &read_iterator, KmerCountingBloomFilter &bf, int k, int lower_coverage, int upper_coverage) {
    this->k = k;
    this->reader = &read_iterator;
    this->filter = &bf;
    this->lower_coverage = lower_coverage;
    this->upper_coverage = upper_coverage;

    auto r = timeMeasureMemberFunc(&BaseReadClusteringEngine::construct_indices, this, "Construct indices")();
    construct_read_category_map();
}

void BaseReadClusteringEngine::construct_indices_thread(){
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        std::set<Kmer> in_read_characteristic;

        while (it.next_kmer()) {
            KmerCount count = filter->get_count(it.current_kmer);
            if (lower_coverage <= count && count <= upper_coverage){
                in_read_characteristic.insert(it.current_kmer);
            }
        }

        if (in_read_characteristic.size() >= 50){
            std::set<KmerID> in_read_characteristic_ids;

            index_merge.lock();

            ClusterID new_cluster_id = cluster_index.size();

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (Kmer kmer : in_read_characteristic){
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer, new_kmer_id));
                if (insert_result.second){
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_characteristic_ids.insert(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].insert(new_cluster_id);
            }

            cluster_index.insert(ClusterIndex::value_type(new_cluster_id, new GenomeReadCluster(new_cluster_id, read->header, in_read_characteristic_ids, read->category_id)));

            index_merge.unlock();
        }
    }
}


int BaseReadClusteringEngine::construct_indices() {
    reader->rewind();
    auto runner = ThreadRunner(&BaseReadClusteringEngine::construct_indices_thread, this);
    return 0;
}


void BaseReadClusteringEngine::get_connections_thread(std::vector<ClusterID> &cluster_indices, std::pair<int, int> range, std::vector<ClusterConnection> &accumulator, ProcessedClusters &processed){
    for (int i = range.first; i < range.second; i++){
        ClusterID pivot_id = cluster_indices[i];
        std::map<ClusterID , ConnectionScore > shared_kmer_counts;

        for (KmerID kmer_id : cluster_index[pivot_id]->characteristic_kmer_ids){
            for (ClusterID candidate_id : kmer_cluster_index[kmer_id]){
                shared_kmer_counts.insert( std::map<ClusterID , ConnectionScore >::value_type(candidate_id, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        connections_mut.lock();
        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < MIN_SCORE) continue;

            bool is_newly_processed;
            if (pivot_id < iter->second){
                is_newly_processed = processed.insert(ProcessedClusters::value_type({pivot_id, iter->second})).second;
            } else {
                is_newly_processed = processed.insert(ProcessedClusters::value_type({iter->second, pivot_id})).second;
            }
            if (is_newly_processed){
                accumulator.push_back({iter->first, pivot_id, iter->second, cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories});
            }
        }
        connections_mut.unlock();
    }
}


std::vector<ClusterConnection> BaseReadClusteringEngine::get_connections(){
    std::vector<ClusterConnection> connections;

    std::vector<ClusterID> cluster_ids;
    std::transform(cluster_index.begin(), cluster_index.end(), std::back_inserter(cluster_ids), [](ClusterIndex::value_type &val) -> ClusterID { return val.first; });

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];

    ProcessedClusters processed;
    int clusters_per_thread = (int)ceil(cluster_ids.size() / (double)num_threads);
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(
                &BaseReadClusteringEngine::get_connections_thread,
                this,
                std::ref(cluster_ids),
                std::make_pair(i * clusters_per_thread, std::min((i + 1) * clusters_per_thread, (int)cluster_ids.size())),
                std::ref(connections),
                std::ref(processed)
                );
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}


ClusterID get_parent(ClusterID cluster_id, tsl::robin_map<ClusterID, ClusterID> &parents){
    if (cluster_id == parents[cluster_id]){
        return cluster_id;
    }
    parents[cluster_id] = get_parent(parents[cluster_id], parents);
    return parents[cluster_id];
}


void BaseReadClusteringEngine::merge_clusters_thread(std::queue<IDComponent> &component_queue){
    while (true){
        merge_mut.lock();
        if (component_queue.empty()){
            merge_mut.unlock();
            return;
        }
        IDComponent component = component_queue.front();
        component_queue.pop();
        merge_mut.unlock();

        //TODO faster merging using minheap?
        GenomeReadCluster* survivor = cluster_index[component[0]];
        for (int i = 1; i < component.size(); i++){
            survivor->absorb(*cluster_index[component[i]]);
        }

        component_erase_mut.lock();
        for (int i = 1; i < component.size(); i++){
            for (KmerID kmer_id : cluster_index[component[i]]->characteristic_kmer_ids){
                kmer_cluster_index[kmer_id].erase(component[i]);
                kmer_cluster_index[kmer_id].insert(survivor->reference_id);
            }
            cluster_index.erase(component[i]);
        }
        component_erase_mut.unlock();
    }
}

int BaseReadClusteringEngine::merge_clusters(const tsl::robin_map<ClusterID, IDComponent> &components){
    // Now do the actual merging of clusters in threads
    std::queue<IDComponent> component_queue;
    for (auto component_it = begin(components); component_it != end(components); component_it++){
        if (component_it->second.size() > 1){
            component_queue.push(component_it->second);
        }
    }

    unsigned int num_threads = std::thread::hardware_concurrency();
    std::thread t[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        t[i] = std::thread(&BaseReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue));
    }
    for (int i = 0; i < num_threads; ++i) t[i].join();

    return 0;
}


int BaseReadClusteringEngine::clustering_round(){
    tsl::robin_map<ClusterID, ClusterID>parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++){
        parents[cluster_iter->second->reference_id] = cluster_iter->second->reference_id;
        components[cluster_iter->second->reference_id] = {cluster_iter->second->reference_id};
    }

    std::vector<ClusterConnection> cluster_connections = timeMeasureMemberFunc(&BaseReadClusteringEngine::get_connections, this, "Cluster connections")();
    plot_connection_quality(cluster_connections);
    ConnectionScore min_score;
    std::cin >> min_score;

    int merge_operations = 0;
    for (ClusterConnection &conn : cluster_connections){
        if (conn.score < min_score) continue;
        //if (!conn.is_good) continue;

        ClusterID parent_x = get_parent(conn.cluster_x_id, parents);
        ClusterID parent_y = get_parent(conn.cluster_y_id, parents);
        if (parent_x == parent_y) continue;

        ClusterID bigger, smaller;
        if (components[parent_x].size() > components[parent_y].size()){
            bigger = parent_x;
            smaller = parent_y;
        } else {
            bigger = parent_y;
            smaller = parent_x;
        }

        for (ClusterID cluster_id : components[smaller]){
            parents[cluster_id] = bigger;
        }
        components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
        components.erase(smaller);

        merge_operations++;
    }

    timeMeasureMemberFunc(&BaseReadClusteringEngine::merge_clusters, this, "Merging of clusters")(components);

    return merge_operations;
}

int BaseReadClusteringEngine::export_clusters(int min_size) {
    tsl::robin_map<std::string, GenomeReadData> header_to_read;
    for (auto cluster_it = begin(cluster_index); cluster_it != end(cluster_index); ++cluster_it){
        if (cluster_it->second->size() >= min_size){
            for (const auto& header : cluster_it->second->read_headers){
                header_to_read[header] = {};
            }
        }
    }

    auto path_concat = [](std::string acc, ReadFileMetaData &m){
        return std::move(acc) + std::string("__") +  std::move(m.filename);
    };

    fs::path directory_path = "./clusters/" + std::accumulate(std::next(reader->file_meta.begin()), reader->file_meta.end(), reader->file_meta[0].filename, path_concat);
    fs::remove_all(directory_path);
    fs::create_directories(directory_path);

    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        if (header_to_read.contains(read->header)){
            header_to_read[read->header] = *read;
        }
    }

    int exported_clusters = 0;
    for (auto cluster_it = begin(cluster_index); cluster_it != end(cluster_index); ++cluster_it){
        if (cluster_it->second->size() >= min_size){
            std::ofstream cluster_file;
            cluster_file.open(fmt::format("{}/{}_reads:cluster#{}.fq", directory_path.string(), cluster_it->second->size(), cluster_it->second->reference_id));
            for (const auto& header: cluster_it->second->read_headers){
                cluster_file << header_to_read[header].fastq_string() << std::endl;
            }
            exported_clusters++;
            cluster_file.close();
        }
    }

    return exported_clusters;
}

void BaseReadClusteringEngine::run_clustering(){
    while (clustering_round() > 0){
        print_clusters(100);
    }
    print_clusters(-1);
    std::cout << fmt::format("Exported {} clusters\n", timeMeasureMemberFunc(&BaseReadClusteringEngine::export_clusters, this, "Exporting clusters")(10));
}