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
std::mutex conn_index_mut;

uint64_t MIN_SCORE = 1;


void plot_connection_quality(std::vector<ClusterConnection> &connections) {
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching;
    for (auto &conn : connections) {
        score_to_matching[conn.is_good].insert(std::map<ConnectionScore, int>::value_type(conn.score, 0)).first->second += 1;
    }
    std::vector<std::string> good_connections, bad_connections;
    for (auto it = begin(score_to_matching[false]); it != end(score_to_matching[false]); it++) {
        bad_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    for (auto it = begin(score_to_matching[true]); it != end(score_to_matching[true]); it++) {
        good_connections.push_back(fmt::format("{}: {}", it->first, it->second));
    }
    std::string hist_input = fmt::format("{{{}}}\n{{{}}}", algo::join(good_connections, ", "), algo::join(bad_connections, ", "));
    run_command_with_input("python scripts/plotting.py --plot connection_histogram", hist_input);
}


void plot_cluster_coverage(std::vector<GenomeReadCluster*> &clusters){
    std::map<CategoryID, std::vector<std::string> > category_intervals = {{0, {}}, {1, {}}};
    for (auto cluster : clusters){
        for (int i = 0; i < cluster->endpoints.size(); i += 2){
            category_intervals[*cluster->categories.begin()].push_back(fmt::format("({}, {})", cluster->endpoints[i].first, cluster->endpoints[i + 1].first));
        }
    }
    std::string input;
    for (const auto& k_v : category_intervals){
        input += ("[" + algo::join(k_v.second, ", ") + "]\n");
    }
    run_command_with_input("python scripts/plotting.py --plot cluster_coverage", input);
}


std::string BaseReadClusteringEngine::cluster_consistency(GenomeReadCluster *cluster) {
    std::map<CategoryID, int> category_counts = {{0, 0}, {1, 0}};
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

    int size_one_clusters = 0;
    for (auto cluster : cluster_pointers) {
        if (cluster->size() > 1){
            std::cout << cluster_consistency(cluster) << " " << cluster->intervals() << std::endl;
        } else {
            size_one_clusters++;
        }

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << std::endl;
    std::cout << fmt::format("As well as {} clusters of size 1\n", size_one_clusters);
}

void BaseReadClusteringEngine::construct_read_category_map() {
    this->reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = this->reader->get_next_record()) != std::nullopt) {
        read_category_map.insert(tsl::robin_map<std::string, CategoryID>::value_type(read->header, read->category_id));
    }
}

BaseReadClusteringEngine::BaseReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers) {
    this->k = k;
    this->reader = &read_iterator;
    this->kmers = &kmers;

    timeMeasureMemberFunc(&BaseReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", cluster_index.size(), reader->meta.records);

    construct_read_category_map();
}

void BaseReadClusteringEngine::construct_indices_thread() {
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        tsl::robin_set<Kmer> in_read_discriminative;
        while (it.next_kmer()) {
            if (kmers->contains(it.current_kmer)) {
                in_read_discriminative.insert(it.current_kmer);
            }
        }

        if (in_read_discriminative.size() >= 5) {
            std::vector<KmerID> in_read_discriminative_ids;

            index_merge.lock();

            ClusterID new_cluster_id = cluster_index.size();

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (Kmer kmer : in_read_discriminative) {
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer, new_kmer_id));
                if (insert_result.second) {
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_discriminative_ids.push_back(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].insert(new_cluster_id);
            }

            std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
            cluster_index.insert(ClusterIndex::value_type(new_cluster_id, new GenomeReadCluster(new_cluster_id, *read, in_read_discriminative_ids, read->category_id)));

            index_merge.unlock();
        }
    }
}


int BaseReadClusteringEngine::construct_indices() {
    reader->rewind();
    auto runner = ThreadRunner(&BaseReadClusteringEngine::construct_indices_thread, this);
    return 0;
}


void BaseReadClusteringEngine::get_connections_thread(std::vector<ClusterID> &cluster_indices, std::vector<ClusterConnection> &accumulator, int &index) {
    uint64_t total_calculated = 0;
    while (true){
        conn_index_mut.lock();
        int pivot_index = index++;
        conn_index_mut.unlock();

        if (pivot_index >= cluster_indices.size()) break;
        ClusterID pivot_id = cluster_indices[pivot_index];

        std::map<ClusterID, ConnectionScore> shared_kmer_counts;

        for (KmerID kmer_id : cluster_index[pivot_id]->discriminative_kmer_ids) {
            for (auto candidate_it = kmer_cluster_index[kmer_id].begin(); candidate_it != kmer_cluster_index[kmer_id].lower_bound(pivot_id); ++candidate_it){
                shared_kmer_counts.insert(std::map<ClusterID, ConnectionScore>::value_type(*candidate_it, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        total_calculated += shared_kmer_counts.size();

        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < MIN_SCORE) continue;
            connections_mut.lock();
            accumulator.push_back({iter->first, pivot_id, iter->second, cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories});
            connections_mut.unlock();
        }
    }
    std::cout << total_calculated << std::endl;
}

std::vector<ClusterConnection> BaseReadClusteringEngine::get_connections(std::vector<ClusterID> &cluster_ids) {
    std::vector<ClusterConnection> connections;

    int index = 0;
    auto r = ThreadRunner(&BaseReadClusteringEngine::get_connections_thread, this, std::ref(cluster_ids), std::ref(connections), std::ref(index));

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

ClusterID get_parent(ClusterID cluster_id, tsl::robin_map<ClusterID, ClusterID> &parents) {
    if (cluster_id == parents[cluster_id]) {
        return cluster_id;
    }
    parents[cluster_id] = get_parent(parents[cluster_id], parents);
    return parents[cluster_id];
}

void BaseReadClusteringEngine::merge_clusters_thread(std::queue<IDComponent> &component_queue) {
    while (true) {
        merge_mut.lock();
        if (component_queue.empty()) {
            merge_mut.unlock();
            return;
        }
        IDComponent component = component_queue.front();
        component_queue.pop();
        merge_mut.unlock();

        GenomeReadCluster *survivor = cluster_index[component[0]];

        std::vector<std::vector<KmerID> *> arrays_to_merge = {&survivor->discriminative_kmer_ids};
        std::vector<std::vector<Endpoint>* > endpoints_to_merge = {&survivor->endpoints};
        for (int i = 1; i < component.size(); i++) {
            survivor->read_headers.insert(survivor->read_headers.end(), cluster_index[component[i]]->read_headers.begin(), cluster_index[component[i]]->read_headers.end());
            survivor->categories.insert(cluster_index[component[i]]->categories.begin(), cluster_index[component[i]]->categories.end());

            arrays_to_merge.push_back(&cluster_index[component[i]]->discriminative_kmer_ids);
            endpoints_to_merge.push_back(&cluster_index[component[i]]->endpoints);
        }

        std::vector<KmerID> merged_discriminative_ids;
        merge_n_vectors(arrays_to_merge, merged_discriminative_ids, true);

        std::vector<Endpoint> merged_endpoints;
        merge_n_vectors(endpoints_to_merge, merged_endpoints, false);

        std::vector<Endpoint> final_endpoints;
        int opened_intervals = 0;
        for (auto endpoint : merged_endpoints){
            if (endpoint.second){
                opened_intervals++;
                if (opened_intervals == 1){
                    final_endpoints.push_back(endpoint);
                }
            } else {
                opened_intervals--;
                if (opened_intervals == 0){
                    final_endpoints.push_back(endpoint);
                }
            }
        }

        survivor->endpoints = final_endpoints;
        survivor->discriminative_kmer_ids = merged_discriminative_ids;

        component_erase_mut.lock();
        for (int i = 1; i < component.size(); i++) {
            for (KmerID kmer_id : cluster_index[component[i]]->discriminative_kmer_ids) {
                kmer_cluster_index[kmer_id].erase(component[i]);
                kmer_cluster_index[kmer_id].insert(survivor->reference_id);
            }
            cluster_index.erase(component[i]);
        }
        component_erase_mut.unlock();
    }
}

int BaseReadClusteringEngine::merge_clusters(const tsl::robin_map<ClusterID, IDComponent> &components) {
    std::queue<IDComponent> component_queue;
    for (auto component_it = begin(components); component_it != end(components); component_it++) {
        if (component_it->second.size() > 1) {
            component_queue.push(component_it->second);
        }
    }

    auto r = ThreadRunner(&BaseReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue));
    return 0;
}

int BaseReadClusteringEngine::clustering_round() {
    tsl::robin_map<ClusterID, ClusterID> parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++) {
        parents[cluster_iter->second->reference_id] = cluster_iter->second->reference_id;
        components[cluster_iter->second->reference_id] = {cluster_iter->second->reference_id};
    }

    std::vector<ClusterID> cluster_ids;
    std::transform(cluster_index.begin(), cluster_index.end(), std::back_inserter(cluster_ids), [](ClusterIndex::value_type &val) -> ClusterID { return val.first; });

    std::vector<ClusterConnection> cluster_connections = timeMeasureMemberFunc(&BaseReadClusteringEngine::get_connections, this, "Cluster connections")(cluster_ids);
    plot_connection_quality(cluster_connections);

    ConnectionScore min_score;
    std::cout << "Enter minimal score\n";
    std::cin >> min_score;

    int merge_operations = 0;
    for (ClusterConnection &conn : cluster_connections) {
        if (conn.score < min_score) continue;
        //if (!conn.is_good) continue;

        ClusterID parent_x = get_parent(conn.cluster_x_id, parents);
        ClusterID parent_y = get_parent(conn.cluster_y_id, parents);
        if (parent_x == parent_y) continue;

        ClusterID bigger, smaller;
        if (components[parent_x].size() > components[parent_y].size()) {
            bigger = parent_x;
            smaller = parent_y;
        } else {
            bigger = parent_y;
            smaller = parent_x;
        }

        for (ClusterID cluster_id : components[smaller]) {
            parents[cluster_id] = bigger;
        }
        components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
        components.erase(smaller);

        merge_operations++;
    }

    timeMeasureMemberFunc(&BaseReadClusteringEngine::merge_clusters, this, "Merging of clusters")(components);

    return merge_operations;
}

std::map<ClusterID, std::string> BaseReadClusteringEngine::export_clusters(std::vector<ClusterID> &cluster_ids, fs::path &directory_path) {
    tsl::robin_map<std::string, GenomeReadData> header_to_read;
    for (auto id : cluster_ids){
        for (const auto &header : cluster_index[id]->read_headers) {
            header_to_read[header] = {};
        }
    }

    fs::remove_all(directory_path);
    fs::create_directories(directory_path);

    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        if (header_to_read.contains(read->header)) {
            header_to_read[read->header] = *read;
        }
    }

    std::map<ClusterID, std::string> mapping;
    for (auto id : cluster_ids){
        std::ofstream cluster_file;
        std::string cluster_file_path = fmt::format("{}/#{}.fq", directory_path.string(), id);

        cluster_file.open(cluster_file_path);
        for (const auto &header : cluster_index[id]->read_headers) {
            cluster_file << header_to_read[header].fastX_string() << std::endl;
        }
        cluster_file.close();

        mapping[id] = cluster_file_path;
    }

    return mapping;
}

void BaseReadClusteringEngine::export_connections(std::vector<ClusterID> &cluster_ids) {
    MIN_SCORE = 10;
    std::vector<ClusterConnection> connections = timeMeasureMemberFunc(&BaseReadClusteringEngine::get_connections, this, "Cluster connections")(cluster_ids);

    std::ofstream conn_file;
    conn_file.open("connections");

    conn_file << cluster_ids.size() << "\n";
    for (ClusterID id : cluster_ids) {
        conn_file << id << " " << *begin(cluster_index[id]->categories) << "\n";
    }

    conn_file << connections.size() << "\n";
    for (auto conn : connections) {
        conn_file << conn.cluster_x_id << " " << conn.cluster_y_id << " " << conn.score << " " << (int) conn.is_good << "\n";
    }

    conn_file.close();
}

void BaseReadClusteringEngine::assemble_clusters(std::vector<ClusterID> &cluster_ids) {
    sort(cluster_ids.rbegin(), cluster_ids.rend(), [this](ClusterID x, ClusterID y) -> bool { return cluster_index[x]->size() < cluster_index[y]->size(); });

    fs::path export_path = "./data/assembly_" + reader->meta.filename;
    auto mapping = export_clusters(cluster_ids, export_path);

    int progress_counter = 1;
    for (ClusterID id : cluster_ids){
        std::string assembled_fasta_path = mapping[id] + "_assembled.fasta";

        std::string cmd = fmt::format("./scripts/assemble_pacbio_cluster.sh {} {}", mapping[id], assembled_fasta_path);
        run_command_with_input(cmd.c_str(), "");

        std::cout << cluster_index[id]->repr() << " assembly: \n";
        try{
            SequenceRecordIterator r(assembled_fasta_path);
            r.show_progress = false;
            std::optional<GenomeReadData> contig;
            while ((contig = r.get_next_record()) != std::nullopt) {
                std::cout << fmt::format("\tContig size {}\n", contig->sequence.length());
            }
        } catch (const std::logic_error &e){
            std::cout << "\t No contigs assembled\n";
        }
    }
}

void BaseReadClusteringEngine::run_clustering() {
    while (clustering_round() > 0){
        print_clusters(-1);

        std::vector<GenomeReadCluster*> clusters_for_display;
        for (auto it = begin(cluster_index); it != end(cluster_index); ++it) {
            if (it->second->size() > 5){
                clusters_for_display.push_back(it->second);
            }
        }
        plot_cluster_coverage(clusters_for_display);
    };

    print_clusters(-1);

    std::vector<ClusterID> cluster_ids;
    std::vector<ClusterID> clusters_for_assembly;
    std::vector<GenomeReadCluster*> clusters_for_display;

    for (auto it = begin(cluster_index); it != end(cluster_index); ++it) {
        if (it->second->size() == 1 && it->second->discriminative_kmer_ids.size() > 5) cluster_ids.push_back(it->first);

        if (it->second->size() > 50) {
            clusters_for_assembly.push_back(it->first);
        }
    }

    //export_connections(cluster_ids);

    assemble_clusters(clusters_for_assembly);

    //std::cout << fmt::format("Exported {} clusters\n", timeMeasureMemberFunc(&BaseReadClusteringEngine::export_clusters, this, "Exporting clusters")(10));
}