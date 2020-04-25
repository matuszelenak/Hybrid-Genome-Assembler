#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <numeric>
#include <experimental/filesystem>
#include <boost/algorithm/string/join.hpp>

#include "KmerIterator.h"
#include "ReadClusteringEngine.h"
#include "Utils.h"

namespace fs = std::experimental::filesystem;
namespace algo = boost::algorithm;

std::mutex connections_mut;
std::mutex merge_mut;
std::mutex component_erase_mut;
std::mutex index_merge;
std::mutex conn_index_mut;
std::mutex index_update;

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


void plot_cluster_coverage(std::vector<GenomeReadCluster *> &clusters) {
    std::map<CategoryID, std::vector<std::string> > category_intervals = {{0, {}},
                                                                          {1, {}}};
    for (auto cluster : clusters) {
        for (int i = 0; i < cluster->endpoints.size(); i += 2) {
            category_intervals[*cluster->categories.begin()].push_back(fmt::format("({}, {})", cluster->endpoints[i].first, cluster->endpoints[i + 1].first));
        }
    }
    std::string input;
    for (const auto &k_v : category_intervals) {
        input += ("[" + algo::join(k_v.second, ", ") + "]\n");
    }
    run_command_with_input("python scripts/plotting.py --plot cluster_coverage", input);
}


std::string ReadClusteringEngine::cluster_consistency(GenomeReadCluster *cluster) {
    std::map<CategoryID, int> category_counts = {{0, 0},
                                                 {1, 0}};
    for (auto read_id : cluster->read_ids) {
        category_counts.insert({read_category_map[read_id], 0}).first->second += 1;
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


void ReadClusteringEngine::print_clusters(int first_n) {
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
        if (cluster->size() > 1) {
            std::cout << fmt::format("#{} {} {} {} disc. kmers\n", cluster->id, cluster_consistency(cluster), cluster->intervals(), cluster->discriminative_kmer_ids.size());
        } else {
            size_one_clusters++;
        }

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << fmt::format("As well as {} clusters of size 1\n", size_one_clusters);
}

void ReadClusteringEngine::construct_read_category_map() {
    this->reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = this->reader->get_next_record()) != std::nullopt) {
        read_category_map.insert({read->id, read->category_id});
    }
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers) {
    this->k = k;
    this->reader = &read_iterator;
    this->kmers = &kmers;

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", cluster_index.size(), reader->meta.records);

    construct_read_category_map();
}

void ReadClusteringEngine::construct_indices_thread() {
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
            ClusterID cluster_id = read->id;

            index_merge.lock();

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (Kmer kmer : in_read_discriminative) {
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer, new_kmer_id));
                if (insert_result.second) {
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_discriminative_ids.push_back(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].push_back(cluster_id);
            }

            std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
            cluster_index.insert(ClusterIndex::value_type(cluster_id, new GenomeReadCluster(*read, in_read_discriminative_ids)));

            index_merge.unlock();
        }
    }
}


int ReadClusteringEngine::construct_indices() {
    reader->rewind();
    auto runner = ThreadRunner(&ReadClusteringEngine::construct_indices_thread, this);

    uint64_t entries = 0;
    for (auto& arr : kmer_cluster_index){
        std::sort(arr.begin(), arr.end());
        entries += arr.size();
    }
    std::cout << fmt::format("Kmer cluster index size {}\n", kmer_cluster_index.size());
    std::cout << fmt::format("{} entries in kmer-cluster index\n", entries);
    return 0;
}


ClusterConnection ReadClusteringEngine::get_connection(GenomeReadCluster *x, GenomeReadCluster *y) {
    int i = 0, j = 0;
    ConnectionScore score = 0;
    while (i < x->discriminative_kmer_ids.size() && j < y->discriminative_kmer_ids.size()) {
        if (x->discriminative_kmer_ids[i] < y->discriminative_kmer_ids[j])
            i++;
        else if (y->discriminative_kmer_ids[j] < x->discriminative_kmer_ids[i])
            j++;
        else
        {
            score++;
            i++;
            j++;
        }
    }
    return {x->id, y->id, score, x->categories == y->categories};
}


void ReadClusteringEngine::get_connections_thread(std::vector<ClusterID> &cluster_indices, std::vector<ClusterConnection> &accumulator, tsl::robin_set<ClusterID> &restricted,
                                                  int &index) {
    while (true) {
        conn_index_mut.lock();
        int pivot_index = index++;
        conn_index_mut.unlock();

        if (pivot_index >= cluster_indices.size()) break;
        ClusterID pivot_id = cluster_indices[pivot_index];

        std::map<ClusterID, ConnectionScore> shared_kmer_counts;

        for (KmerID kmer_id : cluster_index[pivot_id]->discriminative_kmer_ids) {

            auto stop = std::lower_bound(kmer_cluster_index[kmer_id].begin(), kmer_cluster_index[kmer_id].end(), pivot_id);
            for (auto candidate_it = kmer_cluster_index[kmer_id].begin(); candidate_it != stop; ++candidate_it) {
                if (restricted.contains(*candidate_it)) continue;
                shared_kmer_counts.insert(std::map<ClusterID, ConnectionScore>::value_type(*candidate_it, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < MIN_SCORE) continue;
            connections_mut.lock();
            accumulator.push_back({iter->first, pivot_id, iter->second, cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories});
            connections_mut.unlock();
        }
    }
}

std::vector<ClusterConnection> ReadClusteringEngine::get_connections(std::vector<ClusterID> &cluster_ids, tsl::robin_set<ClusterID> &restricted) {
    std::vector<ClusterConnection> connections;

    int index = 0;
    auto r = ThreadRunner(&ReadClusteringEngine::get_connections_thread, this, std::ref(cluster_ids), std::ref(connections), std::ref(restricted), std::ref(index));

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

void ReadClusteringEngine::merge_clusters_thread(std::queue<IDComponent> &component_queue, IndexRemovalMap &for_removal) {
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
        std::vector<std::vector<Endpoint> *> endpoints_to_merge = {&survivor->endpoints};
        for (int i = 1; i < component.size(); i++) {
            survivor->read_ids.insert(survivor->read_ids.end(), cluster_index[component[i]]->read_ids.begin(), cluster_index[component[i]]->read_ids.end());
            survivor->categories.insert(cluster_index[component[i]]->categories.begin(), cluster_index[component[i]]->categories.end());

            arrays_to_merge.push_back(&cluster_index[component[i]]->discriminative_kmer_ids);
            endpoints_to_merge.push_back(&cluster_index[component[i]]->endpoints);
        }

        std::vector<KmerID> merged_discriminative_ids;
        merge_n_vectors(arrays_to_merge, merged_discriminative_ids, true);

        std::vector<Endpoint> merged_endpoints;
        merge_n_vectors(endpoints_to_merge, merged_endpoints, false);

        std::set<Endpoint> final_endpoints;
        int opened_intervals = 0;
        for (auto endpoint : merged_endpoints) {
            if (endpoint.second) {
                opened_intervals++;
                if (opened_intervals == 1) {
                    final_endpoints.insert(endpoint);
                }
            } else {
                opened_intervals--;
                if (opened_intervals == 0) {
                    final_endpoints.insert(endpoint);
                }
            }
        }

        survivor->endpoints.clear();
        for (auto endpoint : final_endpoints) survivor->endpoints.push_back(endpoint);
        survivor->discriminative_kmer_ids = merged_discriminative_ids;

        component_erase_mut.lock();
        for (int i = 1; i < component.size(); i++){
            for (KmerID kmer_id : cluster_index[component[i]]->discriminative_kmer_ids) {
                for_removal.insert({kmer_id, {}}).first.value().push_back(component[i]);
            }
            cluster_index.erase(component[i]);
        }
        component_erase_mut.unlock();
    }
}

void ReadClusteringEngine::kmer_cluster_index_update(IndexRemovalMap::iterator &removal_it, IndexRemovalMap::iterator &removal_end){
    while (true){
        IndexRemovalMap::iterator it;

        index_update.lock();
        if (removal_it != removal_end){
            it = removal_it++;
            index_update.unlock();
        } else {
            index_update.unlock();
            break;
        }

        std::vector<ClusterID>* removal_list = &it.value();
        KmerID kmer_id = it->first;
        std::sort(removal_list->begin(), removal_list->end());

        std::vector<ClusterID> updated_list;
        int i = 0, j = 0;
        while (i < removal_list->size() && j < kmer_cluster_index[kmer_id].size()) {
            if ((*removal_list)[i] < kmer_cluster_index[kmer_id][j]){
                i++;
            }
            else if (kmer_cluster_index[kmer_id][j] < (*removal_list)[i]){
                updated_list.push_back(kmer_cluster_index[kmer_id][j]);
                j++;
            }
            else {  // Should be removed
                i++;
                j++;
            }
        }
        kmer_cluster_index[kmer_id] = updated_list;
    }
}


int ReadClusteringEngine::merge_clusters(const tsl::robin_map<ClusterID, IDComponent> &components) {
    std::queue<IDComponent> component_queue;
    for (auto component_it = begin(components); component_it != end(components); component_it++) {
        if (component_it->second.size() > 1) {
            component_queue.push(component_it->second);
        }
    }
    IndexRemovalMap for_removal;
    auto r = ThreadRunner(&ReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue), std::ref(for_removal));

    auto removal_it = for_removal.begin();
    auto removal_end = for_removal.end();
    auto r2 = ThreadRunner(&ReadClusteringEngine::kmer_cluster_index_update, this, std::ref(removal_it), std::ref(removal_end));
    return 0;
}

ClusterID get_parent(ClusterID cluster_id, tsl::robin_map<ClusterID, ClusterID> &parents) {
    if (cluster_id == parents[cluster_id]) {
        return cluster_id;
    }
    parents[cluster_id] = get_parent(parents[cluster_id], parents);
    return parents[cluster_id];
}

int ReadClusteringEngine::clustering_round() {
    std::vector<ClusterID> cluster_ids;
    std::transform(cluster_index.begin(), cluster_index.end(), std::back_inserter(cluster_ids), [](ClusterIndex::value_type &val) -> ClusterID { return val.first; });

    tsl::robin_set<ClusterID> restricted;
    return clustering_round(cluster_ids, restricted);
}

int ReadClusteringEngine::clustering_round(std::vector<ClusterID> &cluster_ids, tsl::robin_set<ClusterID> &restricted) {
    tsl::robin_map<ClusterID, ClusterID> parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++) {
        parents[cluster_iter->second->id] = cluster_iter->second->id;
        components[cluster_iter->second->id] = {cluster_iter->second->id};
    }

    std::vector<ClusterConnection> cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Cluster connections")(cluster_ids, restricted);
    plot_connection_quality(cluster_connections);

    ConnectionScore min_score;
    std::cout << "Enter minimal score\n";
    std::cin >> min_score;

    int merge_operations = 0;
    for (ClusterConnection &conn : cluster_connections) {
        if (conn.score < min_score) continue;

        ClusterID parent_x = get_parent(conn.cluster_x_id, parents);
        ClusterID parent_y = get_parent(conn.cluster_y_id, parents);
        if (parent_x == parent_y) continue;
        if (restricted.contains(parent_x) && restricted.contains(parent_y)) continue;

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

    timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Merging of clusters")(components);

    return merge_operations;
}

std::map<ClusterID, std::string> ReadClusteringEngine::export_clusters(std::vector<ClusterID> &cluster_ids, fs::path &directory_path) {
    tsl::robin_map<ReadID, GenomeReadData> id_to_read;
    for (auto cluster_id : cluster_ids) {
        for (auto read_id: cluster_index[cluster_id]->read_ids) {
            id_to_read[read_id] = {};
        }
    }

    fs::remove_all(directory_path);
    fs::create_directories(directory_path);

    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        if (id_to_read.contains(read->id)) {
            id_to_read[read->id] = *read;
        }
    }

    std::map<ClusterID, std::string> mapping;
    for (auto id : cluster_ids) {
        std::ofstream cluster_file;
        std::string cluster_file_path = fmt::format("{}/#{}.fq", directory_path.string(), id);

        cluster_file.open(cluster_file_path);
        for (auto read_id : cluster_index[id]->read_ids) {
            cluster_file << id_to_read[read_id].fastX_string() << std::endl;
        }
        cluster_file.close();

        mapping[id] = cluster_file_path;
    }

    return mapping;
}

void ReadClusteringEngine::assemble_clusters(std::vector<ClusterID> &cluster_ids) {
    sort(cluster_ids.rbegin(), cluster_ids.rend(), [this](ClusterID x, ClusterID y) -> bool { return cluster_index[x]->size() < cluster_index[y]->size(); });

    fs::path export_path = "./data/assembly_" + reader->meta.filename;
    auto mapping = export_clusters(cluster_ids, export_path);

    int progress_counter = 1;
    for (ClusterID id : cluster_ids) {
        std::string assembled_fasta_path = mapping[id] + "_assembled.fasta";

        std::string cmd = fmt::format("./scripts/assemble_pacbio_cluster.sh {} {}", mapping[id], assembled_fasta_path);
        run_command_with_input(cmd.c_str(), "");

        std::cout << cluster_index[id]->repr() << " assembly: \n";
        try {
            SequenceRecordIterator r(assembled_fasta_path);
            r.show_progress = false;
            std::optional<GenomeReadData> contig;
            while ((contig = r.get_next_record()) != std::nullopt) {
                std::cout << fmt::format("\tContig size {}\n", contig->sequence.length());
            }
        } catch (const std::logic_error &e) {
            std::cout << "\t No contigs assembled\n";
        }
    }
}

void ReadClusteringEngine::run_clustering() {
    clustering_round();

    std::vector<ClusterID> chain_ids;
    tsl::robin_set<ClusterID> restricted;
    std::vector<GenomeReadCluster *> clusters_for_display;
    for (auto it = begin(cluster_index); it != end(cluster_index); ++it) {
        if (it->second->size() > 5) {
            clusters_for_display.push_back(it->second);
            chain_ids.push_back(it->first);
            restricted.insert(it->first);
        }
    }
    print_clusters(-1);

    plot_cluster_coverage(clusters_for_display);
    while (clustering_round(chain_ids, restricted) > 0) {
        print_clusters(-1);
    }

    plot_cluster_coverage(clusters_for_display);
    std::vector<ClusterConnection> inter_chain_connections;
    for (auto x_id : chain_ids){
        for (auto y_id : chain_ids){
            if (x_id < y_id){
                inter_chain_connections.push_back(get_connection(cluster_index[x_id], cluster_index[y_id]));
            }
        }
    }
    std::sort(inter_chain_connections.rbegin(), inter_chain_connections.rend());
    for (auto conn : inter_chain_connections){
        std::cout << fmt::format("{} {} {} {}\n", conn.cluster_x_id, conn.cluster_y_id, conn.score, (int)conn.is_good);
    }

    std::vector<ClusterID> clusters_for_assembly;
    for (auto it = begin(cluster_index); it != end(cluster_index); ++it) {
        if (it->second->size() > 50) {
            clusters_for_assembly.push_back(it->first);
        }
    }
    assemble_clusters(clusters_for_assembly);
}