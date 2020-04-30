#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <experimental/filesystem>
#include <boost/algorithm/string/join.hpp>

#include "../common/KmerIterator.h"

#include "ReadClusteringEngine.h"

namespace fs = std::experimental::filesystem;
namespace algo = boost::algorithm;

std::mutex connections_mut;
std::mutex component_erase_mut;
std::mutex index_merge;

//DEBUG
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

//DEBUG
void plot_cluster_coverage(std::vector<GenomeReadCluster *> &clusters) {
//    std::map<CategoryID, std::vector<std::string> > category_intervals = {{0, {}},
//                                                                          {1, {}}};
//    for (auto cluster : clusters) {
//        for (int i = 0; i < cluster->endpoints.size(); i += 2) {
//            category_intervals[*cluster->categories.begin()].push_back(fmt::format("({}, {})", cluster->endpoints[i].first, cluster->endpoints[i + 1].first));
//        }
//    }
//    std::string input;
//    for (const auto &k_v : category_intervals) {
//        input += ("[" + algo::join(k_v.second, ", ") + "]\n");
//    }
//    run_command_with_input("python scripts/plotting.py --plot cluster_coverage", input);
}

//DEBUG
std::vector<Interval> ReadClusteringEngine::get_read_coverage_intervals(std::vector<ReadID> &read_ids) {
    std::vector<Interval> result;

    std::vector<Endpoint> endpoints;
    for (auto read_id : read_ids) {
        endpoints.emplace_back(read_intervals[read_id].first, true);
        endpoints.emplace_back(read_intervals[read_id].second, false);
    }
    std::sort(endpoints.begin(), endpoints.end());

    uint32_t current_interval_start;
    int opened_intervals = 0;
    for (auto endpoint : endpoints) {
        if (endpoint.second) {
            opened_intervals++;
            if (opened_intervals == 1) {
                current_interval_start = endpoint.first;
            }
        } else {
            opened_intervals--;
            if (opened_intervals == 0) {
                result.emplace_back(current_interval_start, endpoint.first);
            }
        }
    }
    return result;
}

//DEBUG
void ReadClusteringEngine::print_clusters(std::vector<ClusterID> &ids) {
    sort(ids.rbegin(), ids.rend(), [this](ClusterID x, ClusterID y) -> bool { return cluster_index[x]->size() < cluster_index[y]->size(); });
    for (auto id : ids) {
        std::cout << cluster_index[id]->to_string() << std::endl;
    }
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator) {
    this->reader = &read_iterator;
}

void ReadClusteringEngine::set_kmers(tsl::robin_set<Kmer> &_kmers, int _k){
    this->k = _k;
    this->kmers = new bloom::BloomFilter<Kmer>(_kmers.size(), 0.01);
    for (Kmer kmer : _kmers){
        this->kmers->add(kmer);
    }
}

void ReadClusteringEngine::set_kmers(bloom::BloomFilter<Kmer> *_kmers, int _k){
    this->k = _k;
    this->kmers = _kmers;
}

void ReadClusteringEngine::construct_indices_thread(KmerIndex &kmer_index) {
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        tsl::robin_set<Kmer> in_read_discriminative;
        while (it.next_kmer()) {
            if (kmers->contains(it.current_kmer)) {
                in_read_discriminative.insert(it.current_kmer);
            }
        }

        if (in_read_discriminative.size() > 5) {
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
            ReadMetaData meta = {read->id, read->category_id, read->start, read->end};
            cluster_index.insert(ClusterIndex::value_type(cluster_id, new GenomeReadCluster(meta, in_read_discriminative_ids)));

            read_intervals[read->id] = {read->start, read->end};

            index_merge.unlock();
        }
    }
}

int ReadClusteringEngine::construct_indices() {
    reader->rewind();

    KmerIndex kmer_index;
    run_in_threads(&ReadClusteringEngine::construct_indices_thread, this, std::ref(kmer_index));

    for (auto &arr : kmer_cluster_index) {
        std::sort(arr.begin(), arr.end());
    }
    std::cout << fmt::format("Kmer cluster index size {}\n", kmer_cluster_index.size());
    return 0;
}

void ReadClusteringEngine::get_connections_thread(ConnectionScore min_score, ConcurrentQueue<ClusterID> &cluster_id_queue, std::vector<ClusterConnection> &accumulator) {
    while (true) {
        auto pivot = cluster_id_queue.pop();
        if (pivot == std::nullopt) break;
        ClusterID pivot_id = *pivot;

        std::map<ClusterID, ConnectionScore> shared_kmer_counts;
        for (KmerID kmer_id : cluster_index[pivot_id]->discriminative_kmer_ids) {
            auto stop = std::lower_bound(kmer_cluster_index[kmer_id].begin(), kmer_cluster_index[kmer_id].end(), pivot_id);
            for (auto candidate_it = kmer_cluster_index[kmer_id].begin(); candidate_it != stop; ++candidate_it) {
                shared_kmer_counts.insert(std::map<ClusterID, ConnectionScore>::value_type(*candidate_it, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < min_score) continue;

            connections_mut.lock();
            accumulator.push_back({pivot_id, iter->first, iter->second,
                                   cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories && cluster_index[iter->first]->categories.size() == 1});
            connections_mut.unlock();
        }
    }
}

std::vector<ClusterConnection> ReadClusteringEngine::get_connections(std::vector<ClusterID> &cluster_ids, ConnectionScore min_score) {
    std::vector<ClusterConnection> connections;

    auto queue = ConcurrentQueue<ClusterID>(cluster_ids.begin(), cluster_ids.end());
    run_in_threads(&ReadClusteringEngine::get_connections_thread, this, min_score, std::ref(queue), std::ref(connections));

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

void ReadClusteringEngine::merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal) {
    while (true) {
        auto component_opt = component_queue.pop();
        if (component_opt == std::nullopt) break;
        auto component = *component_opt;

        if (component.size() < 2) continue;
        GenomeReadCluster *survivor = cluster_index[component[0]];

        std::vector<std::vector<KmerID> *> kmer_merge_queue;
        std::set<CategoryID> survivor_categories;
        std::vector<ReadMetaData> survivor_contained_reads;

        for (auto cluster_id : component) {
            auto cluster = cluster_index[cluster_id];
            kmer_merge_queue.push_back(&cluster->discriminative_kmer_ids);
            survivor_contained_reads.insert(survivor_contained_reads.end(), cluster->contained_reads.begin(), cluster->contained_reads.end());
            survivor_categories.insert(cluster->categories.begin(), cluster->categories.end());

            cluster->contained_reads.clear();
        }

        survivor->discriminative_kmer_ids = merge_n_vectors(kmer_merge_queue, true);
        survivor->categories = survivor_categories;
        survivor->contained_reads = survivor_contained_reads;

        component_erase_mut.lock();
        for (auto cluster_id : component) {
            for (KmerID kmer_id : cluster_index[cluster_id]->discriminative_kmer_ids) {
                for_removal.insert({kmer_id, {}}).first.value().push_back(cluster_id);
            }
        }
        component_erase_mut.unlock();
    }
}

void ReadClusteringEngine::kmer_cluster_index_update(ConcurrentQueue<IndexRemovalMap::value_type> &removal_list_queue) {
    while (true) {
        auto removal_pair_opt = removal_list_queue.pop();
        if (removal_pair_opt == std::nullopt) break;

        std::vector<ClusterID> *removal_list = &removal_pair_opt->second;
        KmerID kmer_id = removal_pair_opt->first;
        std::sort(removal_list->begin(), removal_list->end());

        std::vector<ClusterID> updated_list;
        int i = 0, j = 0;
        while (i < removal_list->size() && j < kmer_cluster_index[kmer_id].size()) {
            if ((*removal_list)[i] < kmer_cluster_index[kmer_id][j]) {
                i++;
            } else if (kmer_cluster_index[kmer_id][j] < (*removal_list)[i]) {
                updated_list.push_back(kmer_cluster_index[kmer_id][j]);
                j++;
            } else {  // Should be removed
                i++;
                j++;
            }
        }
        kmer_cluster_index[kmer_id] = updated_list;
    }
}

int ReadClusteringEngine::merge_clusters(std::vector<IDComponent> &components) {
    auto component_queue = ConcurrentQueue<IDComponent>(components.begin(), components.end());

    IndexRemovalMap for_removal;
    run_in_threads(&ReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue), std::ref(for_removal));

    auto removal_list_queue = ConcurrentQueue<IndexRemovalMap::value_type>(for_removal.begin(), for_removal.end());
    run_in_threads(&ReadClusteringEngine::kmer_cluster_index_update, this, std::ref(removal_list_queue));

    return 0;
}

std::vector<std::pair<IDComponent, SpanningTree>>
ReadClusteringEngine::union_find(std::vector<ClusterConnection> &connections, std::set<ClusterID> &restricted, int min_component_size) {
    tsl::robin_set<ClusterID> affected_vertices;
    for (auto &conn : connections) {
        affected_vertices.insert(conn.cluster_x_id);
        affected_vertices.insert(conn.cluster_y_id);
    }

    tsl::robin_map<ClusterID, ClusterID> parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    tsl::robin_map<ClusterID, std::vector<std::pair<ClusterID, ClusterID>>> component_edges;
    tsl::robin_map<ClusterID, bool> component_contains_restricted;
    for (ClusterID vertex : affected_vertices) {
        parents[vertex] = vertex;
        components[vertex] = {vertex};
        component_edges[vertex] = {};
        component_contains_restricted[vertex] = restricted.contains(vertex);
    }

    std::function<ClusterID(ClusterID)> get_parent = [&parents, &get_parent](ClusterID cluster_id) {
        if (cluster_id == parents[cluster_id]) {
            return cluster_id;
        }
        parents[cluster_id] = get_parent(parents[cluster_id]);
        return parents[cluster_id];
    };

    for (ClusterConnection &conn : connections) {
        ClusterID parent_x = get_parent(conn.cluster_x_id);
        ClusterID parent_y = get_parent(conn.cluster_y_id);
        if (parent_x == parent_y) continue;
        if (component_contains_restricted[parent_x] && component_contains_restricted[parent_y]) continue;

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

        component_edges[bigger].emplace_back(conn.cluster_x_id, conn.cluster_y_id);
        component_edges[bigger].insert(component_edges[bigger].end(), component_edges[smaller].begin(), component_edges[smaller].end());
        component_edges.erase(smaller);

        component_contains_restricted[bigger] |= component_contains_restricted[smaller];
        component_contains_restricted.erase(smaller);
    }

    std::vector<std::pair<IDComponent, SpanningTree>> result;
    for (const auto &comp : components) {
        if (comp.second.size() >= min_component_size) {
            result.emplace_back(comp.second, component_edges[comp.first]);
        }
    }
    return result;
}

std::map<ClusterID, std::string> ReadClusteringEngine::export_clusters(std::vector<ClusterID> &cluster_ids, fs::path &directory_path) {
    tsl::robin_map<ReadID, GenomeReadData> id_to_read;

    fs::remove_all(directory_path);
    fs::create_directories(directory_path);
    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        id_to_read[read->id] = *read;
    }

    std::map<ClusterID, std::string> mapping;
    for (auto id : cluster_ids) {
        std::ofstream cluster_file;
        std::string cluster_file_path = fmt::format("{}/#{}.fa", directory_path.string(), id);

        cluster_file.open(cluster_file_path);
        for (auto read_meta : cluster_index[id]->contained_reads) {
            cluster_file << id_to_read[read_meta.id].fasta_string() << std::endl;
        }

        cluster_file.close();
        mapping[id] = cluster_file_path;
    }

    return mapping;
}

std::map<ClusterID, std::vector<std::string>> ReadClusteringEngine::assemble_clusters(std::vector<ClusterID> &cluster_ids) {
    sort(cluster_ids.rbegin(), cluster_ids.rend(), [this](ClusterID x, ClusterID y) -> bool { return cluster_index[x]->size() < cluster_index[y]->size(); });

    fs::path export_path = "./data/assembly_" + reader->meta.filename;
    auto mapping = export_clusters(cluster_ids, export_path);

    std::map<ClusterID, std::vector<std::string>> assembly_mapping;
    for (ClusterID id : cluster_ids) {
        std::string assembled_fasta_path = mapping[id] + "_assembled.fa";

        std::string cmd = fmt::format("./scripts/assemble_cluster.sh {} {}", mapping[id], assembled_fasta_path);
        run_command_with_input(cmd.c_str(), "");

        std::cout << cluster_index[id]->to_string() << " : \n";
        try {
            SequenceRecordIterator r(assembled_fasta_path);

            assembly_mapping[id] = {};
            r.show_progress = false;
            std::optional<GenomeReadData> contig;
            while ((contig = r.get_next_record()) != std::nullopt) {
                std::cout << fmt::format("\tContig size {}\n", contig->sequence.length());
                assembly_mapping[id].push_back(contig->sequence);
            }
        } catch (const std::logic_error &e) {
            std::cout << "\t No contigs assembled\n";
        }
    }
    return assembly_mapping;
}

tsl::robin_map<ReadID, int> distance_bfs(ReadID starting_vertex, NeighborMap &neighbor_map) {
    std::queue<ReadID> bfs_queue;
    bfs_queue.push(starting_vertex);
    tsl::robin_set<ReadID> visited = {starting_vertex};
    tsl::robin_map<ReadID, int> distances;
    int dist = 0;
    while (!bfs_queue.empty()) {
        auto vertex = bfs_queue.front();
        distances[vertex] = dist;
        bfs_queue.pop();
        if (neighbor_map.contains(vertex)) {
            for (ReadID neighbor : neighbor_map[vertex]) {
                if (!visited.contains(neighbor)) {
                    visited.insert(neighbor);
                    bfs_queue.push(neighbor);
                }
            }
        }
        dist++;
    }
    return distances;
}

std::pair<std::vector<ReadID>, std::vector<ReadID>> ReadClusteringEngine::get_spanning_tree_boundary_reads(SpanningTree &tree) {
    // convert the edges into adjacency map
    NeighborMap neighbor_map;
    for (auto edge : tree) {
        neighbor_map.insert({edge.first, {}}).first.value().insert(edge.second);
        neighbor_map.insert({edge.second, {}}).first.value().insert(edge.first);
    }

    // Run BFS to figure out the vertex that has the greatest distance from other vertex in the graph
    auto distances = distance_bfs(neighbor_map.begin()->first, neighbor_map);
    auto left_side_vertex_and_dist = *std::max_element(distances.begin(), distances.end(),
                                                       [](const std::pair<ReadID, int> &p1, const std::pair<ReadID, int> &p2) { return p1.second < p2.second; });

    // Run BFS again to figure out the opposing vertex to the first one
    distances = distance_bfs(left_side_vertex_and_dist.first, neighbor_map);
    auto right_side_vertex_and_dist = *std::max_element(distances.begin(), distances.end(),
                                                        [](const std::pair<ReadID, int> &p1, const std::pair<ReadID, int> &p2) { return p1.second < p2.second; });

    int tail_length = 15;
    // Gather vertices that are about as distant as the two vertices found before and merge their kmers into one vector
    std::vector<ReadID> right_side_vertices;
    for (auto vertex_dist_pair : distances) {
        if (vertex_dist_pair.second + tail_length > right_side_vertex_and_dist.second) {
            right_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    distances = distance_bfs(right_side_vertex_and_dist.first, neighbor_map);
    left_side_vertex_and_dist = *std::max_element(distances.begin(), distances.end(),
                                                  [](const std::pair<ClusterID, int> &p1, const std::pair<ReadID, int> &p2) { return p1.second < p2.second; });

    std::vector<ReadID> left_side_vertices;
    for (auto vertex_dist_pair : distances) {
        if (vertex_dist_pair.second + tail_length > left_side_vertex_and_dist.second) {
            left_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

//    for (auto interval : get_read_coverage_intervals(left_side_vertices)){
//        std::cout << fmt::format("{}-{}  ", interval.first, interval.second);
//    }
//    std::cout << " / ";
//    for (auto interval : get_read_coverage_intervals(right_side_vertices)){
//        std::cout << fmt::format("{}-{}  ", interval.first, interval.second);
//    }
//    std::cout << std::endl;

    return {left_side_vertices, right_side_vertices};
}

std::pair<std::vector<KmerID>, std::vector<KmerID>> ReadClusteringEngine::get_core_cluster_kmer_tails(SpanningTree &tree) {
    auto core_cluster_tails = get_spanning_tree_boundary_reads(tree);

    auto left_tail_connections = get_connections(core_cluster_tails.first, 20);
    //plot_connection_quality(left_tail_connections);
    std::set<ReadID> left_tail_vertex_set(core_cluster_tails.first.begin(), core_cluster_tails.first.end());
    for (auto conn : left_tail_connections) {
        left_tail_vertex_set.insert(conn.cluster_x_id);
        left_tail_vertex_set.insert(conn.cluster_y_id);
    }

    auto right_tail_connections = get_connections(core_cluster_tails.second, 20);
    //plot_connection_quality(right_tail_connections);
    std::set<ReadID> right_tail_vertex_set(core_cluster_tails.second.begin(), core_cluster_tails.second.end());
    for (auto conn : right_tail_connections) {
        right_tail_vertex_set.insert(conn.cluster_x_id);
        right_tail_vertex_set.insert(conn.cluster_y_id);
    }

    std::vector<std::vector<KmerID> *> left_side_kmer_merge_queue;
    for (auto cluster_id : left_tail_vertex_set) {
        left_side_kmer_merge_queue.push_back(&cluster_index[cluster_id]->discriminative_kmer_ids);
    }
    std::vector<std::vector<KmerID> *> right_side_kmer_merge_queue;
    for (auto cluster_id : right_tail_vertex_set) {
        right_side_kmer_merge_queue.push_back(&cluster_index[cluster_id]->discriminative_kmer_ids);
    }
    return {merge_n_vectors(left_side_kmer_merge_queue, true), merge_n_vectors(right_side_kmer_merge_queue, true)};
}

std::vector<ClusterConnection> ReadClusteringEngine::get_core_cluster_connections(std::vector<std::pair<IDComponent, SpanningTree>> &components_with_trees) {
    std::multimap<ClusterID, std::vector<KmerID>> core_cluster_twin_tails;
    for (auto &component_and_tree : components_with_trees) {
        auto twin_tail_kmers = get_core_cluster_kmer_tails(component_and_tree.second);
        core_cluster_twin_tails.insert({component_and_tree.first[0], twin_tail_kmers.first});
        core_cluster_twin_tails.insert({component_and_tree.first[0], twin_tail_kmers.second});
    }

    std::vector<ClusterConnection> tail_connecting_edges;
    for (auto c_id_kmers_pair : core_cluster_twin_tails) {
        for (auto c_id_kmers_pair_2 : core_cluster_twin_tails) {
            if (c_id_kmers_pair.first < c_id_kmers_pair_2.first) {
                auto kmer_intersection = get_vectors_intersection(c_id_kmers_pair.second, c_id_kmers_pair_2.second);

                bool is_good = cluster_index[c_id_kmers_pair.first]->categories == cluster_index[c_id_kmers_pair_2.first]->categories;
                tail_connecting_edges.push_back({c_id_kmers_pair.first, c_id_kmers_pair_2.first, kmer_intersection.size(), is_good});
            }
        }
    }
    std::sort(tail_connecting_edges.rbegin(), tail_connecting_edges.rend());
    for (auto edge : tail_connecting_edges) {
        std::cout << fmt::format("{} : {} between {} and {}\n", edge.is_good ? '+' : '-', edge.score, cluster_index[edge.cluster_x_id]->to_string(),
                                 cluster_index[edge.cluster_y_id]->to_string());
        if (!edge.is_good) break;
    }

    ConnectionScore min_score;
    std::cin >> min_score;
    return filter_connections(tail_connecting_edges, [min_score](ClusterConnection &conn) { return conn.score >= min_score; });
}

std::vector<ClusterID> ReadClusteringEngine::run_clustering() {
    auto extract_components = [](std::vector<std::pair<IDComponent, SpanningTree>> &union_find_output) {
        std::vector<IDComponent> result;
        std::transform(union_find_output.begin(), union_find_output.end(), std::back_inserter(result),
                       [](std::pair<IDComponent, SpanningTree> &p) { return p.first; }
        );
        return result;
    };

    auto filter_clusters = [this](const std::function<bool(GenomeReadCluster *)> &func) {
        std::vector<ClusterID> result;
        for (auto cluster_pair : cluster_index) {
            if (func(cluster_pair.second)) {
                result.push_back(cluster_pair.first);
            }
        }
        return result;
    };

    if (kmers == nullptr || k == 0) throw std::logic_error("You need to set the discriminative kmers first");

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", cluster_index.size(), reader->meta.records);

    int magic = 70;
    auto strong_initial_clusters = filter_clusters([magic](GenomeReadCluster *c) { return c->discriminative_kmer_ids.size() > magic; });
    auto cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Strong connections")(strong_initial_clusters, magic);

    std::set<ClusterID> restricted;
    auto components_and_trees = union_find(cluster_connections, restricted, 30);
    std::vector<IDComponent> components = extract_components(components_and_trees);
    timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Creation of strong clusters")(components);

    std::vector<ClusterID> core_cluster_ids;
    for (const IDComponent &component : components) {
        core_cluster_ids.push_back(component[0]);
    }
    print_clusters(core_cluster_ids);


    auto core_cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_core_cluster_connections, this, "Core connections")(components_and_trees);
    components_and_trees = union_find(core_cluster_connections, restricted, 1);
    components = extract_components(components_and_trees);
    timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Merging of core clusters")(components);

    remove_merged_clusters();
    core_cluster_ids = filter_clusters([](GenomeReadCluster *c) { return c->size() > 10; });

    print_clusters(core_cluster_ids);

    for (int i = 0; i < 3; i++) {
        cluster_connections = get_connections(core_cluster_ids, 1);

        plot_connection_quality(cluster_connections);

        ConnectionScore min_score;
        std::cin >> min_score;
        if (min_score == 0) break;
        cluster_connections = filter_connections(cluster_connections, [min_score](ClusterConnection &conn) { return conn.score >= min_score; });

        restricted.clear();
        restricted.insert(core_cluster_ids.begin(), core_cluster_ids.end());
        components_and_trees = union_find(cluster_connections, restricted, 1);
        components = extract_components(components_and_trees);
        merge_clusters(components);
        remove_merged_clusters();

        core_cluster_ids = filter_clusters([](GenomeReadCluster *c) { return c->size() > 30; });
    }

    print_clusters(core_cluster_ids);

    return core_cluster_ids;
}

ReadClusteringEngine::~ReadClusteringEngine() {
    for (auto cluster_pair : cluster_index) {
        free(cluster_pair.second);
    }
}
