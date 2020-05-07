#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <experimental/filesystem>
#include <boost/algorithm/string/join.hpp>
#include <stack>

#include "../common/KmerIterator.h"

#include "ReadClusteringEngine.h"

namespace fs = std::experimental::filesystem;
namespace algo = boost::algorithm;

std::mutex connections_mut;
std::mutex component_merging_mut;
std::mutex index_merge;

//DEBUG
void plot_connection_quality(std::vector<ClusterConnection> &connections) {
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching = {{true, {}}, {false, {}}};
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

void plot_connection_quality(std::vector<ClusterConnection> &connections, std::vector<ClusterID> &core_ids){
    std::map<ClusterID, std::pair<ConnectionScore, bool> > best_per_vertex;
    for (auto & conn : connections){
        auto insert_pair = best_per_vertex.insert({conn.cluster_y_id, {0, false}}).first;
        if (conn.score > insert_pair->second.first){
            best_per_vertex[conn.cluster_y_id] = {conn.score, conn.is_good};
        }
    }
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching = {{true, {}}, {false, {}}};
    for (auto _score_and_quality : best_per_vertex){
        score_to_matching[_score_and_quality.second.second].insert(std::map<ConnectionScore, int>::value_type(_score_and_quality.second.first, 0)).first->second += 1;
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

void ReadClusteringEngine::plot_read_overlaps(std::vector<ClusterConnection> &connections){
    std::map<ConnectionScore, std::map<int, int> > overlaps_for_strength;
    for (auto & conn : connections){
        auto a = read_intervals[conn.cluster_x_id];
        auto b = read_intervals[conn.cluster_y_id];

        int overlap = std::max(0u, std::min(a.second, b.second) - std::max(a.first, b.first));
        if (overlap < 1) continue;

        if (!overlaps_for_strength.contains(conn.score)) overlaps_for_strength[conn.score] = {};
        if (!overlaps_for_strength[conn.score].contains(overlap)) overlaps_for_strength[conn.score][overlap] = 0;
        overlaps_for_strength[conn.score][overlap] += 1;
    }

    std::vector<std::string> per_score_overlaps;
    for (auto score_overlaps : overlaps_for_strength){
        std::vector<std::string> overlap_counts;
        for (auto overlap_count : score_overlaps.second){
            overlap_counts.push_back(fmt::format("{}: {}", overlap_count.first, overlap_count.second));
        }
        per_score_overlaps.push_back(fmt::format("({}, {{{}}})", score_overlaps.first, boost::algorithm::join(overlap_counts, ", ")));
    }
    std::string plot_input = fmt::format("{}\n{}", per_score_overlaps.size(), boost::algorithm::join(per_score_overlaps, "\n"));
    run_command_with_input("python scripts/plotting.py --plot connection_overlaps", plot_input);
}

void ReadClusteringEngine::export_spanning_tree(SpanningTree &tree){
    std::set<ClusterID> nodes;
    std::vector<std::string> node_strings;
    std::vector<std::string> edge_strings;
    for (auto edge : tree){
        nodes.insert(edge.first);
        nodes.insert(edge.second);
        edge_strings.push_back(fmt::format("{{source: {}, target: {}}}", edge.first, edge.second));
    }
    for (ClusterID node : nodes){
        node_strings.push_back(fmt::format("{{ id: {}, start : {}, end : {} }}", node, read_intervals[node].first, read_intervals[node].second));
    }

    std::ofstream output;
    output.open(fmt::format("./data/trees/{}_tree_{}.js", reader->meta.filename, *nodes.begin()));
    output << fmt::format("nodes = [{}]\nedges = [{}]\n", boost::algorithm::join(node_strings, ", "), boost::algorithm::join(edge_strings, ", "));
    output.close();
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
    std::cout << "Clusters:\n";
    for (auto id : ids) {
        std::cout << cluster_index[id]->to_string() << std::endl;
    }
    std::cout << "End clusters\n";
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

        std::vector<std::pair<Kmer, int>> in_read_discriminative;
        while (it.next_kmer()) {
            if (kmers->contains(it.current_kmer)) {
                in_read_discriminative.push_back({it.current_kmer, it.position_in_sequence});
            }
        }

        if (in_read_discriminative.size() > 1) {
            std::vector<KmerID> in_read_discriminative_ids;
            ClusterID cluster_id = read->id;

            index_merge.lock();
            kmer_positions.insert({read->id, {}});

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (auto kmer_position_pair : in_read_discriminative) {
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer_position_pair.first, new_kmer_id));
                if (insert_result.second) {
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_discriminative_ids.push_back(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].push_back(cluster_id);
                kmer_positions[read->id].insert({insert_result.first->second, kmer_position_pair.second});
            }

            std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
            ReadMetaData meta = {read->id, read->category_id, read->start, read->end};
            cluster_index.insert(ClusterIndex::value_type(cluster_id, new GenomeReadCluster(meta, in_read_discriminative_ids)));
            read_intervals[read->id] = {read->start, read->end};

            read_lengths.insert({read->id, read->sequence.length()});
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
            for (auto candidate : kmer_cluster_index[kmer_id]) {
                shared_kmer_counts.insert(std::map<ClusterID, ConnectionScore>::value_type(candidate, 0)).first->second += 1;
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

void ReadClusteringEngine::merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal, std::vector<ClusterID> &merged_ids) {
    while (true) {
        auto component_opt = component_queue.pop();
        if (component_opt == std::nullopt) break;
        auto component = *component_opt;

        if (component.size() == 1){
            component_merging_mut.lock();
            merged_ids.push_back(component[0]);
            component_merging_mut.unlock();
            continue;
        }
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

        component_merging_mut.lock();
        merged_ids.push_back(component[0]);
        for (auto cluster_id : component) {
            for (KmerID kmer_id : cluster_index[cluster_id]->discriminative_kmer_ids) {
                for_removal.insert({kmer_id, {}}).first.value().push_back(cluster_id);
            }
        }
        component_merging_mut.unlock();
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

std::vector<ClusterID> ReadClusteringEngine::merge_clusters(std::vector<IDComponent> &components) {
    auto component_queue = ConcurrentQueue<IDComponent>(components.begin(), components.end());

    IndexRemovalMap for_removal;
    std::vector<ClusterID> merged_ids;
    run_in_threads(&ReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue), std::ref(for_removal), std::ref(merged_ids));

    auto removal_list_queue = ConcurrentQueue<IndexRemovalMap::value_type>(for_removal.begin(), for_removal.end());
    run_in_threads(&ReadClusteringEngine::kmer_cluster_index_update, this, std::ref(removal_list_queue));

    return merged_ids;
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

int ReadClusteringEngine::approximate_read_overlap(ReadID x, ReadID y){
    auto shared_kmers = get_vectors_intersection(cluster_index[x]->discriminative_kmer_ids, cluster_index[y]->discriminative_kmer_ids);

    std::vector<int> x_positions, y_positions;
    for (auto kmer_id : shared_kmers){
        x_positions.push_back(kmer_positions[x][kmer_id]);
        y_positions.push_back(kmer_positions[y][kmer_id]);
    }
    int max_x, max_y, min_x, min_y;
    max_x = *std::max_element(x_positions.begin(), x_positions.end());
    max_y = *std::max_element(y_positions.begin(), y_positions.end());
    min_x = *std::min_element(x_positions.begin(), x_positions.end());
    min_y = *std::min_element(y_positions.begin(), y_positions.end());

    auto overlap = std::max(max_x - min_x, max_y - min_y);

//    auto a = read_intervals[x];
//    auto b = read_intervals[y];
//    int real_overlap = std::max(0u, std::min(a.second, b.second) - std::max(a.first, b.first));
//    if (real_overlap < 0){
//        std::cout << fmt::format("Estimated {} real {} for #{}({})[{}-{}] and #{}({})[{}-{}]\n",
//                overlap, real_overlap,
//                x, read_lengths[x], read_intervals[x].first, read_intervals[x].second,
//                y, read_lengths[y], read_intervals[y].first, read_intervals[y].second);
//    }
    return overlap;
}

std::pair<std::vector<ReadID>, std::vector<ReadID>> ReadClusteringEngine::get_spanning_tree_boundary_reads(SpanningTree &tree) {
    // convert the edges into adjacency map
    tsl::robin_set<ReadID> vertices;
    tsl::robin_map<ReadID, tsl::robin_map<ReadID, int>> adjacency_map;
    for (auto edge : tree) {
        auto dist = approximate_read_overlap(edge.first, edge.second);
        adjacency_map.insert({edge.first, {}}).first.value().insert({edge.second, dist});
        adjacency_map.insert({edge.second, {}}).first.value().insert({edge.first, dist});
        vertices.insert(edge.first);
        vertices.insert(edge.second);
    }

    auto distance_bfs = [this, &adjacency_map](ReadID starting_vertex){
        std::stack<ReadID> stack;
        stack.push(starting_vertex);

        std::queue<ReadID> bfs_queue;
        bfs_queue.push(starting_vertex);
        tsl::robin_set<ReadID> visited;
        tsl::robin_map<ReadID, int> distances = {{starting_vertex, read_lengths[starting_vertex]}};
        while (!bfs_queue.empty()) {
            auto vertex = bfs_queue.front();
            visited.insert(vertex);
            bfs_queue.pop();
            if (adjacency_map.contains(vertex)) {
                for (auto adjacency : adjacency_map[vertex]) {
                    if (!visited.contains(adjacency.first)) {
                        distances[adjacency.first] = distances[vertex] + read_lengths[adjacency.first] - adjacency.second;
                        bfs_queue.push(adjacency.first);
                    }
                }
            }
        }
        return distances;
    };

    auto get_max_distance_pair = [](tsl::robin_map<ReadID, int> &distances){
        return *std::max_element(distances.begin(), distances.end(),[](const std::pair<ReadID, int> &p1, const std::pair<ReadID, int> &p2) {
            return p1.second < p2.second;
        });
    };

    ReadID edge_most_vertex = *std::max_element(vertices.begin(), vertices.end(), [this](const ReadID &x, const ReadID &y){ return read_intervals[x].second < read_intervals[y].second; });
    auto true_distances = distance_bfs(edge_most_vertex);
    auto edge_most_vertex_pair = get_max_distance_pair(true_distances);

    std::cout << fmt::format("Distances {} for #{}({})[{}-{}] and #{}({})[{}-{}]\n",
                edge_most_vertex_pair.second, edge_most_vertex, read_lengths[edge_most_vertex], read_intervals[edge_most_vertex].first, read_intervals[edge_most_vertex].second,
                edge_most_vertex_pair.first, read_lengths[edge_most_vertex_pair.first], read_intervals[edge_most_vertex_pair.first].first, read_intervals[edge_most_vertex_pair.first].second);


    // Run BFS to figure out the vertex that has the greatest distance from other vertex in the graph
    auto distances = distance_bfs(adjacency_map.begin()->first);
    auto left_side_vertex_and_dist = get_max_distance_pair(distances);

    // Run BFS again to figure out the opposing vertex to the first one

    // Gather vertices that are about as distant as the two vertices found before and merge their kmers into one vector
    std::vector<ReadID> left_side_vertices, right_side_vertices;

    distances = distance_bfs(left_side_vertex_and_dist.first);
    auto right_side_vertex_and_dist = get_max_distance_pair(distances);
    int tail_length = std::max(right_side_vertex_and_dist.second * 0.03, 1000.0);
    for (auto vertex_dist_pair : distances) {
        if (vertex_dist_pair.second + tail_length > right_side_vertex_and_dist.second) {
            right_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    distances = distance_bfs(right_side_vertex_and_dist.first);
    left_side_vertex_and_dist = get_max_distance_pair(distances);
    tail_length = std::max(left_side_vertex_and_dist.second * 0.03, 1000.0);
    for (auto vertex_dist_pair : distances) {
        if (vertex_dist_pair.second + tail_length > left_side_vertex_and_dist.second) {
            left_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    std::cout << fmt::format("Computed distance {} for [{}-{}] and [{}-{}]\n",
                             right_side_vertex_and_dist.second,
                             read_intervals[right_side_vertex_and_dist.first].first, read_intervals[right_side_vertex_and_dist.first].second,
                             read_intervals[left_side_vertex_and_dist.first].first, read_intervals[left_side_vertex_and_dist.first].second);


    for (auto interval : get_read_coverage_intervals(right_side_vertices)){
        std::cout << fmt::format("{}-{} ", interval.first, interval.second);
    } std::cout << " / ";
    for (auto interval : get_read_coverage_intervals(left_side_vertices)){
        std::cout << fmt::format("{}-{} ", interval.first, interval.second);
    } std::cout << "\n\n";

    return {left_side_vertices, right_side_vertices};
}

std::pair<std::vector<KmerID>, std::vector<KmerID>> ReadClusteringEngine::get_core_cluster_kmer_tails(SpanningTree &tree) {
    auto core_cluster_tails = get_spanning_tree_boundary_reads(tree);

    auto left_tail_connections = get_connections(core_cluster_tails.first, 40);
    //plot_connection_quality(left_tail_connections);
    std::set<ReadID> left_tail_vertex_set(core_cluster_tails.first.begin(), core_cluster_tails.first.end());
    for (auto conn : left_tail_connections) {
        left_tail_vertex_set.insert(conn.cluster_x_id);
        left_tail_vertex_set.insert(conn.cluster_y_id);
    }

    auto right_tail_connections = get_connections(core_cluster_tails.second, 40);
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
    for (auto &component_and_tree : components_with_trees){
        export_spanning_tree(component_and_tree.second);
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
        //if (!edge.is_good) break;
        if (edge.score < 100) break;
    }

    ConnectionScore min_score;
    std::cin >> min_score;
    return filter_connections(tail_connecting_edges, [min_score](ClusterConnection &conn) { return conn.score >= min_score; });
}

std::vector<ClusterID> ReadClusteringEngine::run_clustering() {
    if (kmers == nullptr || k == 0) throw std::logic_error("You need to set the discriminative kmers first");

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

    auto get_core_component_ids = [this](int threshold_size){
        std::vector<ClusterID> result;
        for (auto id_ptr_pair : cluster_index){
            if (id_ptr_pair.second->size() >= threshold_size){
                result.push_back(id_ptr_pair.first);
            }
        }
        return result;
    };

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", cluster_index.size(), reader->meta.records);

    int magic = 5;
    auto strong_initial_clusters = filter_clusters([magic](GenomeReadCluster *c) { return c->discriminative_kmer_ids.size() > magic; });
    auto cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Strong connections")(strong_initial_clusters, magic);
    plot_connection_quality(cluster_connections);
    //plot_read_overlaps(cluster_connections);

    ConnectionScore min;
    std::cin >> min;
    auto strong_connections = filter_connections(cluster_connections, [min](ClusterConnection &conn) { return conn.score >= min; });

    std::set<ClusterID> restricted;
    auto components_and_trees = union_find(strong_connections, restricted, 5);
    std::vector<IDComponent> components = extract_components(components_and_trees);
    auto core_cluster_ids = timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Creation of strong clusters")(components);
    print_clusters(core_cluster_ids);


    auto core_cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_core_cluster_connections, this, "Core connections")(components_and_trees);
    if (!core_cluster_connections.empty()){
        components_and_trees = union_find(core_cluster_connections, restricted, 1);
        components = extract_components(components_and_trees);
        timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Merging of core clusters")(components);
        remove_merged_clusters();
        core_cluster_ids = get_core_component_ids(30);
    } else {
        remove_merged_clusters();
    }

    print_clusters(core_cluster_ids);

    for (int i = 0; i < 3; i++) {
        cluster_connections = get_connections(core_cluster_ids, 1);
        plot_connection_quality(cluster_connections, core_cluster_ids);

        ConnectionScore min_score;
        std::cin >> min_score;
        cluster_connections = filter_connections(cluster_connections, [min_score](ClusterConnection &conn) { return conn.score >= min_score; });
        if (cluster_connections.empty()) break;

        restricted.clear();
        restricted.insert(core_cluster_ids.begin(), core_cluster_ids.end());
        components_and_trees = union_find(cluster_connections, restricted, 1);
        components = extract_components(components_and_trees);
        merge_clusters(components);
        remove_merged_clusters();

        core_cluster_ids = get_core_component_ids(30);

        print_clusters(core_cluster_ids);
    }

    return core_cluster_ids;
}

ReadClusteringEngine::~ReadClusteringEngine() {
    for (auto cluster_pair : cluster_index) {
        free(cluster_pair.second);
    }
}
