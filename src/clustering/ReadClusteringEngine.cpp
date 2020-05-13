#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <filesystem>
#include <boost/algorithm/string/join.hpp>
#include <stack>
#include <unordered_set>
#include <eigen2/Eigen/Core>

#include "../common/KmerIterator.h"

#include "ReadClusteringEngine.h"
#include "../lib/clustering/SpectralClustering.h"

namespace algo = boost::algorithm;


//DEBUG
void plot_connection_quality(std::vector<ComponentConnection> &connections) {
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching = {{true,  {}},
                                                                        {false, {}}};
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
    //std::cout << hist_input << std::endl;
    run_command_with_input("python scripts/plotting.py --plot connection_histogram", hist_input);
}


void plot_connection_quality(std::vector<ComponentConnection> &connections, std::vector<ComponentID> &core_ids) {
    std::map<ComponentID, std::pair<ConnectionScore, bool> > best_per_vertex;
    for (auto &conn : connections) {
        auto insert_pair = best_per_vertex.insert({conn.component_y_id, {0, false}}).first;
        if (conn.score > insert_pair->second.first) {
            best_per_vertex[conn.component_y_id] = {conn.score, conn.is_good};
        }
    }
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching = {{true,  {}},
                                                                        {false, {}}};
    for (auto _score_and_quality : best_per_vertex) {
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

std::vector<ComponentConnection> filter_connections(std::vector<ComponentConnection> &original, const std::function<bool(ComponentConnection &)> &func) {
    std::vector<ComponentConnection> result;
    std::copy_if(original.begin(), original.end(), std::back_inserter(result), func);
    return result;
}

//DEBUG
void ReadClusteringEngine::plot_read_overlaps(std::vector<ComponentConnection> &connections) {
    std::map<ConnectionScore, std::map<int, int> > overlaps_for_strength;
    for (auto &conn : connections) {
        auto a = read_intervals[conn.component_x_id];
        auto b = read_intervals[conn.component_y_id];

        int overlap = std::max(0u, std::min(a.second, b.second) - std::max(a.first, b.first));
        if (overlap < 1) continue;

        if (!overlaps_for_strength.contains(conn.score)) overlaps_for_strength[conn.score] = {};
        if (!overlaps_for_strength[conn.score].contains(overlap)) overlaps_for_strength[conn.score][overlap] = 0;
        overlaps_for_strength[conn.score][overlap] += 1;
    }

    std::vector<std::string> per_score_overlaps;
    for (auto score_overlaps : overlaps_for_strength) {
        std::vector<std::string> overlap_counts;
        for (auto overlap_count : score_overlaps.second) {
            overlap_counts.push_back(fmt::format("{}: {}", overlap_count.first, overlap_count.second));
        }
        per_score_overlaps.push_back(fmt::format("({}, {{{}}})", score_overlaps.first, boost::algorithm::join(overlap_counts, ", ")));
    }
    std::string plot_input = fmt::format("{}\n{}", per_score_overlaps.size(), boost::algorithm::join(per_score_overlaps, "\n"));
    run_command_with_input("python scripts/plotting.py --plot connection_overlaps", plot_input);
}

//DEBUG
void ReadClusteringEngine::plot_overlap_differences(std::vector<ComponentConnection> &connections) {
    auto get_real_overlap = [this](ReadID x, ReadID y) {
        auto a = read_intervals[x];
        auto b = read_intervals[y];
        int real_overlap = std::max(0u, std::min(a.second, b.second) - std::max(a.first, b.first));
        return real_overlap;
    };

    std::ofstream dump;
    dump.open("overlaps");
    for (auto &conn : connections) {
        auto overlap = approximate_read_overlap(conn.component_x_id, conn.component_y_id);
        int real_overlap = get_real_overlap(conn.component_x_id, conn.component_y_id);

        if (real_overlap > 1) {
            dump << fmt::format("{} {}\n", overlap, real_overlap);
        }
    }
    dump.close();
}

std::vector<std::pair<ComponentID, ConnectionScore>> score_from_kmer_occurrences(std::vector<std::vector<ComponentID>*> &occurrence_arrays, ConnectionScore min_value){
    auto q_compare = [](std::pair<ComponentID, ConnectionScore> &x, std::pair<ComponentID, ConnectionScore> &y){
        return x.first > y.first;
    };
    std::priority_queue<
            std::pair<ComponentID, ConnectionScore>,
            std::vector<std::pair<ComponentID, ConnectionScore>>,
            std::function<bool(std::pair<ComponentID, ConnectionScore>&, std::pair<ComponentID, ConnectionScore>&)>> q(q_compare);

    for (int i = 0; i < occurrence_arrays.size(); i++) {
        q.push({(*occurrence_arrays[i])[0], i});
    }
    std::vector<int> next_index(occurrence_arrays.size(), 1);

    std::vector<std::pair<ComponentID, ConnectionScore>> result;

    ConnectionScore current_component_count = 0;
    ComponentID current_component = q.top().first;

    while (!q.empty()){
        auto top_pair = q.top();
        q.pop();
        if (top_pair.first != current_component){
            if (current_component_count >= min_value){
                result.emplace_back(current_component, current_component_count);
            }
            current_component_count = 1;
            current_component = top_pair.first;
        } else {
            current_component_count++;
        }
        if (next_index[top_pair.second] < occurrence_arrays[top_pair.second]->size()){
            q.push({(*occurrence_arrays[top_pair.second])[next_index[top_pair.second]], top_pair.second});
            next_index[top_pair.second]++;
        }
    }
    if (current_component_count >= min_value){
        result.emplace_back(current_component, current_component_count);
    }
    return result;
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
void ReadClusteringEngine::print_components(std::vector<ComponentID> &ids) {
    if (!debug) return;

    sort(ids.rbegin(), ids.rend(), [this](ComponentID x, ComponentID y) -> bool { return component_index[x]->size() < component_index[y]->size(); });
    std::cout << fmt::format("### Printing {} components ###\n", ids.size());
    for (auto id : ids) {
        std::cout << component_index[id]->to_string() << std::endl;
    }
    std::cout << "### ###\n\n";
}

void ReadClusteringEngine::export_spanning_tree(SpanningTree &tree, std::pair<std::vector<ReadID>, std::vector<ReadID>> &tails) {
    std::set<ComponentID> nodes;
    std::vector<std::string> node_strings;
    std::vector<std::string> edge_strings;
    for (auto edge : tree) {
        nodes.insert(edge.first);
        nodes.insert(edge.second);
        edge_strings.push_back(fmt::format("{{source: {}, target: {}}}", edge.first, edge.second));
    }

    std::set<ReadID> left_tail(tails.first.begin(), tails.first.end());
    std::set<ReadID> right_tail(tails.second.begin(), tails.second.end());
    for (ComponentID node : nodes) {
        char tail_marker;
        if (left_tail.contains(node)) {
            tail_marker = 'L';
        } else if (right_tail.contains(node)) {
            tail_marker = 'R';
        } else tail_marker = 'N';
        node_strings.push_back(fmt::format("{{ id: {}, start : {}, end : {}, tail: '{}' }}", node, read_intervals[node].first, read_intervals[node].second, tail_marker));
    }

    std::ofstream output;
    output.open(fmt::format("./data/trees/{}_tree_{}.js", reader->meta.filename, *nodes.begin()));
    output << fmt::format("nodes = [{}]\nedges = [{}]\n", boost::algorithm::join(node_strings, ", "), boost::algorithm::join(edge_strings, ", "));
    output.close();
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, ReadClusteringConfig config) {
    debug = read_iterator.categories > 1;
    this->reader = &read_iterator;
    this->config = config;
}

int ReadClusteringEngine::construct_indices(std::unordered_set<Kmer> &discriminative_kmers, int k) {
    reader->rewind();

    KmerIndex kmer_index;
    std::mutex index_merge_mut;
    auto construct_indices_thread = [this, &kmer_index, &index_merge_mut, &discriminative_kmers, k]() {
        std::optional<GenomeReadData> read;
        while ((read = reader->get_next_record()) != std::nullopt) {
            KmerIterator it = KmerIterator(read->sequence, k);
            std::vector<std::pair<Kmer, int>> in_read_discriminative;
            while (it.next_kmer()) {
                if (discriminative_kmers.contains(it.current_kmer)) {
                    in_read_discriminative.emplace_back(it.current_kmer, it.position_in_sequence);
                }
            }

            if (!in_read_discriminative.empty()) {
                std::vector<KmerID> in_read_discriminative_ids;
                ComponentID component_id = read->id;

                index_merge_mut.lock();
                kmer_positions.insert({read->id, {}});

                KmerID new_kmer_id = kmer_index.size();
                for (auto kmer_position_pair : in_read_discriminative) {
                    auto insert_result = kmer_index.insert({kmer_position_pair.first, new_kmer_id});
                    if (insert_result.second) {
                        kmer_component_index.push_back({});
                        ++new_kmer_id;
                    }
                    in_read_discriminative_ids.push_back(insert_result.first->second);
                    kmer_component_index[insert_result.first->second].push_back(component_id);
                    kmer_positions[read->id].insert({insert_result.first->second, kmer_position_pair.second});
                }
                index_merge_mut.unlock();

                std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
                ReadMetaData meta = {read->id, read->category_id, read->start, read->end};

                index_merge_mut.lock();
                component_index.insert({component_id, new ReadComponent(meta, in_read_discriminative_ids)});
                read_lengths.insert({read->id, read->sequence.length()});

                if (debug){
                    read_intervals[read->id] = {read->start, read->end};
                }
                index_merge_mut.unlock();
            }
        }
    };
    run_in_threads(construct_indices_thread);

    for (auto &arr : kmer_component_index) {
        std::sort(arr.begin(), arr.end());
    }
    return 0;
}

std::vector<ComponentConnection> ReadClusteringEngine::get_connections(std::vector<ComponentID> &component_ids, ConnectionScore min_score) {
    std::vector<ComponentConnection> connections;
    std::mutex connections_mut;

    auto queue = ConcurrentQueue<ComponentID>(component_ids.begin(), component_ids.end());
    auto get_connections_thread = [this, &queue, &connections_mut, &connections, min_score]() {
        std::optional<ComponentID> pivot;
        while ((pivot = queue.pop()) != std::nullopt) {
            ComponentID pivot_id = *pivot;

            tsl::robin_map<ComponentID, ConnectionScore> shared_kmer_counts;
            for (KmerID kmer_id : component_index[pivot_id]->discriminative_kmer_ids) {
                for (auto candidate : kmer_component_index[kmer_id]) {
                    shared_kmer_counts.insert({candidate, 0}).first.value()++;
                }
            }
            shared_kmer_counts.erase(pivot_id);

            for (auto id_count_pair : shared_kmer_counts){
                if (id_count_pair.second >= min_score){
                    bool is_good = debug ? component_index[id_count_pair.first]->categories == component_index[pivot_id]->categories : false;
                    connections_mut.lock();
                    connections.push_back({pivot_id, id_count_pair.first, id_count_pair.second, is_good});
                    connections_mut.unlock();
                }
            }
        }
    };
    run_in_threads(get_connections_thread);

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

std::vector<ComponentConnection> ReadClusteringEngine::get_all_connections() {
    std::vector<ComponentID> component_ids;
    for (auto id_ptr_pair : component_index) component_ids.push_back(id_ptr_pair.first);
    return get_connections(component_ids, 1);
}

std::vector<KmerID> ReadClusteringEngine::accumulate_kmer_ids(std::vector<ComponentID> &component_ids) {
    std::vector<std::vector<KmerID> *> merge_queue;
    for (auto component_id : component_ids) {
        merge_queue.push_back(&component_index[component_id]->discriminative_kmer_ids);
    }
    return merge_n_vectors(merge_queue, true);
}

std::vector<ComponentID> ReadClusteringEngine::merge_components(std::vector<Component> &components) {
    auto component_queue = ConcurrentQueue<Component>(components.begin(), components.end());
    IndexRemovalMap for_removal;
    std::vector<ComponentID> merged_ids;
    std::mutex component_merging_mut;

    auto merge_components_thread = [this, &merged_ids, &component_queue, &for_removal, &component_merging_mut]() {
        std::optional<Component> next_component;
        while ((next_component = component_queue.pop()) != std::nullopt) {
            auto component_ids = *next_component;

            if (component_ids.size() == 1) {
                component_merging_mut.lock();
                merged_ids.push_back(component_ids[0]);
                component_merging_mut.unlock();
                continue;
            }
            ReadComponent *survivor = component_index[component_ids[0]];

            std::set<CategoryID> survivor_categories;
            std::vector<ReadMetaData> survivor_contained_reads;

            for (auto component_id : component_ids) {
                auto component = component_index[component_id];
                survivor_contained_reads.insert(survivor_contained_reads.end(), component->contained_reads.begin(), component->contained_reads.end());
                component->contained_reads.clear();

                survivor_categories.insert(component->categories.begin(), component->categories.end());
            }

            survivor->discriminative_kmer_ids = accumulate_kmer_ids(component_ids);
            survivor->categories = survivor_categories;
            survivor->contained_reads = survivor_contained_reads;

            component_merging_mut.lock();
            merged_ids.push_back(component_ids[0]);
            for (auto component_id : component_ids) {
                for (KmerID kmer_id : component_index[component_id]->discriminative_kmer_ids) {
                    for_removal.insert({kmer_id, {}}).first.value().push_back(component_id);
                }
            }
            component_merging_mut.unlock();
        }
    };
    run_in_threads(merge_components_thread);

    auto removal_list_queue = ConcurrentQueue<IndexRemovalMap::value_type>(for_removal.begin(), for_removal.end());
    auto kmer_component_index_update = [this, &removal_list_queue]() {
        std::optional<IndexRemovalMap::value_type> next_removal;
        while ((next_removal = removal_list_queue.pop()) != std::nullopt) {
            KmerID kmer_id = next_removal->first;
            std::vector<ComponentID> *removal_list = &next_removal->second;
            std::sort(removal_list->begin(), removal_list->end());

            std::vector<ComponentID> updated_list;
            int i = 0, j = 0;
            while (i < removal_list->size() && j < kmer_component_index[kmer_id].size()) {
                if ((*removal_list)[i] < kmer_component_index[kmer_id][j]) {
                    i++;
                } else if (kmer_component_index[kmer_id][j] < (*removal_list)[i]) {
                    updated_list.push_back(kmer_component_index[kmer_id][j]);
                    j++;
                } else {
                    i++;
                    j++;
                }
            }
            kmer_component_index[kmer_id] = updated_list;
        }
    };
    run_in_threads(kmer_component_index_update);

    return merged_ids;
}

std::vector<std::pair<Component, SpanningTree>>
ReadClusteringEngine::union_find(std::vector<ComponentConnection> &connections, std::set<ComponentID> &restricted, int min_component_size) {
    tsl::robin_set<ComponentID> affected_vertices;
    for (auto &conn : connections) {
        affected_vertices.insert(conn.component_x_id);
        affected_vertices.insert(conn.component_y_id);
    }

    tsl::robin_map<ComponentID, ComponentID> parents;
    tsl::robin_map<ComponentID, Component> components;
    tsl::robin_map<ComponentID, std::vector<std::pair<ComponentID, ComponentID>>> component_edges;
    tsl::robin_map<ComponentID, bool> component_contains_restricted;
    for (ComponentID vertex : affected_vertices) {
        parents[vertex] = vertex;
        components[vertex] = {vertex};
        component_edges[vertex] = {};
        component_contains_restricted[vertex] = restricted.contains(vertex);
    }

    std::function<ComponentID(ComponentID)> get_parent = [&parents, &get_parent](ComponentID component_id) {
        if (component_id == parents[component_id]) {
            return component_id;
        }
        parents[component_id] = get_parent(parents[component_id]);
        return parents[component_id];
    };

    for (ComponentConnection &conn : connections) {
        ComponentID parent_x = get_parent(conn.component_x_id);
        ComponentID parent_y = get_parent(conn.component_y_id);
        if (parent_x == parent_y) continue;
        if (component_contains_restricted[parent_x] && component_contains_restricted[parent_y]) continue;

        ComponentID bigger, smaller;
        if (components[parent_x].size() > components[parent_y].size()) {
            bigger = parent_x;
            smaller = parent_y;
        } else {
            bigger = parent_y;
            smaller = parent_x;
        }

        for (ComponentID component_id : components[smaller]) {
            parents[component_id] = bigger;
        }
        components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
        components.erase(smaller);

        component_edges[bigger].emplace_back(conn.component_x_id, conn.component_y_id);
        component_edges[bigger].insert(component_edges[bigger].end(), component_edges[smaller].begin(), component_edges[smaller].end());
        component_edges.erase(smaller);

        component_contains_restricted[bigger] |= component_contains_restricted[smaller];
        component_contains_restricted.erase(smaller);
    }

    std::vector<std::pair<Component, SpanningTree>> result;
    for (const auto &comp : components) {
        if (comp.second.size() >= min_component_size) {
            result.emplace_back(comp.second, component_edges[comp.first]);
        }
    }
    return result;
}

void ReadClusteringEngine::export_components(std::vector<ComponentID> &component_ids, std::string &directory_path) {
    std::filesystem::remove_all(directory_path);
    std::filesystem::create_directories(directory_path);

    std::map<ComponentID, std::ofstream> export_files;
    tsl::robin_map<ReadID, ComponentID> read_id_to_component;
    for (auto component_id : component_ids){
        export_files[component_id] = std::ofstream(fmt::format("{}/#{}.fa", directory_path, component_id));
        for (auto read_meta : component_index[component_id]->contained_reads){
            read_id_to_component[read_meta.id] = component_id;
        }
    }

    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        export_files[read_id_to_component[read->id]] << read->fastX_string() << std::endl;
    }
    for (auto it = begin(export_files); it != end(export_files); ++it){
        it->second.close();
    }
    std::cout << fmt::format("Exported {} components\n", component_ids.size());
}

int ReadClusteringEngine::approximate_read_overlap(ReadID x, ReadID y) {
    auto shared_kmers = get_vectors_intersection(component_index[x]->discriminative_kmer_ids, component_index[y]->discriminative_kmer_ids);

    std::vector<int> x_positions, y_positions;
    for (auto kmer_id : shared_kmers) {
        x_positions.push_back(kmer_positions[x][kmer_id]);
        y_positions.push_back(kmer_positions[y][kmer_id]);
    }
    int max_x, max_y, min_x, min_y;
    max_x = *std::max_element(x_positions.begin(), x_positions.end());
    max_y = *std::max_element(y_positions.begin(), y_positions.end());
    min_x = *std::min_element(x_positions.begin(), x_positions.end());
    min_y = *std::min_element(y_positions.begin(), y_positions.end());

    auto overlap = std::max(max_x - min_x, max_y - min_y);
    return overlap;
}

std::pair<std::vector<ReadID>, std::vector<ReadID>> ReadClusteringEngine::get_spanning_tree_tails(SpanningTree &tree) {
    tsl::robin_map<ReadID, tsl::robin_map<ReadID, int>> adjacency_map;
    for (auto edge : tree) {
        auto dist = approximate_read_overlap(edge.first, edge.second);
        adjacency_map.insert({edge.first, {}}).first.value().insert({edge.second, dist});
        adjacency_map.insert({edge.second, {}}).first.value().insert({edge.first, dist});
    }

    auto distance_bfs = [this, &adjacency_map](ReadID starting_vertex) {
        std::stack<ReadID> stack;
        stack.push(starting_vertex);

        std::queue<ReadID> bfs_queue;
        bfs_queue.push(starting_vertex);
        tsl::robin_set<ReadID> visited;
        tsl::robin_map<ReadID, uint64_t> distances = {{starting_vertex, read_lengths[starting_vertex]}};
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

    auto get_max_distance_pair = [](tsl::robin_map<ReadID, uint64_t> &distances) {
        return *std::max_element(distances.begin(), distances.end(), [](const std::pair<ReadID, uint64_t> &p1, const std::pair<ReadID, uint64_t> &p2) {
            return p1.second < p2.second;
        });
    };

    auto initial_distances = distance_bfs(adjacency_map.begin()->first);
    auto farthest_vertex_and_dist = get_max_distance_pair(initial_distances);

    std::vector<ReadID> left_side_vertices, right_side_vertices;

    auto distances_to_right = distance_bfs(farthest_vertex_and_dist.first);
    auto farthest_right_vertex_and_dist = get_max_distance_pair(distances_to_right);
    uint64_t tail_length = reader->meta.avg_read_length * 3;
    for (auto vertex_dist_pair : distances_to_right) {
        if (vertex_dist_pair.second + tail_length > farthest_right_vertex_and_dist.second) {
            right_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    auto distances_to_left = distance_bfs(farthest_right_vertex_and_dist.first);
    auto farthest_left_vertex_and_dist = get_max_distance_pair(distances_to_left);
    for (auto vertex_dist_pair : distances_to_left) {
        if (vertex_dist_pair.second + tail_length > farthest_left_vertex_and_dist.second) {
            left_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    if (debug){
        std::cout << fmt::format("Initial left {} - {}\n", read_intervals[farthest_left_vertex_and_dist.first].first, read_intervals[farthest_left_vertex_and_dist.first].second);
        std::cout << fmt::format("Initial right {} - {}\n", read_intervals[farthest_right_vertex_and_dist.first].first, read_intervals[farthest_right_vertex_and_dist.first].second);
    }


    return {left_side_vertices, right_side_vertices};
}

std::vector<ComponentID> ReadClusteringEngine::amplify_component(Component &component, int min_score) {
    auto connections = get_connections(component, min_score);
    std::set<ReadID> id_set(component.begin(), component.end());
    for (auto conn : connections) {
        id_set.insert(conn.component_x_id);
        id_set.insert(conn.component_y_id);
    }
    std::vector<ComponentID> result(id_set.begin(), id_set.end());
    return result;
}

std::vector<ComponentConnection> ReadClusteringEngine::get_core_component_connections(std::vector<std::pair<Component, SpanningTree>> &components_with_trees) {
    std::map<ComponentID, std::pair<std::vector<KmerID>, std::vector<KmerID>>> core_component_tails;
    for (auto &component_and_tree : components_with_trees) {
        auto tails = get_spanning_tree_tails(component_and_tree.second);
        auto left_vertices = amplify_component(tails.first, config.tail_amplification_min_score);
        auto right_vertices = amplify_component(tails.second, config.tail_amplification_min_score);

        auto left_kmers = accumulate_kmer_ids(left_vertices);
        auto right_kmers = accumulate_kmer_ids(right_vertices);
        auto component_kmers = accumulate_kmer_ids(component_and_tree.first);

        core_component_tails.insert({component_and_tree.first[0], {left_kmers, right_kmers}});

        if (debug){
            export_spanning_tree(component_and_tree.second, tails);

            auto left_tail_intervals = get_read_coverage_intervals(tails.first);
            for (auto i : left_tail_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << std::endl;
            auto right_tail_intervals = get_read_coverage_intervals(tails.second);
            for (auto i : right_tail_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << std::endl;
            auto component_intervals = get_read_coverage_intervals(component_and_tree.first);
            for (auto i : component_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << "\n\n";
        }
    }

    std::vector<ComponentConnection> tail_connecting_edges;
    for (auto c_id_kmers_pair : core_component_tails) {
        for (auto c_id_kmers_pair_2 : core_component_tails) {
            if (c_id_kmers_pair.first < c_id_kmers_pair_2.first) {
                auto left_1_left_2 = get_vectors_intersection(c_id_kmers_pair.second.first, c_id_kmers_pair_2.second.first).size();
                auto left_1_right_2 = get_vectors_intersection(c_id_kmers_pair.second.first, c_id_kmers_pair_2.second.second).size();
                auto right_1_left_1 = get_vectors_intersection(c_id_kmers_pair.second.second, c_id_kmers_pair_2.second.first).size();
                auto right_1_right_2 = get_vectors_intersection(c_id_kmers_pair.second.second, c_id_kmers_pair_2.second.second).size();

                std::vector<ConnectionScore> scores = {left_1_left_2, left_1_right_2, right_1_left_1, right_1_right_2};
                auto best_score = *std::max_element(scores.begin(), scores.end());

                bool is_good = debug ? component_index[c_id_kmers_pair.first]->categories == component_index[c_id_kmers_pair_2.first]->categories : false;
                tail_connecting_edges.push_back({c_id_kmers_pair.first, c_id_kmers_pair_2.first, best_score, is_good});
            }
        }
    }
    std::sort(tail_connecting_edges.rbegin(), tail_connecting_edges.rend());

    if (debug){
        for (auto edge : tail_connecting_edges) {
            std::cout << fmt::format(
                    "{} : {} between {} and {}\n",
                    edge.is_good ? '+' : '-',
                    edge.score,
                    component_index[edge.component_x_id]->to_string(),
                    component_index[edge.component_y_id]->to_string()
                    );

            if (edge.score == 0) break;
        }
    }

    return filter_connections(tail_connecting_edges, [](ComponentConnection &conn) { return conn.score > 0; });
}

std::vector<Component> spectral_clustering(std::vector<ComponentConnection> &connections){
    std::unordered_map<ComponentID, unsigned int> component_to_id;
    std::unordered_map<unsigned int, ComponentID> id_to_component;
    for (auto &conn : connections){
        auto id = component_to_id.insert({conn.component_x_id, component_to_id.size()}).first->second;
        id_to_component[id] = conn.component_x_id;

        id = component_to_id.insert({conn.component_y_id, component_to_id.size()}).first->second;
        id_to_component[id] = conn.component_y_id;
    }

    unsigned int size = component_to_id.size();
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size,size);

    for (auto &conn : connections){
        auto x = component_to_id[conn.component_x_id];
        auto y = component_to_id[conn.component_y_id];
        m(x, y) = exp(sqrt(conn.score));
        m(y, x) = exp(sqrt(conn.score));
    }
    SpectralClustering c = SpectralClustering(m, size);
    auto clusters = c.clusterRotate();

    std::vector<Component> result(clusters.size(), Component());
    for (unsigned int i=0; i < clusters.size(); i++) {
        for (auto id : clusters[i]){
            result[i].push_back(id_to_component[id]);
        }
    }
    return result;
}

std::vector<ComponentID> ReadClusteringEngine::run_clustering(std::unordered_set<Kmer> &discriminative_kmers, int k) {
    auto extract_components = [](std::vector<std::pair<Component, SpanningTree>> &union_find_output) {
        std::vector<Component> result;
        std::transform(union_find_output.begin(), union_find_output.end(), std::back_inserter(result),
                       [](std::pair<Component, SpanningTree> &p) { return p.first; }
        );
        return result;
    };

    auto filter_components = [this](const std::function<bool(ReadComponent *)> &func) {
        std::vector<ComponentID> result;
        for (auto component_id_and_ptr : component_index) {
            if (func(component_id_and_ptr.second)) {
                result.push_back(component_id_and_ptr.first);
            }
        }
        return result;
    };

    auto remove_merged_components = [this](){
        for(auto it = component_index.begin(); it != component_index.end();){
            if (it->second->size() == 0){
                free(it->second);
                component_index.erase(it++);
            } else ++it;
        }
    };

    auto get_core_component_ids = [this](int threshold_size) {
        std::vector<ComponentID> result;
        for (auto id_ptr_pair : component_index) {
            if (id_ptr_pair.second->size() >= threshold_size) {
                result.push_back(id_ptr_pair.first);
            }
        }
        return result;
    };

    if (debug){
        timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")(discriminative_kmers, k);
    } else {
        construct_indices(discriminative_kmers, k);
    }

    auto component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_all_connections, this, "Strong connections")();
    std::vector<ComponentConnection> strong_connections(
            component_connections.begin(),
            component_connections.begin() + component_connections.size() * config.core_forming_connections
            );

    if (debug){
        plot_connection_quality(component_connections);
        auto threshold_value = std::find_if(component_connections.begin(), component_connections.end(), [](ComponentConnection &conn) { return !conn.is_good; })->score + 1;
        std::cout << fmt::format("Threshold value is {}\n", threshold_value);
        std::cout << fmt::format("Selected threshold is {}\n", strong_connections[strong_connections.size() - 1].score);
    }

    std::set<ComponentID> restricted;
    auto components_and_trees = union_find(strong_connections, restricted, config.core_component_min_size);
    std::vector<Component> components = extract_components(components_and_trees);
    auto core_component_ids = timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Creation of strong clusters")(components);
    print_components(core_component_ids);

    if (core_component_ids.size() > 2){
        auto core_component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_core_component_connections, this, "Core connections")(components_and_trees);
        std::vector<ComponentConnection> strong_core_connections(core_component_connections.begin(), core_component_connections.begin() + components_and_trees.size() - 1);
        if (!strong_core_connections.empty()) {
            auto spectral_components = spectral_clustering(strong_core_connections);
            timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Creation of spectral")(spectral_components);
        }
        remove_merged_components();
        core_component_ids = get_core_component_ids(config.core_component_min_size);
    }

    print_components(core_component_ids);

    {
        component_connections = get_connections(core_component_ids, config.enrichment_connections_min_score);
        restricted.insert(core_component_ids.begin(), core_component_ids.end());
        components_and_trees = union_find(component_connections, restricted, 2);
        components = extract_components(components_and_trees);
        merge_components(components);
        remove_merged_components();
        core_component_ids = get_core_component_ids(config.core_component_min_size);
    }

    print_components(core_component_ids);

    return core_component_ids;
}
