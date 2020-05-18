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
    run_command_with_input("python scripts/plotting.py --plot connection_histogram", hist_input);
}


void ReadClusteringEngine::plot_read_coverage(std::vector<ComponentID> &components){
    auto get_coverage_runs = [this](std::vector<ReadID> &reads){
        std::vector<Endpoint> endpoints;
        for (auto read_id : reads) {
            endpoints.emplace_back(read_metas[read_id].start, true);
            endpoints.emplace_back(read_metas[read_id].end, false);
        }
        std::sort(endpoints.begin(), endpoints.end());

        std::map<int, int> result;
        int current_coverage = 0;
        int prev_pos = -1;
        for (auto endpoint : endpoints) {
            if (endpoint.second) {
                current_coverage++;
                result.insert({current_coverage, 0}).first->second += (endpoint.first - prev_pos);
                prev_pos = endpoint.first;
            } else {
                result.insert({current_coverage, 0}).first->second += (endpoint.first - prev_pos);
                prev_pos = endpoint.first;
                current_coverage--;
            }
        }
        return result;
    };

    std::vector<std::string> component_strings;
    for (const auto & component_id : components){
        std::vector<std::string> coverage_strings;
        std::map<int, int> coverage_runs = get_coverage_runs(component_index[component_id]->contained_read_ids);
        for (auto cov_bases : coverage_runs){
            coverage_strings.push_back(fmt::format("({}, {})", cov_bases.first, cov_bases.second));
        }
        component_strings.push_back(fmt::format("[{}]", boost::algorithm::join(coverage_strings, ", ")));
    }
    std::string input = fmt::format("[{}]", boost::algorithm::join(component_strings, ", "));
    run_command_with_input("python scripts/plotting.py --plot component_coverage", input);
}


void plot_spectral(std::vector<ComponentConnection> &connections, std::vector<Component> &result, ComponentIndex &index) {
    std::string input;
    std::vector<std::string> conn_strings;
    std::vector<std::string> category_strings;
    for (const auto &conn : connections) {
        conn_strings.push_back(fmt::format("({}, {}, {})", conn.component_x_id, conn.component_y_id, conn.score));
        category_strings.push_back(fmt::format("({}, {})", conn.component_x_id, *index[conn.component_x_id]->categories.begin()));
        category_strings.push_back(fmt::format("({}, {})", conn.component_y_id, *index[conn.component_y_id]->categories.begin()));
    }
    input += fmt::format("[{}]\n", boost::algorithm::join(conn_strings, ", "));
    input += fmt::format("[{}]\n", boost::algorithm::join(category_strings, ", "));

    std::vector<std::string> component_strings;
    for (const auto & comp : result){
        std::vector<std::string> temp;
        for (auto id : comp){
            temp.push_back(fmt::format("{}", id));
        }
        component_strings.push_back("[" + boost::algorithm::join(temp, ", ") + "]");
    }
    input += fmt::format("[{}]\n", boost::algorithm::join(component_strings, ", "));
    run_command_with_input("python scripts/plotting.py --plot spectral_graph", input);
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
        auto a = std::make_pair(read_metas[conn.component_x_id].start, read_metas[conn.component_x_id].end);
        auto b = std::make_pair(read_metas[conn.component_y_id].start, read_metas[conn.component_y_id].end);

        int overlap = std::max(0u, std::min(a.second, b.second) - std::max(a.first, b.first));
        if (overlap < 1) continue;

        if (!overlaps_for_strength.contains(conn.score)) overlaps_for_strength[conn.score] = {};
        if (!overlaps_for_strength[conn.score].contains(overlap)) overlaps_for_strength[conn.score][overlap] = 0;
        overlaps_for_strength[conn.score][overlap] += 1;
    }

    std::vector<std::string> per_score_overlaps;
    for (const auto &score_overlaps : overlaps_for_strength) {
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
std::vector<Interval> ReadClusteringEngine::get_read_coverage_intervals(std::vector<ReadID> &read_ids) {
    std::vector<Interval> result;

    std::vector<Endpoint> endpoints;
    for (auto read_id : read_ids) {
        endpoints.emplace_back(read_metas[read_id].start, true);
        endpoints.emplace_back(read_metas[read_id].end, false);
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
        std::cout << component_index[id]->to_string(read_metas) << std::endl;
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
        node_strings.push_back(fmt::format("{{ id: {}, start : {}, end : {}, tail: '{}' }}", node, read_metas[node].start, read_metas[node].end, tail_marker));
    }

    std::ofstream output;
    output.open(fmt::format("./data/trees/{}_tree_{}.js", reader->meta.filename, *nodes.begin()));
    output << fmt::format("nodes = [{}]\nedges = [{}]\n", boost::algorithm::join(node_strings, ", "), boost::algorithm::join(edge_strings, ", "));
    output.close();
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, ReadClusteringConfig config) {
    debug = read_iterator.file_meta.size() == read_iterator.categories;
    this->reader = &read_iterator;
    this->config = config;
}

int ReadClusteringEngine::construct_indices(std::unordered_set<Kmer> &discriminative_kmers, int k) {
    reader->rewind();

    KmerIndex kmer_index;
    KmerID kmer_id = 0;
    for (auto kmer : discriminative_kmers){
        kmer_index[kmer] = kmer_id++;
    }
    kmer_component_index.resize(kmer_index.size(), {});

    std::mutex comp_index_mut, kmer_comp_index_mut;
    auto construct_indices_thread = [this, &kmer_index, &comp_index_mut, &kmer_comp_index_mut, &discriminative_kmers, k]() {
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

                ReadMetaData meta = {read->id, (uint32_t) read->sequence.length(), read->category_id, read->start, read->end, {}};

                for (auto kmer_position_pair : in_read_discriminative) {
                    auto kmer_id = kmer_index[kmer_position_pair.first];
                    in_read_discriminative_ids.push_back(kmer_id);
                    kmer_comp_index_mut.lock();
                    kmer_component_index[kmer_id].push_back(component_id);
                    meta.kmer_positions.insert({kmer_id, kmer_position_pair.second});
                    kmer_comp_index_mut.unlock();
                }

                std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());

                comp_index_mut.lock();
                read_metas[read->id] = meta;
                component_index.insert({component_id, new ReadComponent(meta, in_read_discriminative_ids)});
                comp_index_mut.unlock();
            }
        }
    };
    run_in_threads(construct_indices_thread, config.threads);

    for (auto &arr : kmer_component_index) {
        std::sort(arr.begin(), arr.end());
    }
    if (debug) {
        uint32_t discriminative = 0;
        uint32_t total = 0;
        for (const auto &occurrences : kmer_component_index) {
            std::set<CategoryID> categories;
            for (auto component_id : occurrences) {
                categories.insert(read_metas[component_id].category_id);
            }
            if (categories.size() == 1) discriminative++;
            if (!categories.empty()) total++;
        }
        std::cout << fmt::format("{} out of {} kmers are discriminative \n", discriminative, total);
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

            for (auto id_count_pair : shared_kmer_counts) {
                if (id_count_pair.second >= min_score) {
                    bool is_good = debug ? component_index[id_count_pair.first]->categories == component_index[pivot_id]->categories : false;
                    connections_mut.lock();
                    connections.push_back({pivot_id, id_count_pair.first, id_count_pair.second, is_good});
                    connections_mut.unlock();
                }
            }
        }
    };
    run_in_threads(get_connections_thread, config.threads);

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

std::vector<ComponentConnection> ReadClusteringEngine::get_all_connections(ConnectionScore min_score) {
    std::vector<ComponentID> component_ids;
    for (auto id_ptr_pair : component_index) component_ids.push_back(id_ptr_pair.first);
    return get_connections(component_ids, min_score);
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
            std::vector<ReadID> survivor_contained_reads;

            for (auto component_id : component_ids) {
                auto component = component_index[component_id];
                survivor_contained_reads.insert(survivor_contained_reads.end(), component->contained_read_ids.begin(), component->contained_read_ids.end());
                component->contained_read_ids.clear();

                survivor_categories.insert(component->categories.begin(), component->categories.end());
            }

            survivor->discriminative_kmer_ids = accumulate_kmer_ids(component_ids);
            survivor->categories = survivor_categories;
            survivor->contained_read_ids = survivor_contained_reads;

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
    run_in_threads(merge_components_thread, config.threads);

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
    run_in_threads(kmer_component_index_update, config.threads);

    return merged_ids;
}

std::vector<std::pair<Component, SpanningTree>>
union_find(std::vector<ComponentConnection> &connections, std::set<ComponentID> &restricted, int min_component_size, int max_component_size) {
    if (max_component_size == -1) max_component_size = INT_MAX;
    tsl::robin_set<ComponentID> affected_vertices;
    for (auto &conn : connections) {
        affected_vertices.insert(conn.component_x_id);
        affected_vertices.insert(conn.component_y_id);
    }

    tsl::robin_map<ComponentID, ComponentID> parents;
    tsl::robin_map<ComponentID, Component> components;
    tsl::robin_map<ComponentID, bool> component_contains_restricted;
    tsl::robin_map<ComponentID, std::vector<std::pair<ComponentID, ComponentID>>> spanning_trees;
    for (ComponentID vertex : affected_vertices) {
        parents[vertex] = vertex;
        components[vertex] = {vertex};
        spanning_trees[vertex] = {};
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
        if (components[parent_x].size() + components[parent_y].size() > max_component_size) continue;

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

        spanning_trees[bigger].emplace_back(conn.component_x_id, conn.component_y_id);
        spanning_trees[bigger].insert(spanning_trees[bigger].end(), spanning_trees[smaller].begin(), spanning_trees[smaller].end());
        spanning_trees.erase(smaller);

        component_contains_restricted[bigger] |= component_contains_restricted[smaller];
        component_contains_restricted.erase(smaller);
    }

    std::vector<std::pair<Component, SpanningTree>> result;
    for (const auto &comp : components) {
        if (comp.second.size() >= min_component_size) {
            result.emplace_back(comp.second, spanning_trees[comp.first]);
        }
    }
    return result;
}

int ReadClusteringEngine::approximate_read_overlap(ReadID x, ReadID y) {
    auto shared_kmers = get_vectors_intersection(component_index[x]->discriminative_kmer_ids, component_index[y]->discriminative_kmer_ids);

    std::vector<int> x_positions, y_positions;
    for (auto kmer_id : shared_kmers) {
        x_positions.push_back(read_metas[x].kmer_positions[kmer_id]);
        y_positions.push_back(read_metas[y].kmer_positions[kmer_id]);
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
        std::queue<ReadID> bfs_queue;
        bfs_queue.push(starting_vertex);
        tsl::robin_set<ReadID> visited;
        tsl::robin_map<ReadID, uint64_t> distances = {{starting_vertex, read_metas[starting_vertex].length}};
        while (!bfs_queue.empty()) {
            auto vertex = bfs_queue.front();
            visited.insert(vertex);
            bfs_queue.pop();
            if (adjacency_map.contains(vertex)) {
                for (auto adjacency : adjacency_map[vertex]) {
                    if (!visited.contains(adjacency.first)) {
                        distances[adjacency.first] = distances[vertex] + read_metas[adjacency.first].length - adjacency.second;
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
    uint64_t tail_length = reader->meta.avg_read_length * 2;

    auto distances_to_right = distance_bfs(farthest_vertex_and_dist.first);
    auto farthest_right_vertex_and_dist = get_max_distance_pair(distances_to_right);
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

    if (debug) {
        std::cout << fmt::format("Initial left: {} - {}\n", read_metas[farthest_left_vertex_and_dist.first].start, read_metas[farthest_left_vertex_and_dist.first].end);
        std::cout << fmt::format("Initial right: {} - {}\n", read_metas[farthest_right_vertex_and_dist.first].start, read_metas[farthest_right_vertex_and_dist.first].end);
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

        if (debug) {
            export_spanning_tree(component_and_tree.second, tails);

            auto left_tail_intervals = get_read_coverage_intervals(tails.first);
            std::cout << "Left tail: ";
            for (auto i : left_tail_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second);
            std::cout << std::endl;
            auto right_tail_intervals = get_read_coverage_intervals(tails.second);
            std::cout << "Right tail: ";
            for (auto i : right_tail_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second);
            std::cout << std::endl;
            auto component_intervals = get_read_coverage_intervals(component_and_tree.first);
            std::cout << "Component: ";
            for (auto i : component_intervals) std::cout << fmt::format("{}-{} , ", i.first, i.second);
            std::cout << "\n\n";
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

    if (debug) {
        for (auto edge : tail_connecting_edges) {
            std::cout << fmt::format(
                    "{} : {} between {} and {}\n",
                    edge.is_good ? '+' : '-',
                    edge.score,
                    component_index[edge.component_x_id]->to_string(read_metas),
                    component_index[edge.component_y_id]->to_string(read_metas)
            );

            if (edge.score == 0) break;
        }
    }

    return filter_connections(tail_connecting_edges, [](ComponentConnection &conn) { return conn.score > 0; });
}

std::vector<Component> spectral_clustering(std::vector<ComponentConnection> &connections, int dims) {
    std::unordered_map<ComponentID, unsigned int> component_to_id;
    std::unordered_map<unsigned int, ComponentID> id_to_component;
    std::vector<ConnectionScore> scores;
    for (auto &conn : connections) {
        auto id = component_to_id.insert({conn.component_x_id, component_to_id.size()}).first->second;
        id_to_component[id] = conn.component_x_id;

        id = component_to_id.insert({conn.component_y_id, component_to_id.size()}).first->second;
        id_to_component[id] = conn.component_y_id;

        scores.push_back(conn.score);
    }

    int max_exponent = 20;

    int size = component_to_id.size();
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size, size);

    ConnectionScore max_score = *std::max_element(scores.begin(), scores.end());
    ConnectionScore min_score = *std::min_element(scores.begin(), scores.end());
    auto scale_strength = [max_exponent, max_score, min_score](ConnectionScore score){
        return ((double)(max_exponent - 0.3)*(double)(score - min_score)) / (double)(max_score - min_score) + 0.3;
    };

    for (auto &conn : connections) {
        auto x = component_to_id[conn.component_x_id];
        auto y = component_to_id[conn.component_y_id];
        double scaled_score = scale_strength(conn.score);
        m(x, y) = exp(scaled_score);
        m(y, x) = exp(scaled_score);
    }
    dims = std::min(size, dims);

    auto c = SpectralClustering(m, dims);
    auto clusters = c.clusterRotate();

    std::vector<Component> result(clusters.size(), Component());
    for (unsigned int i = 0; i < clusters.size(); i++) {
        for (auto id : clusters[i]) {
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

    auto remove_merged_components = [this]() {
        for (auto it = component_index.begin(); it != component_index.end();) {
            if (it->second->size() == 0) {
                free(it->second);
                component_index.erase(it++);
            } else ++it;
        }
    };

    auto get_component_ids = [this](int threshold_size) {
        std::vector<ComponentID> result;
        for (auto id_ptr_pair : component_index) {
            if (id_ptr_pair.second->size() >= threshold_size) {
                result.push_back(id_ptr_pair.first);
            }
        }
        return result;
    };

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Index construction")(discriminative_kmers, k);

    if (config.force_spectral){
        auto component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_all_connections, this, "Calculation of connections between reads")(5);
        auto spectral = timeMeasure(spectral_clustering, "Forced spectral clustering")(component_connections, config.spectral_dims);
        merge_components(spectral);
        auto ids = get_component_ids(config.scaffold_component_min_size);
        print_components(ids);
        return ids;
    }

    std::vector<ComponentConnection> scaffold_forming_connections, component_connections;
    if (config.scaffold_forming_score > 0){
        auto component_ids = filter_components([this](ReadComponent* c){ return c->discriminative_kmer_ids.size() >= config.scaffold_forming_score; });
        component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Calculation of connections between reads")(component_ids, config.scaffold_forming_score);
        scaffold_forming_connections = filter_connections(component_connections, [this](ComponentConnection &conn){ return conn.score > config.scaffold_forming_score; });
    } else {
        component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_all_connections, this, "Calculation of connections between reads")(1);
        scaffold_forming_connections = std::vector<ComponentConnection>(component_connections.begin(), component_connections.begin() + component_connections.size() * config.scaffold_forming_fraction);
    }

    if (debug) {
        plot_connection_quality(component_connections);
    }

    std::set<ComponentID> restricted;
    auto components_and_trees = timeMeasure(union_find, "Union-find")(scaffold_forming_connections, restricted, config.scaffold_component_min_size, config.scaffold_component_max_size);
    std::vector<Component> scaffold_components = extract_components(components_and_trees);
    auto scaffold_component_ids = timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Merging of initial components")(scaffold_components);
    print_components(scaffold_component_ids);

    if (scaffold_component_ids.size() > 2) {
        auto core_forming_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_core_component_connections, this, "Calculation of tail connections")(components_and_trees);
        auto strong_core_connections = filter_connections(core_forming_connections, [](ComponentConnection &conn) { return conn.score > 5; });
        if (!strong_core_connections.empty()) {
            auto spectral_components = timeMeasure(spectral_clustering, "Spectral clustering")(strong_core_connections, config.spectral_dims);
            if (debug) plot_spectral(strong_core_connections, spectral_components, component_index);
            timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Merging of scaffold components")(spectral_components);
        }
        remove_merged_components();
    }

    auto core_component_ids = get_component_ids(config.scaffold_component_min_size);
    if (debug){
        print_components(core_component_ids);
        plot_read_coverage(core_component_ids);
    }

    {
        component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Calculation of enrichment connections")(core_component_ids, config.enrichment_connections_min_score);
        if (debug) plot_connection_quality(component_connections, core_component_ids);
        restricted.insert(core_component_ids.begin(), core_component_ids.end());
        components_and_trees = union_find(component_connections, restricted, 2, -1);
        scaffold_components = extract_components(components_and_trees);
        timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Merging into core components")(scaffold_components);
        remove_merged_components();
        core_component_ids = get_component_ids(config.scaffold_component_min_size);
    }

    if (debug){
        print_components(core_component_ids);
        plot_read_coverage(core_component_ids);
    }

    return core_component_ids;
}

void ReadClusteringEngine::export_components(std::vector<ComponentID> &component_ids, std::string &directory_path) {
    std::filesystem::remove_all(directory_path);
    std::filesystem::create_directories(directory_path);

    std::map<ComponentID, std::ofstream> export_files;
    tsl::robin_map<ReadID, ComponentID> read_id_to_component;
    for (auto component_id : component_ids) {
        export_files[component_id] = std::ofstream(fmt::format("{}/#{}.fa", directory_path, component_id));
        for (auto read_id : component_index[component_id]->contained_read_ids) {
            read_id_to_component[read_id] = component_id;
        }
    }

    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        export_files[read_id_to_component[read->id]] << read->fastX_string() << std::endl;
    }
    for (auto it = begin(export_files); it != end(export_files); ++it) {
        it->second.close();
    }
    std::cout << fmt::format("Exported {} components\n", component_ids.size());
}