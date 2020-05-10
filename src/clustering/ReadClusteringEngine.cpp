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

//DEBUG
void plot_cluster_coverage(std::vector<ReadComponent *> &clusters) {
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
    sort(ids.rbegin(), ids.rend(), [this](ComponentID x, ComponentID y) -> bool { return component_index[x]->size() < component_index[y]->size(); });
    std::cout << "Clusters:\n";
    for (auto id : ids) {
        std::cout << component_index[id]->to_string() << std::endl;
    }
    std::cout << "End clusters\n";
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

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator) {
    this->reader = &read_iterator;
}

void ReadClusteringEngine::set_kmers(tsl::robin_set<Kmer> &_kmers, int _k) {
    this->k = _k;
    this->kmers = new bloom::BloomFilter<Kmer>(_kmers.size(), 0.01);
    for (Kmer kmer : _kmers) {
        this->kmers->add(kmer);
    }
}

void ReadClusteringEngine::set_kmers(bloom::BloomFilter<Kmer> *_kmers, int _k) {
    this->k = _k;
    this->kmers = _kmers;
}

int ReadClusteringEngine::construct_indices() {
    reader->rewind();

    KmerIndex kmer_index;
    std::mutex index_merge_mut;
    auto construct_indices_thread = [this, &kmer_index, &index_merge_mut]() {
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
                ComponentID component_id = read->id;

                index_merge_mut.lock();
                kmer_positions.insert({read->id, {}});

                std::pair<KmerIndex::iterator, bool> insert_result;
                KmerID new_kmer_id = kmer_index.size();
                for (auto kmer_position_pair : in_read_discriminative) {
                    insert_result = kmer_index.insert(KmerIndex::value_type(kmer_position_pair.first, new_kmer_id));
                    if (insert_result.second) {
                        kmer_component_index.push_back({});
                        ++new_kmer_id;
                    }
                    in_read_discriminative_ids.push_back(insert_result.first->second);
                    kmer_component_index[insert_result.first->second].push_back(component_id);
                    kmer_positions[read->id].insert({insert_result.first->second, kmer_position_pair.second});
                }

                std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
                ReadMetaData meta = {read->id, read->category_id, read->start, read->end};
                component_index.insert(ComponentIndex::value_type(component_id, new ReadComponent(meta, in_read_discriminative_ids)));
                read_intervals[read->id] = {read->start, read->end};

                read_lengths.insert({read->id, read->sequence.length()});
                index_merge_mut.unlock();
            }
        }
    };
    run_in_threads(construct_indices_thread);

    for (auto &arr : kmer_component_index) {
        std::sort(arr.begin(), arr.end());
    }
    std::cout << fmt::format("Kmer cluster index size {}\n", kmer_component_index.size());
    return 0;
}

std::vector<ComponentConnection> ReadClusteringEngine::get_connections(std::vector<ComponentID> &component_ids, ConnectionScore min_score) {
    std::vector<ComponentConnection> connections;
    std::mutex connections_mut;

    auto queue = ConcurrentQueue<ComponentID>(component_ids.begin(), component_ids.end());
    auto get_connections_thread = [this, &queue, &connections_mut, &connections, min_score]() {
        while (true) {
            auto pivot = queue.pop();
            if (pivot == std::nullopt) break;
            ComponentID pivot_id = *pivot;

            std::map<ComponentID, ConnectionScore> shared_kmer_counts;
            for (KmerID kmer_id : component_index[pivot_id]->discriminative_kmer_ids) {
                for (auto candidate : kmer_component_index[kmer_id]) {
                    shared_kmer_counts.insert(std::map<ComponentID, ConnectionScore>::value_type(candidate, 0)).first->second += 1;
                }
            }
            shared_kmer_counts.erase(pivot_id);

            for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
                if (iter->second < min_score) continue;
                connections_mut.lock();
                connections.push_back({pivot_id, iter->first, iter->second,
                                       component_index[iter->first]->categories == component_index[pivot_id]->categories && component_index[iter->first]->categories.size() == 1});
                connections_mut.unlock();
            }
        }
    };
    run_in_threads(get_connections_thread);

    std::sort(connections.rbegin(), connections.rend());
    return connections;
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

std::map<ComponentID, std::string> ReadClusteringEngine::export_components(std::vector<ComponentID> &component_ids, fs::path &directory_path) {
    tsl::robin_map<ReadID, GenomeReadData> id_to_read;

    fs::remove_all(directory_path);
    fs::create_directories(directory_path);
    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        id_to_read[read->id] = *read;
    }

    std::map<ComponentID, std::string> mapping;
    for (auto id : component_ids) {
        std::ofstream component_file;
        std::string component_file_path = fmt::format("{}/#{}.fa", directory_path.string(), id);

        component_file.open(component_file_path);
        for (auto read_meta : component_index[id]->contained_reads) {
            component_file << id_to_read[read_meta.id].fasta_string() << std::endl;
        }

        component_file.close();
        mapping[id] = component_file_path;
    }

    return mapping;
}

std::map<ComponentID, std::vector<std::string>> ReadClusteringEngine::assemble_clusters(std::vector<ComponentID> &component_ids) {
    sort(component_ids.rbegin(), component_ids.rend(), [this](ComponentID x, ComponentID y) -> bool { return component_index[x]->size() < component_index[y]->size(); });

    fs::path export_path = "./data/assembly_" + reader->meta.filename;
    auto mapping = export_components(component_ids, export_path);

    std::map<ComponentID, std::vector<std::string>> assembly_mapping;
    for (ComponentID id : component_ids) {
        std::string assembled_fasta_path = mapping[id] + "_assembled.fa";

        std::string cmd = fmt::format("./scripts/assemble_cluster.sh {} {}", mapping[id], assembled_fasta_path);
        run_command_with_input(cmd.c_str(), "");

        std::cout << component_index[id]->to_string() << " : \n";
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

std::pair<std::vector<ReadID>, std::vector<ReadID>> ReadClusteringEngine::get_spanning_tree_boundary_reads(SpanningTree &tree) {
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

    auto get_max_distance_pair = [](tsl::robin_map<ReadID, int> &distances) {
        return *std::max_element(distances.begin(), distances.end(), [](const std::pair<ReadID, int> &p1, const std::pair<ReadID, int> &p2) {
            return p1.second < p2.second;
        });
    };

    auto initial_distances = distance_bfs(adjacency_map.begin()->first);
    auto farthest_vertex_and_dist = get_max_distance_pair(initial_distances);

    std::vector<ReadID> left_side_vertices, right_side_vertices;

    auto distances_to_right = distance_bfs(farthest_vertex_and_dist.first);
    auto farthest_right_vertex_and_dist = get_max_distance_pair(distances_to_right);
    int tail_length = std::max(farthest_right_vertex_and_dist.second * 0.03, 1000.0);
    for (auto vertex_dist_pair : distances_to_right) {
        if (vertex_dist_pair.second + tail_length > farthest_right_vertex_and_dist.second) {
            right_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    auto distances_to_left = distance_bfs(farthest_right_vertex_and_dist.first);
    auto farthest_left_vertex_and_dist = get_max_distance_pair(distances_to_left);
    tail_length = std::max(farthest_left_vertex_and_dist.second * 0.03, 1000.0);
    for (auto vertex_dist_pair : distances_to_left) {
        if (vertex_dist_pair.second + tail_length > farthest_left_vertex_and_dist.second) {
            left_side_vertices.push_back(vertex_dist_pair.first);
        }
    }

    std::cout << fmt::format("Initial left {} - {}\n", read_intervals[farthest_left_vertex_and_dist.first].first, read_intervals[farthest_left_vertex_and_dist.first].second);
    std::cout << fmt::format("Initial right {} - {}\n", read_intervals[farthest_right_vertex_and_dist.first].first, read_intervals[farthest_right_vertex_and_dist.first].second);

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
    std::ofstream intervals;
    intervals.open(reader->meta.filename + "_intervals");

    std::multimap<ComponentID, std::vector<KmerID>> core_component_twin_tails;
    for (auto &component_and_tree : components_with_trees) {
        auto tails = get_spanning_tree_boundary_reads(component_and_tree.second);
        auto left_vertices = amplify_component(tails.first, 40);
        auto right_vertices = amplify_component(tails.second, 40);

        auto left_kmers = accumulate_kmer_ids(left_vertices);
        auto right_kmers = accumulate_kmer_ids(right_vertices);

        core_component_twin_tails.insert({component_and_tree.first[0], left_kmers});
        core_component_twin_tails.insert({component_and_tree.first[0], right_kmers});

        export_spanning_tree(component_and_tree.second, tails);

        auto left_tail_interval = get_read_coverage_intervals(tails.first);
        for (auto i : left_tail_interval) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << std::endl;
        auto right_tail_interval = get_read_coverage_intervals(tails.second);
        for (auto i : right_tail_interval) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << std::endl;
        auto component_interval = get_read_coverage_intervals(component_and_tree.first);
        for (auto i : component_interval) std::cout << fmt::format("{}-{} , ", i.first, i.second); std::cout << "\n\n";
        intervals << fmt::format("({}, {}), ({}, {}), ({}, {})\n", left_tail_interval[0].first, left_tail_interval[0].second, right_tail_interval[0].first, right_tail_interval[0].second,
                                 component_interval[0].first, component_interval[0].second);
    }
    intervals.close();

    std::vector<ComponentConnection> tail_connecting_edges;
    for (auto c_id_kmers_pair : core_component_twin_tails) {
        for (auto c_id_kmers_pair_2 : core_component_twin_tails) {
            if (c_id_kmers_pair.first < c_id_kmers_pair_2.first) {
                auto kmer_intersection = get_vectors_intersection(c_id_kmers_pair.second, c_id_kmers_pair_2.second);
                bool is_good = component_index[c_id_kmers_pair.first]->categories == component_index[c_id_kmers_pair_2.first]->categories;
                tail_connecting_edges.push_back({c_id_kmers_pair.first, c_id_kmers_pair_2.first, kmer_intersection.size(), is_good});
            }
        }
    }
    std::sort(tail_connecting_edges.rbegin(), tail_connecting_edges.rend());
    for (auto edge : tail_connecting_edges) {
        std::cout << fmt::format("{} : {} between {} and {}\n", edge.is_good ? '+' : '-', edge.score, component_index[edge.component_x_id]->to_string(),
                                 component_index[edge.component_y_id]->to_string());
        if (edge.score < 100) break;
    }

    ConnectionScore min_score;
    std::cin >> min_score;
    return filter_connections(tail_connecting_edges, [min_score](ComponentConnection &conn) { return conn.score >= min_score; });
}

std::vector<ComponentID> ReadClusteringEngine::run_clustering() {
    if (kmers == nullptr || k == 0) throw std::logic_error("You need to set the discriminative kmers first");

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

    auto get_core_component_ids = [this](int threshold_size) {
        std::vector<ComponentID> result;
        for (auto id_ptr_pair : component_index) {
            if (id_ptr_pair.second->size() >= threshold_size) {
                result.push_back(id_ptr_pair.first);
            }
        }
        return result;
    };

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", component_index.size(), reader->meta.records);

    int magic = 5;
    auto strong_initial_clusters = filter_components([magic](ReadComponent *c) { return c->discriminative_kmer_ids.size() > magic; });
    auto component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_connections, this, "Strong connections")(strong_initial_clusters, magic);
    plot_connection_quality(component_connections);

    auto threshold_value = std::find_if(component_connections.begin(), component_connections.end(), [](ComponentConnection &conn) { return !conn.is_good; })->score + 1;
    std::cout << fmt::format("Threshold value is {}\n", threshold_value);
    //plot_read_overlaps(component_connections);

    ConnectionScore min;
    std::cin >> min;
    auto strong_connections = filter_connections(component_connections, [min](ComponentConnection &conn) { return conn.score >= min; });
    //plot_overlap_differences(strong_connections);

    std::set<ComponentID> restricted;
    auto components_and_trees = union_find(strong_connections, restricted, 5);
    std::vector<Component> components = extract_components(components_and_trees);
    auto core_component_ids = timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Creation of strong clusters")(components);
    print_components(core_component_ids);


    auto core_component_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_core_component_connections, this, "Core connections")(components_and_trees);
    if (!core_component_connections.empty()) {
        components_and_trees = union_find(core_component_connections, restricted, 1);
        components = extract_components(components_and_trees);
        timeMeasureMemberFunc(&ReadClusteringEngine::merge_components, this, "Merging of core clusters")(components);
        remove_merged_components();
        core_component_ids = get_core_component_ids(30);
    } else {
        remove_merged_components();
    }

    print_components(core_component_ids);

    for (int i = 0; i < 3; i++) {
        component_connections = get_connections(core_component_ids, 1);
        plot_connection_quality(component_connections, core_component_ids);

        ConnectionScore min_score;
        std::cin >> min_score;
        component_connections = filter_connections(component_connections, [min_score](ComponentConnection &conn) { return conn.score >= min_score; });
        if (component_connections.empty()) break;

        restricted.clear();
        restricted.insert(core_component_ids.begin(), core_component_ids.end());
        components_and_trees = union_find(component_connections, restricted, 1);
        components = extract_components(components_and_trees);
        merge_components(components);
        remove_merged_components();

        core_component_ids = get_core_component_ids(30);

        print_components(core_component_ids);
    }

    return core_component_ids;
}

ReadClusteringEngine::~ReadClusteringEngine() {
    for (auto component_pair : component_index) {
        free(component_pair.second);
    }
}
