#include <numeric>
#include <random>
#include <iostream>
#include "QuickCut.h"

using namespace quick_cut;


Cost QuickCut::get_cost(VertexID a, VertexID b) {
    auto weight_it = adjacencies[a].find(b);
    return weight_it != adjacencies[a].end() ? weight_it->second : 0;
}

Gain QuickCut::get_gain(std::reverse_iterator<PartitionSet::iterator> &a_it, std::reverse_iterator<PartitionSet::iterator> &b_it) {
    return a_it->first + b_it->first - 2 * get_cost(a_it->second, b_it->second);
}


bool QuickCut::is_neighbor(VertexID a, VertexID b) {
    return get_cost(a, b) != 0;
}

Cost QuickCut::cut_cost() {
    Cost cost = 0;
    for (auto edge : edges) {
        if (vertices[edge.x].partition != vertices[edge.y].partition) {
            cost += edge.c;
        }
    }
    return cost;
}


void QuickCut::compute_deltas() {
    std::vector<Delta> deltas(vertices.size(), 0);
    for (auto edge : edges) {
        if (vertices[edge.x].partition != vertices[edge.y].partition) {
            deltas[edge.x] += edge.c;
            deltas[edge.y] += edge.c;
        } else {
            deltas[edge.x] -= edge.c;
            deltas[edge.y] -= edge.c;
        }
    }

    partition_a_set.clear();
    partition_b_set.clear();

    for (VertexID vertex_id = 0; vertex_id < vertices.size(); vertex_id++) {
        vertices[vertex_id].delta = deltas[vertex_id];
        if (vertices[vertex_id].partition == PartitionA) {
            partition_a_set.insert({deltas[vertex_id], vertex_id});
        } else {
            partition_b_set.insert({deltas[vertex_id], vertex_id});
        }
    }
}


void QuickCut::update_deltas(VertexID a, VertexID b) {
    for (auto it = begin(adjacencies[a]); it != end(adjacencies[a]); ++it) {
        if (vertices[it->first].swap_candidate) continue;

        if (vertices[it->first].partition != vertices[a].partition) {
            partition_b_set.erase({vertices[it->first].delta, it->first});
            vertices[it->first].delta -= 2 * it->second;
            partition_b_set.insert({vertices[it->first].delta, it->first});
        } else {
            partition_a_set.erase({vertices[it->first].delta, it->first});
            vertices[it->first].delta += 2 * it->second;
            partition_a_set.insert({vertices[it->first].delta, it->first});
        }

        //std::cout << "Updated " << it->first << " to " << vertices[it->first].delta << std::endl;
    }
    for (auto it = begin(adjacencies[b]); it != end(adjacencies[b]); ++it) {
        if (vertices[it->first].swap_candidate) continue;

        if (vertices[it->first].partition != vertices[b].partition) {
            partition_a_set.erase({vertices[it->first].delta, it->first});
            vertices[it->first].delta -= 2 * it->second;
            partition_a_set.insert({vertices[it->first].delta, it->first});
        } else {
            partition_b_set.erase({vertices[it->first].delta, it->first});
            vertices[it->first].delta += 2 * it->second;
            partition_b_set.insert({vertices[it->first].delta, it->first});
        }

        //std::cout << "Updated " << it->first << " to " << vertices[it->first].delta << std::endl;
    }
}


SwapGain QuickCut::get_best_gain() {
    auto a_iterator_first_k = partition_a_set.rbegin();
    auto b_iterator = partition_b_set.rbegin();
    while (a_iterator_first_k != partition_a_set.rend() && is_neighbor(a_iterator_first_k->second, b_iterator->second)) {
        ++a_iterator_first_k;
    }

    SwapGain max_gain = {get_gain(a_iterator_first_k, b_iterator), a_iterator_first_k->second, b_iterator->second};

    for (auto a_iterator = partition_a_set.rbegin(); a_iterator != a_iterator_first_k; ++a_iterator) {
        b_iterator = partition_b_set.rbegin();
        while (b_iterator != partition_b_set.rend()) {
            Gain g = get_gain(a_iterator, b_iterator);
            if (g > max_gain.gain) {
                max_gain.gain = g;
                max_gain.a_vertex = a_iterator->second;
                max_gain.b_vertex = b_iterator->second;
            }

            if (b_iterator == partition_b_set.rbegin() && a_iterator_first_k->first + b_iterator->first <= max_gain.gain) {
                return max_gain;
            }

            if (!is_neighbor(a_iterator->second, b_iterator->second)) {
                break;
            }
        }
    }

    return max_gain;
}


void QuickCut::print_deltas(){
    for (VertexID v_id = 0; v_id < vertices.size(); v_id++){
        if (vertices[v_id].partition == PartitionA){
            if (!vertices[v_id].swap_candidate && partition_a_set.find({vertices[v_id].delta, v_id}) == partition_a_set.end()){
                throw std::logic_error("WTF");
            }
        } else {
            if (!vertices[v_id].swap_candidate && partition_b_set.find({vertices[v_id].delta, v_id}) == partition_b_set.end()){
                throw std::logic_error("WTF");
            }
        }
        std::cout << "Vertex " << v_id << " delta " << vertices[v_id].delta << std::endl;
    }
}


bool QuickCut::bisection_pass() {
    std::vector<SwapGain> best_gains;

    compute_deltas();

    while (best_gains.size() < vertices.size() / 2) {
        auto best_gain = get_best_gain();
        //std::cout << "Best gain = " << best_gain.gain << " " <<  best_gain.a_vertex << " " << best_gain.b_vertex << std::endl;

        best_gains.push_back(best_gain);

        update_deltas(best_gain.a_vertex, best_gain.b_vertex);

        partition_a_set.erase({vertices[best_gain.a_vertex].delta, best_gain.a_vertex});
        partition_b_set.erase({vertices[best_gain.b_vertex].delta, best_gain.b_vertex});

        vertices[best_gain.a_vertex].swap_candidate = true;
        vertices[best_gain.b_vertex].swap_candidate = true;
    }

    std::vector<SwapGain> partial_gain_sums(best_gains.size());
    std::partial_sum(best_gains.begin(), best_gains.end(), partial_gain_sums.begin(), [](const SwapGain &accumulator, const SwapGain &current) -> SwapGain {
        return {accumulator.gain + current.gain, current.a_vertex, current.b_vertex};
    });

    auto max_gain_it = std::max_element(partial_gain_sums.begin(), partial_gain_sums.end(), [](const SwapGain &g1, const SwapGain &g2) -> bool {
        return g1.gain < g2.gain;
    });
    if ((*max_gain_it).gain <= 0) return false;

    for (auto it = begin(partial_gain_sums); it != std::next(max_gain_it); it++) {
        vertices[it->a_vertex].partition = PartitionB;
        vertices[it->b_vertex].partition = PartitionA;

        partition_a.erase(it->a_vertex);
        partition_b.erase(it->b_vertex);

        partition_a.insert(it->b_vertex);
        partition_b.insert(it->a_vertex);
    }

    for (auto gain : best_gains){
        vertices[gain.a_vertex].swap_candidate = false;
        vertices[gain.b_vertex].swap_candidate = false;
    }
    return true;
}

void QuickCut::bisection(int max_iter) {
    for (int iter_count = 0; iter_count < max_iter; iter_count++) {
        std::cout << "Iteration " << iter_count + 1 << std::endl;
        std::cout << "Cut cost is " << cut_cost() << std::endl;
        if (!bisection_pass()) break;
    }
}

QuickCut::QuickCut(std::vector<Edge> &_edges) {
    edges = std::vector<Edge>(_edges);

    for (auto edge : edges) {
        auto required_slots = std::max(edge.x, edge.y) + 1;
        if (required_slots > vertices.size()) {
            vertices.resize(required_slots, {0, PartitionA, false});
            adjacencies.resize(required_slots, {});
        }
        adjacencies[edge.x][edge.y] = edge.c;
        adjacencies[edge.y][edge.x] = edge.c;
    }

    std::vector<VertexID> indices(vertices.size());
    for (VertexID i = 0; i < vertices.size(); i++) indices[i] = i;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::shuffle(indices.begin(), indices.end(), rng);

    for (auto i = 0; i < indices.size() / 2; i++) {
        vertices[indices[i]].partition = PartitionA;
        partition_a.insert(indices[i]);
    }
    for (auto i = indices.size() / 2; i < vertices.size(); i++) {
        vertices[indices[i]].partition = PartitionB;
        partition_b.insert(indices[i]);
    }
}

QuickCut::QuickCut(std::vector<Edge> &_edges, std::pair<std::vector<VertexID>, std::vector<VertexID> > &partitions) {
    edges = std::vector<Edge>(_edges);

    for (auto edge : edges) {
        auto required_slots = std::max(edge.x, edge.y) + 1;
        if (required_slots > vertices.size()) {
            vertices.resize(required_slots, {0, PartitionA, false});
            adjacencies.resize(required_slots, {});
        }
        adjacencies[edge.x][edge.y] = edge.c;
        adjacencies[edge.y][edge.x] = edge.c;
    }
    for (VertexID id : partitions.first){
        vertices[id].partition = PartitionA;
        partition_a.insert(id);
    }
    for (VertexID id : partitions.second){
        vertices[id].partition = PartitionB;
        partition_b.insert(id);
    }
}
