#include <numeric>
#include <random>
#include <iostream>
#include <climits>
#include "KernighanLin.h"
#include "../common/Utils.h"

KernighanLin::KernighanLin(tsl::robin_map<VertexID, Vertex*> &_vertices, std::vector<Edge*> &_edges){
    edges = _edges;
    vertices = _vertices;

    auto rng = std::default_random_engine {0};
    std::vector<VertexID>vertex_indices;
    for (auto it = begin(vertices); it != end(vertices); it++){
        vertex_indices.push_back(it->first);
    }
    std::shuffle(vertex_indices.begin(), vertex_indices.end(), rng);

    for (auto i = 0; i < vertex_indices.size() / 2; i++){
        partition_a.insert(vertex_indices[i]);
        vertices[vertex_indices[i]]->category = CategoryA;
    }
    for (auto i = vertex_indices.size() / 2; i < vertices.size(); i++){
        partition_b.insert(vertex_indices[i]);
        vertices[vertex_indices[i]]->category = CategoryB;
    }
}

KernighanLin::KernighanLin(tsl::robin_map<VertexID, Vertex*> &_vertices, std::vector<Edge*> &_edges, std::pair<std::vector<VertexID>, std::vector<VertexID>> &partitions){
    edges = _edges;
    vertices = _vertices;

    for (auto id : partitions.first){
        partition_a.insert(id);
        vertices[id]->category = CategoryA;
    }
    for (auto id : partitions.second){
        partition_b.insert(id);
        vertices[id]->category = CategoryB;
    }
}

Gain KernighanLin::get_gain(VertexID vertex_a_id, VertexID vertex_b_id){
    auto weight_it = vertices[vertex_a_id]->neighbors.find(vertex_b_id);
    Weight w = weight_it != vertices[vertex_a_id]->neighbors.end() ? (*weight_it).second : 0;
    return deltas[vertex_a_id] + deltas[vertex_b_id] - 2 * w;
}


SwapGain KernighanLin::get_best_gain(std::set<VertexID> &restricted) {
    std::vector<VertexID> allowed_a, allowed_b;
    std::set_difference(partition_a.begin(), partition_a.end(), restricted.begin(), restricted.end(), std::back_inserter(allowed_a));
    std::set_difference(partition_b.begin(), partition_b.end(), restricted.begin(), restricted.end(), std::back_inserter(allowed_b));

    SwapGain best_gain = {INT_MIN, 0, 0};
    for (auto vertex_a_id : allowed_a) {
        for (auto vertex_b_id: allowed_b) {
            Gain gain = get_gain(vertex_a_id, vertex_b_id);
            if (gain > best_gain.gain){
                best_gain = {gain, vertex_a_id, vertex_b_id};
            }
        }
    }
    return best_gain;
}

tsl::robin_map<VertexID, Delta> KernighanLin::compute_deltas() {
    tsl::robin_map<VertexID, Delta> _deltas;

    auto num_of_edges = edges.size();
    uint32_t edge_cntr = 1;
    for (auto edge : edges) {
        _deltas.insert({edge->vertex_a->id, 0});
        _deltas.insert({edge->vertex_b->id, 0});

        if (edge->vertex_a->category != edge->vertex_b->category) {
            _deltas[edge->vertex_a->id] += edge->weight;
            _deltas[edge->vertex_b->id] += edge->weight;
        } else {
            _deltas[edge->vertex_a->id] -= edge->weight;
            _deltas[edge->vertex_b->id] -= edge->weight;
        }

        show_progress(edge_cntr, num_of_edges, "Calculating deltas");
        edge_cntr++;
    }
    return _deltas;
}

void KernighanLin::update_deltas(tsl::robin_map<VertexID, Delta> &_deltas, Vertex *swapped_from_a, Vertex *swapped_from_b) {
    for (auto it = begin(swapped_from_a->neighbors); it != end(swapped_from_a->neighbors); ++it) {
        if (vertices[it->first]->category != swapped_from_a->category) {
            _deltas[it->first] -= 2 * it->second;
        } else {
            _deltas[it->first] += 2 * it->second;
        }
    }
    for (auto it = begin(swapped_from_b->neighbors); it != end(swapped_from_b->neighbors); ++it) {
        if (vertices[it->first]->category != swapped_from_b->category) {
            _deltas[it->first] -= 2 * it->second;
        } else {
            _deltas[it->first] += 2 * it->second;
        }
    }
}

bool KernighanLin::bisection_pass(){
    std::vector<SwapGain> best_gains;
    std::set<VertexID> swapped;

    deltas = compute_deltas();
    while (swapped.size() < vertices.size()){
        auto best_gain = get_best_gain(swapped);
        if (best_gain.gain == INT_MIN) break;

        swapped.insert(best_gain.vertex_a);
        swapped.insert(best_gain.vertex_b);

        update_deltas(deltas, vertices[best_gain.vertex_a], vertices[best_gain.vertex_b]);

        best_gains.push_back(best_gain);
    }

    std::vector<SwapGain>partial_gain_sums(best_gains.size());
    std::partial_sum(best_gains.begin(), best_gains.end(), partial_gain_sums.begin(), [](const SwapGain &accumulator, const SwapGain &current)-> SwapGain {
        return {accumulator.gain + current.gain, current.vertex_a, current.vertex_b};
    });

    auto max_gain_it = std::max_element(partial_gain_sums.begin(), partial_gain_sums.end(), [](const SwapGain &g1, const SwapGain &g2)-> bool {
        return g1.gain < g2.gain;
    });
    if ((*max_gain_it).gain <= 0) return false;

    for (auto it = begin(partial_gain_sums); it != std::next(max_gain_it); it++){
        VertexID swapped_from_a = it->vertex_a;
        VertexID swapped_from_b = it->vertex_b;

        partition_a.erase(swapped_from_a);
        partition_b.erase(swapped_from_b);

        partition_a.insert(swapped_from_b);
        partition_b.insert(swapped_from_a);

        vertices[swapped_from_a]->category = CategoryB;
        vertices[swapped_from_b]->category = CategoryA;
    }
    return true;
}

void KernighanLin::print_partitions(){
    for (auto v : partition_a){
        std::cout << v << " ";
    }
    std::cout << std::endl;
    for (auto v : partition_b){
        std::cout << v << " ";
    }
    std::cout << std::endl;
}

int KernighanLin::cut_cost(){
    int cost = 0;
    for (auto edge : edges){
        if (edge->vertex_a->category != edge->vertex_b->category){
            cost += edge->weight;
        }
    }
    return cost;
}

void KernighanLin::bisection(int max_iter){
    for (int iter_count = 0; iter_count < max_iter; iter_count++){
        std::cout << "Iteration " << iter_count + 1 << std::endl;
        std::cout << fmt::format("Cut cost is {}\n", cut_cost());
        if (!bisection_pass()) break;
    }
}