#include <numeric>
#include <random>
#include <iostream>
#include <climits>
#include "KernighanLin.h"

using namespace kln;

KernighanLin::KernighanLin(std::vector<Edge> &_edges){
    edges = std::vector<Edge>(_edges);

    for (auto edge : edges) {
        auto required_slots = std::max(edge.x, edge.y) + 1;
        if (required_slots > adjacencies.size()) {
            vertices.resize(required_slots, {PartitionA});
            adjacencies.resize(required_slots, {});
        }
        adjacencies[edge.x][edge.y] = edge.cost;
        adjacencies[edge.y][edge.x] = edge.cost;
    }

    std::vector<VertexID> indices(vertices.size());
    for (VertexID i = 0; i < vertices.size(); i++) indices[i] = i;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::shuffle(indices.begin(), indices.end(), rng);

    for (auto i = 0; i < indices.size() / 2; i++) {
        partition_a.insert(indices[i]);
        vertices[indices[i]].partition = PartitionA;
    }
    for (auto i = indices.size() / 2; i < vertices.size(); i++) {
        partition_b.insert(indices[i]);
        vertices[indices[i]].partition = PartitionB;
    }
}

KernighanLin::KernighanLin(std::vector<Edge> &_edges, std::pair<std::vector<VertexID>, std::vector<VertexID>> &partitions){
    edges = std::vector<Edge>(_edges);

    for (auto edge : edges) {
        auto required_slots = std::max(edge.x, edge.y) + 1;
        if (required_slots > adjacencies.size()) {
            vertices.resize(required_slots, {PartitionA});
            adjacencies.resize(required_slots, {});
        }
        adjacencies[edge.x][edge.y] = edge.cost;
        adjacencies[edge.y][edge.x] = edge.cost;
    }

    for (auto id : partitions.first){
        partition_a.insert(id);
        vertices[id].partition = PartitionA;
    }
    for (auto id : partitions.second){
        partition_b.insert(id);
        vertices[id].partition = PartitionB;
    }
}

Gain KernighanLin::get_gain(VertexID vertex_a_id, VertexID vertex_b_id){
    auto weight_it = adjacencies[vertex_a_id].find(vertex_b_id);
    Cost w = weight_it != adjacencies[vertex_a_id].end() ? (*weight_it).second : 0;
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

void KernighanLin::compute_deltas() {
    deltas.clear();
    deltas.resize(vertices.size());

    for (auto edge : edges) {
        if (vertices[edge.x].partition != vertices[edge.y].partition){
            deltas[edge.x] += edge.cost;
            deltas[edge.y] += edge.cost;
        } else {
            deltas[edge.x] -= edge.cost;
            deltas[edge.y] -= edge.cost;
        }
    }
}

void KernighanLin::print_deltas(){
    for (VertexID v_id = 0; v_id < vertices.size(); v_id++){
        std::cout << "Vertex " << v_id << " delta " << deltas[v_id] << std::endl;
    }
}

void KernighanLin::update_deltas(VertexID swapped_from_a, VertexID swapped_from_b) {
    for (auto it = begin(adjacencies[swapped_from_a]); it != end(adjacencies[swapped_from_a]); ++it) {
        if (vertices[it->first].partition != vertices[swapped_from_a].partition) {
            deltas[it->first] -= 2 * it->second;
        } else {
            deltas[it->first] += 2 * it->second;
        }

        std::cout << "Updated " << it->first << " to " << deltas[it->first] << std::endl;
    }
    for (auto it = begin(adjacencies[swapped_from_b]); it != end(adjacencies[swapped_from_b]); ++it) {
        if (vertices[it->first].partition != vertices[swapped_from_b].partition) {
            deltas[it->first] -= 2 * it->second;
        } else {
            deltas[it->first] += 2 * it->second;
        }

        std::cout << "Updated " << it->first << " to " << deltas[it->first] << std::endl;
    }
}

bool KernighanLin::bisection_pass(){
    std::vector<SwapGain> best_gains;
    std::set<VertexID> swapped;

    compute_deltas();

    while (swapped.size() < vertices.size()){
        auto best_gain = get_best_gain(swapped);
        std::cout << "Best gain = " << best_gain.gain << " " <<  best_gain.vertex_a << " " << best_gain.vertex_b << std::endl;

        if (best_gain.gain == INT_MIN) break;

        swapped.insert(best_gain.vertex_a);
        swapped.insert(best_gain.vertex_b);

        update_deltas(best_gain.vertex_a, best_gain.vertex_b);

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

        vertices[swapped_from_a].partition = PartitionB;
        vertices[swapped_from_b].partition = PartitionA;
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
        if (vertices[edge.x].partition != vertices[edge.y].partition){
            cost += edge.cost;
        }
    }
    return cost;
}

void KernighanLin::bisection(int max_iter){
    for (int iter_count = 0; iter_count < max_iter; iter_count++){
        std::cout << "Iteration " << iter_count + 1 << std::endl;
        std::cout << "Cut cost is " << cut_cost() << std::endl;
        if (!bisection_pass()) break;
    }
}