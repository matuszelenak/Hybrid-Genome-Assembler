#ifndef SRC_KERNIGHANLIN_H
#define SRC_KERNIGHANLIN_H

#include <map>
#include <set>
#include <vector>
#include <tsl/robin_map.h>

typedef uint32_t VertexID;
typedef uint64_t Weight;
typedef int Delta;
typedef int Gain;
typedef tsl::robin_map<VertexID, Weight> NeighborMap;

enum Category {
    CategoryA, CategoryB
};

struct Vertex {
    VertexID id;
    Category category;
    NeighborMap neighbors;
};

struct Edge {
    Vertex *vertex_a;
    Vertex *vertex_b;
    Weight weight;
};

struct SwapGain {
    Gain gain;
    VertexID vertex_a;
    VertexID vertex_b;
};

class KernighanLin {
private:
    tsl::robin_map<VertexID, Delta> deltas;
    tsl::robin_map<VertexID, Vertex *> vertices;
    std::vector<Edge *> edges;

    Gain get_gain(VertexID vertex_a_id, VertexID vertex_b_id);
    SwapGain get_best_gain(std::set<VertexID> &restricted);

    tsl::robin_map<VertexID, Delta> compute_deltas();
    void update_deltas(tsl::robin_map<VertexID, Delta> &_deltas, Vertex *swapped_from_a, Vertex *swapped_from_b);

public:
    KernighanLin(tsl::robin_map<VertexID, Vertex *> &_vertices, std::vector<Edge *> &_edges);
    KernighanLin(tsl::robin_map<VertexID, Vertex*> &_vertices, std::vector<Edge*> &_edges, std::pair<std::vector<VertexID>, std::vector<VertexID>> &partitions);

    std::set<VertexID> partition_a;
    std::set<VertexID> partition_b;

    bool bisection_pass();
    void bisection(int max_iter);

    void print_partitions();
    int cut_cost();
};


#endif //SRC_KERNIGHANLIN_H
