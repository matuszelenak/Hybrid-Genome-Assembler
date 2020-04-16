#ifndef SRC_KERNIGHANLIN_H
#define SRC_KERNIGHANLIN_H

#include <map>
#include <set>
#include <vector>

namespace kln {

    typedef uint32_t VertexID;
    typedef uint64_t Cost;
    typedef int Gain;
    typedef int Delta;

    enum Partition {
        PartitionA, PartitionB
    };

    struct Edge {
        VertexID x;
        VertexID y;
        Cost cost;
    };

    struct Vertex {
        Partition partition;
    };

    struct SwapGain {
        Gain gain;
        VertexID vertex_a;
        VertexID vertex_b;
    };

    class KernighanLin {
    private:
        std::vector<Delta> deltas;
        std::vector<Edge> edges;
        std::vector<Vertex> vertices;

        std::vector<std::map<VertexID, Cost>> adjacencies;

        Gain get_gain(VertexID vertex_a_id, VertexID vertex_b_id);

        SwapGain get_best_gain(std::set<VertexID> &restricted);

        void compute_deltas();

        void update_deltas(VertexID swapped_from_a, VertexID swapped_from_b);

    public:
        explicit KernighanLin(std::vector<Edge> &edges);

        KernighanLin(std::vector<Edge> &edges, std::pair<std::vector<VertexID>, std::vector<VertexID>> &partitions);

        std::set<VertexID> partition_a;
        std::set<VertexID> partition_b;

        bool bisection_pass();

        void bisection(int max_iter);

        void print_partitions();

        int cut_cost();

        void print_deltas();
    };
}

#endif //SRC_KERNIGHANLIN_H
