#ifndef SRC_QUICKCUT_H
#define SRC_QUICKCUT_H

#include <set>
#include <vector>
#include <map>

namespace quick_cut {

    typedef long long Delta;
    typedef uint32_t VertexID;
    typedef long long Gain;
    typedef std::set<std::pair<Delta, VertexID> > PartitionSet;
    typedef long long Cost;

    enum Partition {
        PartitionA, PartitionB
    };

    struct Vertex {
        Delta delta;
        Partition partition;
        bool swap_candidate;
    };

    struct Edge {
        VertexID x;
        VertexID y;
        Cost c;
    };

    struct SwapGain {
        Gain gain;
        VertexID a_vertex;
        VertexID b_vertex;
    };

    class QuickCut {
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;

        std::vector<std::map<VertexID, Cost>> adjacencies;

        PartitionSet partition_a_set, partition_b_set;
        SwapGain get_best_gain();

        Cost get_cost(VertexID a, VertexID b);

        Gain get_gain(std::reverse_iterator<PartitionSet::iterator> &a_it, std::reverse_iterator<PartitionSet::iterator> &b_it);

        bool is_neighbor(VertexID a, VertexID b);

        void compute_deltas();

        void update_deltas(VertexID a, VertexID b);

    public:
        std::set<VertexID> partition_a, partition_b;

        explicit QuickCut(std::vector<Edge> &_edges);
        QuickCut(std::vector<Edge> &_edges, std::pair<std::vector<VertexID>, std::vector<VertexID>> &partitions);

        bool bisection_pass();
        void bisection(int max_iter);

        Cost cut_cost();

        void print_deltas();
    };

}

#endif //SRC_QUICKCUT_H
