#include "KernighanLinClustering.h"
#include <boost/algorithm/string/join.hpp>

#include "../common/Utils.h"
#include "../lib/QuickCut.h"

namespace algo = boost::algorithm;

using namespace quick_cut;

std::string KernighanLinClustering::partition_consistency(std::set<VertexID> &vertex_ids) {
    std::map<CategoryID, int> category_counts;
    for (auto vertex_id : vertex_ids) {
        category_counts.insert(std::map<CategoryID, int>::value_type(cluster_category_map[vertex_id], 0)).first->second += 1;
    }
    std::vector<std::string> category_count_vector;
    std::transform(
            category_counts.begin(),
            category_counts.end(),
            std::back_inserter(category_count_vector),
            [](std::pair<const CategoryID, int> &p) -> std::string { return fmt::format("{}", p.second); }
    );
    return algo::join(category_count_vector, "/");
}


void KernighanLinClustering::construct_cluster_category_map() {
    for (auto it = begin(cluster_index); it != end(cluster_index); it++){
        cluster_category_map.insert({it->first, *(it->second->categories.begin())});
    }
}



void KernighanLinClustering::run_clustering() {
    std::vector<ClusterConnection> cluster_connections = timeMeasureMemberFunc(&KernighanLinClustering::get_connections, this, "Cluster connections")();
    //plot_connection_quality(cluster_connections);
    construct_cluster_category_map();

    std::vector<Edge> edges;

    std::ofstream output;
    output.open("graph_pls");

    for (auto conn : cluster_connections){
        Cost w = conn.score * conn.score;
        edges.push_back({conn.cluster_x_id, conn.cluster_y_id, w});

        output << fmt::format("{} {} {}\n", conn.cluster_x_id, conn.cluster_y_id, w);
    }
    output.close();

    std::cout << fmt::format("{} edges\n", edges.size());

    std::vector<VertexID> a, b;
    for (auto it = begin(cluster_index); it != end(cluster_index); it++){
        if (*(it->second->categories.begin()) == 0) a.push_back(it->second->reference_id); else b.push_back(it->second->reference_id);
    }
    std::cout << a.size() << " " << b.size() << std::endl;
    auto parts = std::make_pair(a, b);
    auto p = QuickCut(edges);
    std::cout << "Initial cost " << p.cut_cost() << std::endl;
    for (int iter = 0; iter < 100; iter++){
        p.bisection_pass();

        std::cout << partition_consistency(p.partition_a) << std::endl;
        std::cout << partition_consistency(p.partition_b) << std::endl;
        std::cout << "Cost after pass " << p.cut_cost() << std::endl << std::endl;
    }
    std::cout << "Final cost" << p.cut_cost() << std::endl;
}
