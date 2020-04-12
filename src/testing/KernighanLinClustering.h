#ifndef SRC_KERNIGHANLINCLUSTERING_H
#define SRC_KERNIGHANLINCLUSTERING_H


#include "../common/BaseReadClusteringEngine.h"
#include "../lib/KernighanLin.h"

class KernighanLinClustering : public BaseReadClusteringEngine {
    tsl::robin_map<ClusterID, CategoryID> cluster_category_map;
    void construct_cluster_category_map();
    std::string partition_consistency(std::set<VertexID> &vertex_ids);
public:
    void run_clustering() override;
    using BaseReadClusteringEngine::BaseReadClusteringEngine;
};


#endif //SRC_KERNIGHANLINCLUSTERING_H
