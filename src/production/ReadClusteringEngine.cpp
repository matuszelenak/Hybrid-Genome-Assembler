#include <iostream>
#include "ReadClusteringEngine.h"
#include "../common/KmerIterator.h"
#include "../common/Utils.h"


void ReadClusteringEngine::run_clustering(){
    while (clustering_round() > 0){}
    std::cout << fmt::format("Exported {} clusters\n", timeMeasureMemberFunc(&ReadClusteringEngine::export_clusters, this, "Exporting clusters")(10));
}