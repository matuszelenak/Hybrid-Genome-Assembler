#include <iostream>
#include <fmt/format.h>
#include <queue>
#include <experimental/filesystem>
#include <boost/algorithm/string/join.hpp>

#include "KmerIterator.h"
#include "ReadClusteringEngine.h"

namespace fs = std::experimental::filesystem;
namespace algo = boost::algorithm;

std::mutex connections_mut;
std::mutex component_erase_mut;
std::mutex index_merge;


void plot_connection_quality(std::vector<ClusterConnection> &connections) {
    std::map<bool, std::map<ConnectionScore, int>> score_to_matching;
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


void plot_cluster_coverage( std::vector<GenomeReadCluster *> &clusters) {
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

std::vector<ClusterConnection> filter_connections(std::vector<ClusterConnection> &original, const std::function<bool(ClusterConnection&)>& func){
    std::vector<ClusterConnection> result;
    for (auto conn : original){
        if (func(conn)) result.push_back(conn);
    }
    return result;
}

ClusterIDPair get_cluster_id_pair(ClusterID x, ClusterID y){
    if (x < y){
        return (((ClusterIDPair)x) << 32u) | y;
    } else {
        return (((ClusterIDPair)y) << 32u) | x;
    }
}

void ReadClusteringEngine::print_clusters(int first_n) {
    std::vector<GenomeReadCluster *> cluster_pointers;
    std::transform(
            cluster_index.begin(),
            cluster_index.end(),
            std::back_inserter(cluster_pointers),
            [](std::pair<const ClusterID, GenomeReadCluster *> &p) -> GenomeReadCluster * { return p.second; }
    );

    sort(cluster_pointers.rbegin(), cluster_pointers.rend(), [](GenomeReadCluster *x, GenomeReadCluster *y) -> bool { return x->size() < y->size(); });

    int iterate_first = first_n;
    if (first_n == -1) iterate_first = cluster_pointers.size();

    int size_one_clusters = 0;
    for (auto cluster : cluster_pointers) {
        if (cluster->size() > 1) {
            std::cout << cluster->to_string() << std::endl;
        } else {
            size_one_clusters++;
        }

        iterate_first--;
        if (iterate_first == 0) break;
    }
    std::cout << fmt::format("As well as {} clusters of size 1\n", size_one_clusters);
}

ReadClusteringEngine::ReadClusteringEngine(SequenceRecordIterator &read_iterator, int k, bloom::BloomFilter<Kmer> &kmers, Platform platform) {
    this->platform = platform;
    this->k = k;
    this->reader = &read_iterator;
    this->kmers = &kmers;

    timeMeasureMemberFunc(&ReadClusteringEngine::construct_indices, this, "Construct indices")();
    std::cout << fmt::format("{}/{} reads converted to clusters\n", cluster_index.size(), reader->meta.records);
}

void ReadClusteringEngine::construct_indices_thread(KmerIndex &kmer_index) {
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        KmerIterator it = KmerIterator(read->sequence, k);

        tsl::robin_set<Kmer> in_read_discriminative;
        while (it.next_kmer()) {
            if (kmers->contains(it.current_kmer)) {
                in_read_discriminative.insert(it.current_kmer);
            }
        }

        if (in_read_discriminative.size() >= 5) {
            std::vector<KmerID> in_read_discriminative_ids;
            ClusterID cluster_id = read->id;

            index_merge.lock();

            std::pair<KmerIndex::iterator, bool> insert_result;
            KmerID new_kmer_id = kmer_index.size();
            for (Kmer kmer : in_read_discriminative) {
                insert_result = kmer_index.insert(KmerIndex::value_type(kmer, new_kmer_id));
                if (insert_result.second) {
                    kmer_cluster_index.push_back({});
                    ++new_kmer_id;
                }
                in_read_discriminative_ids.push_back(insert_result.first->second);
                kmer_cluster_index[insert_result.first->second].push_back(cluster_id);
            }

            std::sort(in_read_discriminative_ids.begin(), in_read_discriminative_ids.end());
            ReadMetaData meta = {read->id, read->category_id, read->start, read->end};
            cluster_index.insert(ClusterIndex::value_type(cluster_id, new GenomeReadCluster(meta, in_read_discriminative_ids)));

            index_merge.unlock();
        }
    }
}


int ReadClusteringEngine::construct_indices() {
    reader->rewind();

    KmerIndex kmer_index;
    run_in_threads(&ReadClusteringEngine::construct_indices_thread, this, std::ref(kmer_index));

    kmer_id_to_kmer.resize(kmer_index.size());
    for (auto kmer_id_pair : kmer_index){
        kmer_id_to_kmer[kmer_id_pair.second] = kmer_id_pair.first;
    }

    for (auto& arr : kmer_cluster_index){
        std::sort(arr.begin(), arr.end());
    }
    std::cout << fmt::format("Kmer cluster index size {}\n", kmer_cluster_index.size());
    return 0;
}


void ReadClusteringEngine::get_connections_thread(ConnectionScore min_score, ConcurrentQueue<ClusterID> &cluster_id_queue, std::vector<ClusterConnection> &accumulator) {
    while (true) {
        auto pivot = cluster_id_queue.pop();
        if (pivot == std::nullopt) break;
        ClusterID pivot_id = *pivot;

        std::map<ClusterID, ConnectionScore> shared_kmer_counts;
        for (KmerID kmer_id : cluster_index[pivot_id]->discriminative_kmer_ids) {
            auto stop = std::lower_bound(kmer_cluster_index[kmer_id].begin(), kmer_cluster_index[kmer_id].end(), pivot_id);
            for (auto candidate_it = kmer_cluster_index[kmer_id].begin(); candidate_it != stop; ++candidate_it) {
                shared_kmer_counts.insert(std::map<ClusterID, ConnectionScore>::value_type(*candidate_it, 0)).first->second += 1;
            }
        }
        shared_kmer_counts.erase(pivot_id);

        for (auto iter = begin(shared_kmer_counts); iter != end(shared_kmer_counts); iter++) {
            if (iter->second < min_score) continue;

            connections_mut.lock();
            accumulator.push_back({iter->first, pivot_id, iter->second, cluster_index[iter->first]->categories == cluster_index[pivot_id]->categories});
            connections_mut.unlock();
        }
    }
}

std::vector<ClusterConnection> ReadClusteringEngine::get_all_connections(ConnectionScore min_score){
    std::vector<ClusterID> cluster_ids;
    for (auto c_pair : cluster_index) cluster_ids.push_back(c_pair.first);
    return get_connections(cluster_ids, min_score);
}


std::vector<ClusterConnection> ReadClusteringEngine::get_connections(std::vector<ClusterID> &cluster_ids, ConnectionScore min_score) {
    std::vector<ClusterConnection> connections;

    auto queue = ConcurrentQueue<ClusterID>(cluster_ids.begin(), cluster_ids.end());
    run_in_threads(&ReadClusteringEngine::get_connections_thread, this, min_score, std::ref(queue), std::ref(connections));

    std::sort(connections.rbegin(), connections.rend());
    return connections;
}

void ReadClusteringEngine::merge_clusters_thread(ConcurrentQueue<IDComponent> &component_queue, IndexRemovalMap &for_removal, std::vector<ClusterID> &merge_result_ids) {
    while (true) {
        auto component_opt = component_queue.pop();
        if (component_opt == std::nullopt) break;
        auto component = *component_opt;

        GenomeReadCluster *survivor = cluster_index[component[0]];

        std::vector<std::vector<KmerID> *> kmer_merge_queue;
        std::set<CategoryID> survivor_categories;
        std::vector<ReadMetaData> survivor_contained_reads;

        for (auto cluster_id : component){
            auto cluster = cluster_index[cluster_id];
            kmer_merge_queue.push_back(&cluster->discriminative_kmer_ids);
            survivor_contained_reads.insert(survivor_contained_reads.end(), cluster->contained_reads.begin(), cluster->contained_reads.end());
            survivor_categories.insert(cluster->categories.begin(), cluster->categories.end());
        }

        survivor->discriminative_kmer_ids = merge_n_vectors(kmer_merge_queue, true);
        survivor->categories = survivor_categories;
        survivor->contained_reads = survivor_contained_reads;

        component_erase_mut.lock();
        merge_result_ids.push_back(survivor->id);
        for (auto cluster_id : component){
            for (KmerID kmer_id : cluster_index[cluster_id]->discriminative_kmer_ids) {
                for_removal.insert({kmer_id, {}}).first.value().push_back(cluster_id);
            }
        }
        for (int i = 1; i < component.size(); i++){
            free(cluster_index[component[i]]);
            cluster_index.erase(component[i]);
        }
        component_erase_mut.unlock();
    }
}

void ReadClusteringEngine::kmer_cluster_index_update(ConcurrentQueue<IndexRemovalMap::value_type> &removal_list_queue){
    while (true){
        auto removal_pair_opt = removal_list_queue.pop();
        if (removal_pair_opt == std::nullopt) break;

        std::vector<ClusterID>* removal_list = &removal_pair_opt->second;
        KmerID kmer_id = removal_pair_opt->first;
        std::sort(removal_list->begin(), removal_list->end());

        std::vector<ClusterID> updated_list;
        int i = 0, j = 0;
        while (i < removal_list->size() && j < kmer_cluster_index[kmer_id].size()) {
            if ((*removal_list)[i] < kmer_cluster_index[kmer_id][j]){
                i++;
            }
            else if (kmer_cluster_index[kmer_id][j] < (*removal_list)[i]){
                updated_list.push_back(kmer_cluster_index[kmer_id][j]);
                j++;
            }
            else {  // Should be removed
                i++;
                j++;
            }
        }
        kmer_cluster_index[kmer_id] = updated_list;
    }
}


std::vector<ClusterID> ReadClusteringEngine::merge_clusters(std::vector<IDComponent> &components) {
    auto component_queue = ConcurrentQueue<IDComponent>(components.begin(), components.end());

    IndexRemovalMap for_removal;
    std::vector<ClusterID> merge_result_ids;
    run_in_threads(&ReadClusteringEngine::merge_clusters_thread, this, std::ref(component_queue), std::ref(for_removal), std::ref(merge_result_ids));

    auto removal_list_queue = ConcurrentQueue<IndexRemovalMap::value_type>(for_removal.begin(), for_removal.end());
    run_in_threads(&ReadClusteringEngine::kmer_cluster_index_update, this, std::ref(removal_list_queue));
    return merge_result_ids;
}

ClusterID get_parent(ClusterID cluster_id, tsl::robin_map<ClusterID, ClusterID> &parents) {
    if (cluster_id == parents[cluster_id]) {
        return cluster_id;
    }
    parents[cluster_id] = get_parent(parents[cluster_id], parents);
    return parents[cluster_id];
}


std::vector<IDComponent> ReadClusteringEngine::union_find(std::vector<ClusterConnection> &connections, tsl::robin_set<ClusterIDPair> &restricted) {
    tsl::robin_map<ClusterID, ClusterID> parents;
    tsl::robin_map<ClusterID, IDComponent> components;
    for (auto cluster_iter = begin(cluster_index); cluster_iter != end(cluster_index); cluster_iter++) {
        parents[cluster_iter->second->id] = cluster_iter->second->id;
        components[cluster_iter->second->id] = {cluster_iter->second->id};
    }

    for (ClusterConnection &conn : connections) {
        ClusterID parent_x = get_parent(conn.cluster_x_id, parents);
        ClusterID parent_y = get_parent(conn.cluster_y_id, parents);
        if (parent_x == parent_y) continue;
        if (restricted.contains(get_cluster_id_pair(parent_x, parent_y))) continue;

        ClusterID bigger, smaller;
        if (components[parent_x].size() > components[parent_y].size()) {
            bigger = parent_x;
            smaller = parent_y;
        } else {
            bigger = parent_y;
            smaller = parent_x;
        }

        for (ClusterID cluster_id : components[smaller]) {
            parents[cluster_id] = bigger;
        }
        components[bigger].insert(components[bigger].end(), components[smaller].begin(), components[smaller].end());
        components.erase(smaller);
    }

    std::vector<IDComponent> result;
    for (const auto& comp : components){
        if (comp.second.size() > 1) result.push_back(comp.second);
    }
    return result;
}

std::map<ClusterID, std::string> ReadClusteringEngine::export_clusters(std::vector<ClusterID> &cluster_ids, fs::path &directory_path) {
    tsl::robin_map<ReadID, GenomeReadData> id_to_read;
    for (auto cluster_id : cluster_ids) {
        for (auto read_meta: cluster_index[cluster_id]->contained_reads) {
            id_to_read[read_meta.id] = {};
        }
    }
    fs::remove_all(directory_path);
    fs::create_directories(directory_path);
    reader->rewind();
    std::optional<GenomeReadData> read;
    while ((read = reader->get_next_record()) != std::nullopt) {
        if (id_to_read.contains(read->id)) {
            id_to_read[read->id] = *read;
        }
    }

    std::map<ClusterID, std::string> mapping;
    for (auto id : cluster_ids) {
        std::ofstream cluster_file;
        std::string cluster_file_path = fmt::format("{}/#{}.fq", directory_path.string(), id);

        cluster_file.open(cluster_file_path);
        for (auto read_meta : cluster_index[id]->contained_reads) {
            cluster_file << id_to_read[read_meta.id].fastX_string() << std::endl;
        }
        cluster_file.close();
        mapping[id] = cluster_file_path;
    }

    return mapping;
}

void ReadClusteringEngine::assemble_clusters(std::vector<ClusterID> &cluster_ids) {
    sort(cluster_ids.rbegin(), cluster_ids.rend(), [this](ClusterID x, ClusterID y) -> bool { return cluster_index[x]->size() < cluster_index[y]->size(); });

    fs::path export_path = "./data/assembly_" + reader->meta.filename;
    auto mapping = export_clusters(cluster_ids, export_path);

    for (ClusterID id : cluster_ids) {
        std::string assembled_fasta_path = mapping[id] + "_assembled.fasta";

        std::string cmd = fmt::format("./scripts/assemble_pacbio_cluster.sh {} {} {}", mapping[id], assembled_fasta_path, platform == PacBio ? "pb" : "ont");
        run_command_with_input(cmd.c_str(), "");

        std::cout << cluster_index[id]->to_string() << " assembly: \n";
        try {
            SequenceRecordIterator r(assembled_fasta_path);
            r.show_progress = false;
            std::optional<GenomeReadData> contig;
            while ((contig = r.get_next_record()) != std::nullopt) {
                std::cout << fmt::format("\tContig size {}\n", contig->sequence.length());
            }
        } catch (const std::logic_error &e) {
            std::cout << "\t No contigs assembled\n";
        }
    }
}

void ReadClusteringEngine::run_clustering() {
    std::vector<ClusterConnection> cluster_connections = timeMeasureMemberFunc(&ReadClusteringEngine::get_all_connections, this, "Connection calculation")(1);
    plot_connection_quality(cluster_connections);
    int first_round_min_score;
    std::cin >> first_round_min_score;
    cluster_connections = filter_connections(cluster_connections, [first_round_min_score](ClusterConnection &conn){ return conn.score >= first_round_min_score; });

    tsl::robin_set<ClusterIDPair> restricted;
    auto chain_components = timeMeasureMemberFunc(&ReadClusteringEngine::union_find, this, "Union find")(cluster_connections, restricted);
    std::vector<ClusterID> chain_ids = timeMeasureMemberFunc(&ReadClusteringEngine::merge_clusters, this, "Cluster merging")(chain_components);

    print_clusters(-1);

    for (auto cluster_id : chain_ids){
        for (auto cluster_id_2 : chain_ids){
            restricted.insert(get_cluster_id_pair(cluster_id, cluster_id_2));
        }
    }

    while (true){
        cluster_connections = get_connections(chain_ids, 1);
        plot_connection_quality(cluster_connections);

        int refine_score;
        std::cin >> refine_score;
        if (refine_score == 0) break;

        cluster_connections = filter_connections(cluster_connections, [refine_score](ClusterConnection &conn){ return conn.score >= refine_score; });
        chain_components = union_find(cluster_connections, restricted);
        merge_clusters(chain_components);

        print_clusters(-1);
    }

    std::vector<GenomeReadCluster *> clusters_for_display;
    for (auto cluster_id : chain_ids) clusters_for_display.push_back(cluster_index[cluster_id]);
    plot_cluster_coverage(clusters_for_display);

    std::vector<ClusterID> clusters_for_assembly = filter_clusters([](GenomeReadCluster* c){ return (c->size() > 20); });
    assemble_clusters(clusters_for_assembly);
}

ReadClusteringEngine::~ReadClusteringEngine() {
    for (auto cluster_pair : cluster_index){
        free(cluster_pair.second);
    }
}
