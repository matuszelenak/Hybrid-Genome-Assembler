Connection = namedtuple(x=ClusterID, y=ClusterID, score=int)

def get_connections(from_cluster: ClusterID):
	shared_kmer_counts : Dict[ClusterID, int]
	for kmer in clusters[from_cluster].discriminative:
		for cluster_id in kmer_cluster_index[kmer]:
			if cluster_id in shared_kmer_counts:
				shared_kmer_counts[cluster_id] += 1
			else:
				shared_kmer_counts[cluster_id] = 1

	del shared_kmer_counts[from_cluster]

	result = []
	for cluster_id, shared_count in shared_kmer_counts.items():
		result.append(from_cluster, cluster_id, shared_count)

	return result
