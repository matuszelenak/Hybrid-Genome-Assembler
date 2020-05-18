Connection = namedtuple(x=ComponentID, y=ComponentID, strength=int)

def get_connections(from_component: ComponentID):
	shared_kmer_counts : Dict[ComponentID, int]
	for kmer in components[from_component].discriminative:
		for component_id in kmer_component_index[kmer]:
			if component_id in shared_kmer_counts:
				shared_kmer_counts[component_id] += 1
			else:
				shared_kmer_counts[component_id] = 1

	del shared_kmer_counts[from_component]

	result = []
	for component_id, shared_count in shared_kmer_counts.items():
		result.append(from_component, component_id, shared_count)

	return result
