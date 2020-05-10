def approximate_overlap_length(Read x, Read y) -> int:
	shared_kmers = intersection(x.discriminative, y.discriminative)

	x_positions, y_positions = [], []
	for kmer in shared_kmers:
		x_positions.append(kmer_positions[x][kmer])
		y_positions.append(kmer_positions[y][kmer])

	max_x = max(x_positions)
	min_y = min(x_positions)
	max_y = max(y_positions)
	min_y = min(y_positions)

	return max(max_x - min_y, max_y - min_y)
