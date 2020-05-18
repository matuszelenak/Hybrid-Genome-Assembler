def get_unique_kmer_length():
	k = 11
	previous_kmer_count = HyperLogLog(k)

	while k < 33:
		next_count = HyperLogLog(k + 2)
		if difference(previous_kmer_count, next_count) < 0.2:
			return k

		previous_kmer_count = next_count
		k += 2

	return k