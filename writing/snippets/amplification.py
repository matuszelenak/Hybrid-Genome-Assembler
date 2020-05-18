def amplified_tail_SDKs(tail: List[Vertex]) -> Set[Kmer]:
	amplified_tail = tail[::]
	for connection in get_connections(from=tail, minimal_score=40):
		amplified_tail.append(connection.to)

	amplified_tail_SDKs = {}
	for vertex in amplified_tail:
		amplified_tail_SDKs.add(vertex.discriminative_kmers)

	return amplified_tail_SDKs