def spanning_tree_tails(adjacency: Dict[Vertex, Dict[Vertex, Distance]]):
	initial_distances = distance_bfs(adjacency.keys()[0], adjacency)
	farthest = max(initial_distances.items(), key=lambda key_value: key_value[1])

	distance_from_left = distance_bfs(farthest[0], adjacency)
	farthest_from_left = max(distance_from_left.items(), key=lambda key_value: key_value[1])

	distance_from_right = distance_bfs(farthest_from_left[0], adjacency)
	farthest_from_right = max(distance_from_right.items(), key=lambda key_value: key_value[1])

	tail_length = farthest_from_right[1] * 0.03

	left_tail_vertices, right_tail_vertices = [], []
	for vertex, distance in distance_from_right.items():
		if distance + tail_length > farthest_from_right:
			left_tail_vertices.append(vertex)

	for vertex, distance in distance_from_left.items():
		if distance + tail_length > farthest_from_left:
			right_tail_vertices.append(vertex)

	return left_tail_vertices, right_tail_vertices