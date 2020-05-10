type Vertex = int
type Distance = int

def distance_bfs(start: Vertex, adjacency: Dict[Vertex, Dict[Vertex, Distance]]) -> Dict[Vertex, Distance]:
	vertex_queue = Queue()
	vertex_queue.push(starting)
	visited = set()
	distances = {start: read_length[start]}

	while (!vertex_queue.empty()):
		vertex = vertex_queue.pop()
		visited.add(vertex)

		if vertex not in adjacency:
			continue

		for neighbor, distance in adjacency[vertex].items():
			if neighbor in visited:
				continue

			distances[neighbor] = distances[vertex] + read_length[neighbor] - 2 * overlap_size(neighbor, vertex)
			vertex_queue.push(neighbor)

	return distances
