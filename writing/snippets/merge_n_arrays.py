def merge_sorted_arrays(arrays: List[List[Element]], unique: bool):
	heap = PriorityQueue()
	next_position_for_array = [0 for _ in range(len(arrays))]
	for i, arr in enumerate(arrays):
		if next_position_for_array[array_index] < len(arrays[array_index]):
			heap.push((arr[next_position_for_array[i]], i))
			next_position_for_array[i] += 1

	previous_element = None
	result = []
	while not heap.empty():
		element, array_index = heap.pop()

		if not unique or element != previous_element:
			result.append(element)
			previous_element = element

		if next_position_for_array[array_index] < len(arrays[array_index]):
			next_element = arrays[array_index][next_position_for_array[array_index]]
			heap.push((next_element, array_index))
			next_position_for_array[array_index] += 1

	return result
