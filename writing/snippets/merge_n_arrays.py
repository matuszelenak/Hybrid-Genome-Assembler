def merge_sorted_arrays(arrays: List[List[Element]], unique: bool):
	heap = PriorityQueue()
	next_position_for_array = []
	for i, arr in enumerate(arrays):
		heap.push((arr[0], i))

	not_exhausted = len(arrays)
	previous_element = None
	result = []
	while not_exhausted > 0:
		element, array_index = heap.pop()

		if not unique or element != previous_element:
			result.append(element)
			previous_element = element

		if next_position_for_array[array_index] == len(arrays[array_index]):
			non_exhausted -= 1
		else:
			next_element = arrays[array_index][next_position_for_array[array_index]]
			heap.push((next_element, array_index))

	return result
