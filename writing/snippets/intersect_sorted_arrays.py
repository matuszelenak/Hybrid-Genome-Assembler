def intersect_sorted_arrays(x, y):
	x_index, y_index = 0, 0
	while x_index < len(x) and y_index < len(y):
		if x[x_index] < y[y_index]:
			x_index += 1
		elif y[y_index] < x[x_index]:
			y_index += 1
		else:
			result.append(x[x_index])
			x_index += 1
			y_index += 1
			
	return result