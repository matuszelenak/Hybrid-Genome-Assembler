type Position = int
type IsStart = bool
type Interval = namedtuple(start=int, end=int)

def overlap_intervals(intervals: List[Interval]) -> List[Interval]:
	interval_endpoints: List[Tuple[Position, IsStart]]
	for interval in intervals:
		interval_endpoints.append((read.start, True))
		interval_endpoints.append((read.start, False))
	interval_endpoints.sort()

	overlapped_intervals = []

	opened_interval_start = -1
	opened_intervals_count = 0
	for position, is_start in interval_endpoints:
		if is_start:
			opened_intervals_count += 1
			if opened_intervals_count == 1:
				opened_interval_start = position
		else:
			opened_intervals_count -= 1
			if opened_intervals_count == 0:
				overlapped_intervals.append(
					Interval(opened_interval_start, position)
					)

	return overlapped_intervals