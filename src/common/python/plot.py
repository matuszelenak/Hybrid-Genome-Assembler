from typing import Dict, List

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks


def find_peak_positions(occurrences: Dict[int, int]):
    occ_values = [occurrences[coverage] for coverage in sorted(occurrences)]
    print(occ_values)
    peak_indices, props = find_peaks(np.array(occ_values))
    print(peak_indices)
    print(props)
    return peak_indices


occurrences: Dict[int, int] = eval(input())
for i in range(max(occurrences.keys())):
    if i not in occurrences:
        occurrences[i] = 0

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
indices = np.arange(max(occurrences.keys()))

data = np.array([occurrences[coverage] for coverage in indices])
bars = plt.bar(indices, data, width=0.8, color='r')

peak_indices = find_peak_positions(occurrences)

for _, index in sorted([(occurrences[i], i) for i in peak_indices], reverse=True)[:3]:
    height = bars[index].get_height()
    plt.text(bars[index].get_x() + bars[index].get_width() / 2.0, height, f'{height}', ha='center', va='bottom')

plt.show()
