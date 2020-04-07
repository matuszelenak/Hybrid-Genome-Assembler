import json
from typing import Dict

import numpy as np
from matplotlib import pyplot as plt

occurrences: Dict[int, int] = eval(input())

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
indices = np.arange(max(occurrences.keys()))

data = np.array([occurrences.get(coverage, 0) for coverage in indices])
plt.bar(indices, data, width=0.8, color='r')

plt.show()
