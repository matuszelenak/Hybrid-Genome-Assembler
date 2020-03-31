import json
from typing import Dict

import numpy as np
from matplotlib import pyplot as plt

occurrences: Dict[int, int] = eval(input())
print(json.dumps(occurrences, indent=4))

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
indices = np.arange(max(occurrences.keys()))

bad_data = np.array([occurrences.get(coverage, 0) for coverage in indices])
plt.bar(indices, bad_data, width=0.8, color='r')

plt.show()
