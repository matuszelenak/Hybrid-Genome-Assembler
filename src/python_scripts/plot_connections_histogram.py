import math
from typing import Dict

import numpy as np
from matplotlib import pyplot as plt

good_conns: Dict[int, int] = eval(input())
bad_conns: Dict[int, int] = eval(input())

strengths = sorted(list(list(good_conns.keys()) + list(bad_conns.keys())), reverse=True)

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
indices = np.array(strengths)
good_data = np.array([(math.log(good_conns[strength]) if strength in good_conns else 0) for strength in strengths])
plt.bar(indices, good_data, width=0.8, color='g')

bad_data = np.array([(math.log(bad_conns[strength]) if strength in bad_conns else 0) for strength in strengths])
plt.bar(indices, bad_data, bottom=good_data, width=0.8, color='r')

plt.show()
