from typing import List, Tuple

import numpy as np
from matplotlib import pyplot as plt

BoundSpecificities = Tuple[int, List[Tuple[int, int]]]
KmerSpecificities = Tuple[int, List[BoundSpecificities]]


get_color = dict((
    (70, 'k'),
    (85, 'r'),
    (90, 'xkcd:orange'),
    (95, 'y'),
    (100, 'b'),
    (101, 'g')
)).get


def render_subplot(subplot, bound_specificities: List[BoundSpecificities], max_coverage: int):
    bars = []
    indices = np.arange(1, max_coverage + 1)
    bottoms = np.zeros((max_coverage,))
    for upper_bound, coverage_data in bound_specificities:
        coverage_dict = dict(coverage_data)
        specificity_level_data = np.array([
            coverage_dict.get(coverage, 0) for coverage in range(1, max_coverage + 1)
        ])
        bars.append(subplot.bar(indices, specificity_level_data, bottom=bottoms, width=0.8, color=get_color(upper_bound)))
        bottoms += specificity_level_data

    bounds = [x[0] for x in bound_specificities]
    subplot.legend(bars, [f'{low} ≤ x < {high}' for low, high in zip([50] + bounds, bounds)])


def plot_histograms(specificities: List[KmerSpecificities], max_coverage: int):
    fig, axs = plt.subplots(len(specificities), sharex='col', sharey='row', num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')

    for row, (k_value, data) in enumerate(specificities):
        axs[row].set_title(f'Occurrences of characteristic {k_value}-mers in reads')
        render_subplot(axs[row], data, max_coverage)

    plt.show()


files = input().split('/')
num_of_k, max_cov = [int(x) for x in input().split()]
specs: List[KmerSpecificities] = [eval(input()) for _ in range(num_of_k)]  # FOR THE GLORY OF SATAN
plot_histograms(specs, max_cov)
# plt.savefig(f'histogram_{",".join(files)}_k=[{",".join(map(str, [k for k, _ in specs]))}].png')
