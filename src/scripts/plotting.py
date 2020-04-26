import argparse
import math
from typing import List, Tuple, Dict

import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

BoundSpecificities = Tuple[int, List[Tuple[int, int]]]
KmerSpecificities = Tuple[int, List[BoundSpecificities]]


class Unimplemented:
    def __init__(self):
        raise NotImplementedError


class KmerHistogram:
    def __init__(self):
        occurrences: Dict[int, int] = eval(input())

        for i in range(max(occurrences.keys())):
            if i not in occurrences:
                occurrences[i] = 0

        indices = np.arange(1, max(occurrences.keys()) + 1)
        data = np.array([occurrences.get(coverage, 0) for coverage in indices])
        bars = plt.bar(indices, data, width=0.8, color='r')

        peak_indices = self.find_peak_positions(occurrences)

        for _, index in sorted([(occurrences[i], i) for i in peak_indices], reverse=True)[:3]:
            height = bars[index].get_height()
            plt.text(bars[index].get_x() + bars[index].get_width() / 2.0, height, f'{height}', ha='center', va='bottom')

        plt.show()

    @staticmethod
    def find_peak_positions(occurrences: Dict[int, int]):
        occ_values = [occurrences[coverage] for coverage in sorted(occurrences)]
        peak_indices, props = find_peaks(np.array(occ_values))
        return peak_indices


class KmerHistogramWithSpec:
    def __init__(self):
        num_of_k, max_cov = [int(x) for x in input().split()]
        specs: List[KmerSpecificities] = [eval(input()) for _ in range(num_of_k)]  # FOR THE GLORY OF SATAN
        self.plot_histograms(specs, max_cov)

    @staticmethod
    def render_subplot(subplot, bound_specificities: List[BoundSpecificities], max_coverage: int):
        get_color = dict((
            (70, 'k'),
            (85, 'r'),
            (90, 'xkcd:orange'),
            (95, 'y'),
            (100, 'b'),
            (101, 'g')
        )).get

        skip_first = 5

        bars = []
        indices = np.arange(1, max_coverage + 1)
        bottoms = np.zeros((max_coverage,))
        for upper_bound, coverage_data in bound_specificities:
            coverage_dict = dict(coverage_data)
            specificity_level_data = np.array([0 for _ in range(skip_first)] + [
                coverage_dict.get(coverage, 0) for coverage in range(skip_first + 1, max_coverage + 1)
            ])
            bars.append(subplot.bar(indices, specificity_level_data, bottom=bottoms, width=0.8, color=get_color(upper_bound)))
            bottoms += specificity_level_data

        bounds = [x[0] for x in bound_specificities]
        subplot.legend(bars, [f'{low} ≤ x < {high}' for low, high in zip([50] + bounds, bounds)])

    @staticmethod
    def plot_histograms(specificities: List[KmerSpecificities], max_coverage: int):
        if len(specificities) > 1:
            fig, axs = plt.subplots(len(specificities), sharex='col', sharey='row', num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')

            for row, (k_value, data) in enumerate(specificities):
                axs[row].set_title(f'Occurrences of characteristic {k_value}-mers in reads')
                KmerHistogramWithSpec.render_subplot(axs[row], data, max_coverage)
        else:
            for row, (k_value, data) in enumerate(specificities):
                plt.title(f'Occurrences of characteristic {k_value}-mers in reads')
                KmerHistogramWithSpec.render_subplot(plt, data, max_coverage)
        plt.show()


class ConnectionHistogram:
    def __init__(self):
        good_conns: Dict[int, int] = eval(input())
        bad_conns: Dict[int, int] = eval(input())

        strengths = sorted(list(list(good_conns.keys()) + list(bad_conns.keys())), reverse=True)

        indices = np.array(strengths)
        good_data = np.array([(math.sqrt(good_conns[strength]) if strength in good_conns else 0) for strength in strengths])
        plt.bar(indices, good_data, width=0.8, color='g')

        bad_data = np.array([(math.sqrt(bad_conns[strength]) if strength in bad_conns else 0) for strength in strengths])
        plt.bar(indices, bad_data, bottom=good_data, width=0.8, color='r')

        plt.show()


class ClusterCoverage:
    def __init__(self):
        a_intervals = sorted(eval(input()))
        b_intervals = sorted(eval(input()))

        def plot_intervals(intervals, y):
            for start, stop in intervals:
                plt.hlines(y, start, stop, (0.0, 1.0, 0.0, 0.3), lw=20)
                plt.vlines(start, y+0.03, y-0.03, 'g')
                plt.vlines(stop, y+0.03, y-0.03, 'r')

        plot_intervals(a_intervals, -0.1)
        plot_intervals(b_intervals, 0.1)

        plt.show()


SUPPORTED_PLOTS = {
    'kmer_histogram': KmerHistogram,
    'kmer_histogram_with_spec': KmerHistogramWithSpec,
    'connection_histogram': ConnectionHistogram,
    'cluster_coverage': ClusterCoverage
}

parser = argparse.ArgumentParser()
parser.add_argument('--plot', dest='plot_choice', type=str, choices=SUPPORTED_PLOTS, required=True)
args = parser.parse_args()

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
SUPPORTED_PLOTS.get(args.plot_choice, Unimplemented)()
