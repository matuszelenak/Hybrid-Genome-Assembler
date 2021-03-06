import argparse
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
            (100.01, 'g')
        )).get

        skip_first = 3

        bars = []
        indices = np.arange(1, max_coverage + 1)
        bottoms = np.zeros((max_coverage,))
        for upper_bound, coverage_data in bound_specificities:
            coverage_dict = dict(coverage_data)
            specificity_level_data = np.array([0 for _ in range(skip_first)] + [
                coverage_dict.get(coverage, 0) for coverage in range(skip_first + 1, max_coverage + 1)
            ])
            bars.append(subplot.bar(indices, specificity_level_data, bottom=bottoms, width=1, color=get_color(upper_bound)))
            bottoms += specificity_level_data

        bounds = [x[0] for x in bound_specificities]
        subplot.legend(bars, [f'{int(low)}% ≤ S < {int(high)}%' for low, high in zip([50] + bounds, bounds)])
        subplot.xlabel("Occurrence of kmers")
        subplot.ylabel("Number of unique kmers with occurrence")

    @staticmethod
    def plot_histograms(specificities: List[KmerSpecificities], max_coverage: int):
        if len(specificities) > 1:
            fig, axs = plt.subplots(len(specificities), sharex='col', sharey='row', num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')

            for row, (k_value, data) in enumerate(specificities):
                axs[row].set_title(f'Occurrence histogram of {k_value}-mers in reads')
                KmerHistogramWithSpec.render_subplot(axs[row], data, max_coverage)
        else:
            for row, (k_value, data) in enumerate(specificities):
                plt.title(f'Occurrence histogram of {k_value}-mers in reads')
                KmerHistogramWithSpec.render_subplot(plt, data, max_coverage)
        plt.show()


class ConnectionHistogram:
    def __init__(self):
        good_conns: Dict[int, int] = eval(input())
        bad_conns: Dict[int, int] = eval(input())

        strengths = list(set(sorted(list(good_conns.keys()) + list(bad_conns.keys()))))

        indices = np.array(strengths[2:400])
        plt.subplot("211")
        good_data = np.array([((good_conns[strength]) if strength in good_conns else 0) for strength in strengths[2:400]])
        g = plt.bar(indices, good_data, width=0.8, color='g')

        bad_data = np.array([((bad_conns[strength]) if strength in bad_conns else 0) for strength in strengths[2:400]])
        b = plt.bar(indices, bad_data, bottom=good_data, width=0.8, color='r')
        plt.legend([g, b], ["Correct connections", "Incorrect connections"])
        plt.xlabel("Connection strength")
        plt.ylabel("Number of connections")

        probabilities = []
        for strength in strengths:
            probabilities.append(
                (good_conns.get(strength, 0)) / (good_conns.get(strength, 0) + bad_conns.get(strength, 0))
            )
        plt.subplot("212")
        plt.hlines(1.0, 0, strengths[99], linestyle='dotted', linewidth=2)
        plt.plot(np.array(strengths[:100]), np.array(probabilities[:100]), color="green", linewidth=2)
        plt.xlabel("Connection strength")
        plt.ylabel("Proportion of correct connections")
        plt.show()

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


class ConnectionIntervalOverlap:
    def __init__(self):
        scores_values = int(input())
        scores = []
        scores_avg_overlap = {}
        scores_min_overlap = {}
        for _ in range(scores_values):
            score, overlap_map = eval(input())
            s = 0
            c = 0
            scores_min_overlap[score] = 10000000000
            for overlap, count in overlap_map.items():
                s += (overlap * count)
                c += count

                scores_min_overlap[score] = min(overlap, scores_min_overlap[score])

            scores.append(score)
            scores_avg_overlap[score] = s / c if c > 0 else 0

        scores.sort()
        index = 0
        step = 5
        y_avg, y_min = [], []
        x_values = []
        while True:
            sum_of_min = 0
            sum_of_avg = 0
            for score_index in range(index, min(len(scores), index + step)):
                sum_of_avg += scores_avg_overlap[scores[score_index]]
                sum_of_min += scores_min_overlap[scores[score_index]]

            y_avg.append(sum_of_avg / min(step, len(scores) - index))
            y_min.append(sum_of_min / min(step, len(scores) - index))

            x_values.append(scores[index])

            if index + step >= len(scores):
                break

            index += step
            if index >= len(scores):
                break

            step = int(step * 1.3)

        plt.plot(np.array(x_values), np.array(y_avg), color="blue", label="Average overlap")
        plt.plot(np.array(x_values), np.array(y_min), color="green", label="Minimal overlap")
        plt.xlabel("Connection strength")
        plt.ylabel("Read Overlap size")
        plt.legend(loc="upper left")
        plt.title("Relation of read connection strength and overlaps")
        plt.show()


class OverlapErrorHisto:
    def __init__(self):
        import math
        errors = []
        while True:
            try:
                i = input()
            except EOFError:
                break

            calculated, real = [int(x) for x in i.split()]
            error = ((calculated - real) / ((calculated + real) / 2))
            errors.append(math.floor(error * 100))

        mi = min(errors)
        ma = max(errors)
        tick_values = np.arange(mi, ma, (ma - mi) / 30)

        plt.hist(errors, bins=300)
        plt.title("Estimated read overlaps vs real read overlaps")
        plt.xticks(tick_values, ['{:.2f}%'.format(err / 100) for err in tick_values])
        plt.xlabel("Difference rate")
        plt.show()


def SpectralGraph():
    edges = eval(input())
    edges = [(u, v, w) for u, v, w in edges if w > 50]
    categories = dict(set(eval(input())))
    components = eval(input())

    import networkx as nx
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    options = {
        'node_size': 250,
        'font_color': "white"
    }
    color_map = []
    for node in G:
        if categories[node] == 0:
            color_map.append("blue")
        elif categories[node] == 1:
            color_map.append("black")
        else:
            color_map.append("red")

    component_colors = ["green", "orange", "yellow", "red", "black", "cyan"]
    in_component_colors = {}
    for i, component in enumerate(components):
        for vertex in component:
            in_component_colors[vertex] = component_colors[i]

    plt.subplot(121)
    layout = nx.spring_layout(G, scale=2)
    nx.draw(G, layout, **options, edges=edges, node_color=color_map)
    nx.draw_networkx_edge_labels(G, layout, edge_labels={(u, v): G[u][v]['weight'] for u, v in G.edges})

    plt.subplot(122)
    nx.draw(G, layout, **options, edges=edges, node_color=[in_component_colors[vertex] for vertex in G])
    nx.draw_networkx_edge_labels(G, layout, edge_labels={(u, v): G[u][v]['weight'] for u, v in G.edges})
    plt.show()


def ComponentCoverage():
    components = eval(input())
    for i, c in enumerate(components):
        plt.subplot(f'{len(components)}1{i + 1}')
        plt.bar(np.array([x[0] for x in c[1:]]), np.array([x[1] for x in c[1:]]), width=0.8, color='g')
        plt.xlabel("Coverage")
        plt.ylabel("Bases with coverage")

    plt.show()


SUPPORTED_PLOTS = {
    'kmer_histogram': KmerHistogram,
    'kmer_histogram_with_spec': KmerHistogramWithSpec,
    'connection_histogram': ConnectionHistogram,
    'cluster_coverage': ClusterCoverage,
    'connection_overlaps': ConnectionIntervalOverlap,
    'overlap_errors': OverlapErrorHisto,
    'spectral_graph': SpectralGraph,
    'component_coverage': ComponentCoverage
}

parser = argparse.ArgumentParser()
parser.add_argument('--plot', dest='plot_choice', type=str, choices=SUPPORTED_PLOTS, required=True)
args = parser.parse_args()

plt.figure(num=None, figsize=(1920 // 80, 1080 // 80), dpi=80, facecolor='w', edgecolor='k')
plt.rcParams.update({'font.size': 22})
SUPPORTED_PLOTS.get(args.plot_choice, Unimplemented)()
