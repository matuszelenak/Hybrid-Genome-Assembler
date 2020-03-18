import argparse
import os

from Bio import SeqIO
import matplotlib.pyplot as plt

from .constants import *


class ReadCategorizer:

    def __init__(self, path, k):
        self.parent = []
        self.components = []
        self.component_kmers = []

        self.categorize_reads(path, k)

    def get_parent(self, vertex):
        if self.parent[vertex] == vertex:
            return vertex

        self.parent[vertex] = self.get_parent(self.parent[vertex])
        return self.parent[vertex]

    def unite(self, x, y):
        parent_x = self.get_parent(x)
        parent_y = self.get_parent(y)
        if parent_x == parent_y:
            return False

        if len(self.components[parent_x]) > len(self.components[parent_y]):
            bigger = parent_x
            smaller = parent_y
        else:
            bigger = parent_y
            smaller = parent_x

        for vertex in self.components[smaller]:
            self.parent[vertex] = bigger

        self.components[bigger] += self.components[smaller]
        self.components[smaller] = []

        self.component_kmers[bigger] = self.component_kmers[bigger] | self.component_kmers[smaller]
        self.component_kmers[smaller] = set()

        return True

    def histogram(self, data):
        plt.figure()
        plt.title('Occurence histogram')
        plt.xlabel('Coverage')
        plt.ylabel('Frequency')
        plt.hist(data, color='blue', bins=len(set(data)))
        plt.show()

    def component_label(self, component, labels, label_kinds):
        percentages = []
        for label_id in label_kinds.values():
            percentages.append((
                int(sum([100 for v in self.components[component] if labels[v] == label_id]) / len(self.components[component])),
                label_id
            ))
        percentages.sort(reverse=True)
        return "{}%".format(percentages[0][0]), percentages[0][1]

    def component_stats(self, labels, label_kinds):
        r = [(len(x), i) for i, x in enumerate(self.components) if len(x) > 0]
        r.sort(reverse=True)

        print('|'.join(
            [str(self.component_label(i, labels, label_kinds)) for _, i in r[:30]]
        ))

    def categorize_reads(self, path, k):
        base_path, extension = os.path.splitext(path)
        if extension in ('.fasta', '.fa'):
            extension = 'fasta'
        else:
            extension = 'fastq'

        kmer_counts = {}
        for seq_record in SeqIO.parse(path, extension):
            sequence = str(seq_record.seq)
            for i in range(k, len(sequence) + 1):
                forward_kmer = sequence[i - k: i]
                backward_kmer = ''.join([complement[x] for x in forward_kmer])[::-1]
                signature = min(forward_kmer, backward_kmer)
                if signature not in kmer_counts:
                    kmer_counts[signature] = 0
                kmer_counts[signature] = kmer_counts[signature] + 1

        self.histogram([x for _, x in kmer_counts.items() if 3 < x < 100])
        print("Enter lower and upper bound")
        lower, upper = [int(x) for x in input().split()]
        kmers = set([sequence for sequence, count in kmer_counts.items() if lower <= count <= upper])

        kmer_to_number = {kmer: i for i, kmer in enumerate(kmers)}

        # for each read calculate its own set of characteristic kmers
        read_kmer_sets = []
        labels = []
        label_kinds = {}
        for seq_record in SeqIO.parse(path, extension):
            sequence = str(seq_record.seq)
            label = seq_record.id
            if label not in label_kinds:
                label_kinds[label] = len(label_kinds)
            labels.append(label_kinds[label])

            kmer_set = set()
            for i in range(k, len(sequence) + 1):
                forward_kmer = sequence[i - k: i]
                backward_kmer = ''.join([complement[x] for x in forward_kmer])[::-1]
                signature = min(forward_kmer, backward_kmer)
                if signature in kmers:
                    kmer_set.add(kmer_to_number[signature])

            read_kmer_sets.append(kmer_set)

        edges = []
        for i, A in enumerate(read_kmer_sets):
            for j, B in enumerate(read_kmer_sets):
                if i == j:
                    continue

                edge_weight = len(A & B)
                if edge_weight == 0:
                    continue

                edges.append((edge_weight, i, j))

        edges.sort(reverse=True)

        print("Label kinds: ", label_kinds)
        print("Labels: ", labels[:30])
        print('Constructed {} nodes'.format(len(read_kmer_sets)))
        print("Constructed {} edges".format(len(edges)))

        self.parent = [x for x in range(len(read_kmer_sets))]
        self.components = [[x] for x in range(len(read_kmer_sets))]
        self.component_kmers = [set(x for x in s) for s in read_kmer_sets]

        comp_size_treshold = 30
        for edge in edges:
            if edge[0] < 80:
                break
            if self.unite(edge[1], edge[2]):
                print(edge[0])
                self.component_stats(labels, label_kinds=label_kinds)
                # i = input()
                # if i:
                #     comp_size_treshold = int(i)
                #     break

        # Second round of union-find on clusters
        edges = []
        filtered_components = [(i, self.component_kmers[i]) for i, x in enumerate(self.components) if len(x) > comp_size_treshold]
        for i, A in filtered_components:
            for j, B in filtered_components:
                if i == j:
                    continue

                edge_weight = len(A & B)
                print("Components of sizes {} {}".format(len(self.components[i]), len(self.components[j])))
                print("Component types {} and {}".format(
                    self.component_label(i, labels, label_kinds),
                    self.component_label(j, labels, label_kinds)
                ))
                print("Components with {} and {} kmers share {} kmers".format(
                    len(A), len(B), edge_weight
                ))
                input()
                if edge_weight == 0:
                    continue

                edges.append((edge_weight, i, j))

        edges.sort(reverse=True)

        for edge in edges:
            if self.unite(edge[1], edge[2]):
                print(edge[0])
                self.component_stats(labels, label_kinds=label_kinds)
                input()

        # for c in self.big_components:
        #     hist_data = []
        #     for v in self.components[c]:
        #         hist_data.append(labels[v])
        #     self.histogram(hist_data)


parser = argparse.ArgumentParser()
parser.add_argument('path', type=str,
                    help='Path to a read file')
parser.add_argument('-k', type=int,
                    help='Length of kmer')

args = parser.parse_args()
c = ReadCategorizer(args.path, args.k)
