import argparse
import os

from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt

from constants import *


class ReadCategorizer:

    def __init__(self, path, k):
        self.parent = []
        self.components = []

        self.big_components = []
        self.max_big = 2

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

        # if len(self.components[bigger]) > 50 and bigger not in self.big_components and len(self.big_components) < self.max_big:
        #    self.big_components.append(bigger)

        return True

    def histogram(self, data):
        plt.figure()
        plt.title('Occurence histogram')
        plt.xlabel('Coverage')
        plt.ylabel('Frequency')
        plt.hist(data, color='blue', bins=len(set(data)))
        plt.show()

    def components_by_size(self):
        r = [(len(x), i) for i, x in enumerate(self.components)]
        r.sort(reverse=True)
        return r

    def categorize_reads(self, path, k):
        print(path)
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

        self.histogram([x for _, x in kmer_counts.items() if 3 < x < 70])
        print("Enter lower and upper bound")
        lower, upper = [int(x) for x in input().split()]
        kmers = set([sequence for sequence, count in kmer_counts.items() if lower <= count <= upper])

        # for each read calculate its own set of characteristic kmers
        read_kmer_sets = []
        labels = []
        for seq_record in SeqIO.parse(path, extension):
            sequence = str(seq_record.seq)
            labels.append(seq_record.id[:30])
            read_kmer_sets.append(set())
            for i in range(k, len(sequence) + 1):
                forward_kmer = sequence[i - k: i]
                backward_kmer = ''.join([complement[x] for x in forward_kmer])[::-1]
                signature = min(forward_kmer, backward_kmer)
                if signature in kmers:
                    read_kmer_sets[-1].add(signature)

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

        print('Constructed {} nodes'.format(len(read_kmer_sets)))
        print("Constructed {} edges".format(len(edges)))

        self.parent = [x for x in range(len(read_kmer_sets))]
        self.components = [[x] for x in range(len(read_kmer_sets))]
        num_of_components = len(read_kmer_sets)

        G = nx.Graph()
        G.add_nodes_from(range(len(read_kmer_sets)))
        for edge in edges:

            if self.get_parent(edge[1]) in self.big_components and self.get_parent(edge[2]) in self.big_components:
                continue

            if self.unite(edge[1], edge[2]):
                # G.add_edge(edge[1], edge[2])
                # nx.draw(G)
                # plt.show()
                print(edge[0])
                print([x for x, _ in self.components_by_size()][:30])
                input()
                # print(self.big_components)
                # input()
                num_of_components -= 1
                if num_of_components == 2:
                    break

        for c in self.big_components:
            hist_data = []
            for v in self.components[c]:
                hist_data.append(labels[v])
            self.histogram(hist_data)


parser = argparse.ArgumentParser()
parser.add_argument('path', type=str,
                    help='Path to a read file')
parser.add_argument('-k', type=int,
                    help='Length of kmer')

args = parser.parse_args()
c = ReadCategorizer(args.path, args.k)
