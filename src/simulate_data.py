import argparse
import numpy as np
from Bio import SeqIO

bases = ['A', 'C', 'G', 'T']
base_positions = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}


def mutate_reference(sequence, difference_rate):
    mutate_indices = np.nonzero(np.random.binomial(1, difference_rate, len(sequence)))[0]
    mutations = np.random.randint(1, 4, size=len(mutate_indices))
    mutated = [x for x in sequence]
    for i, offset in zip(mutate_indices, mutations):
        mutated[i] = bases[(base_positions[mutated[i]] + offset) % 4]
    return ''.join(mutated)


def simulate_reads():
    pass


parser = argparse.ArgumentParser()
parser.add_argument('paths', type=str, nargs='+',
                    help='Path to a read file')
parser.add_argument('h', type=float, default=0.03,
                    help='Ratio of heterozygosity')
parser.add_argument('-o', '--output', help='Path of mixed read output file', type=str, default='mixed_reads')

args = parser.parse_args()

reference_sequence = ""
for path in args.paths[:1]:
    for seq_record in SeqIO.parse(path, 'fasta'):
        reference_sequence += str(seq_record.seq)

print(mutate_reference(reference_sequence, args.h))
