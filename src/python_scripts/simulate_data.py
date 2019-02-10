import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os
import random
import argparse


bases = ['A', 'C', 'G', 'T']
base_positions = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3
}

complement = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


def read_sequence_from_file(path):
    _, extension = os.path.splitext(path)
    if extension in ('.fa', 'fasta'):
        data_format = 'fasta'
    else:
        data_format = 'fastq'
    for seq in SeqIO.parse(path, data_format):
        seq.description = ''
        yield seq


def write_sequences_to_file(sequences, path, data_format='fasta'):
    if data_format in ('fasta', 'plain'):
        with open(path, 'w') as f:
            for i, s in enumerate(sequences):
                if data_format == 'fasta':
                    f.write(">" + str(s.id) + '\n')
                f.write(str(s.seq) + '\n')
        return

    with open(path, 'w') as f:
        SeqIO.write(sequences, f, data_format)


def mutate_sequence(sequence_record, difference_rate):
    sequence = str(sequence_record.seq)
    mutate_indices = np.nonzero(np.random.binomial(1, difference_rate, len(sequence)))[0]
    mutations = np.random.randint(1, 4, size=len(mutate_indices))
    mutated = [x for x in sequence]
    for i, offset in zip(mutate_indices, mutations):
        mutated[i] = bases[(base_positions[mutated[i]] + offset) % 4]
    return SeqRecord(
        seq=''.join(mutated),
        id=sequence_record.id + '|mutated_diff={}'.format(difference_rate),
        description=''
    )


def simulate_perfect_reads(sequence_record, read_length=200, coverage=30, circular=False):
    sequence = str(sequence_record.seq)
    num_of_reads = int(coverage * len(sequence) / read_length)
    if circular:
        read_beginnings = np.random.randint(0, len(sequence), size=num_of_reads)
    else:
        read_beginnings = np.random.randint(0, len(sequence) - read_length, size=num_of_reads)
    complementary = ''.join([complement[x] for x in sequence[::-1]])
    strand_choices = np.random.binomial(1, 0.5, num_of_reads)
    strands = [sequence, complementary]

    reads = []
    for i, strand_choice, beginning in zip(range(num_of_reads), strand_choices, read_beginnings):
        seq = strands[strand_choice][beginning:beginning + read_length]
        if beginning + read_length > len(strands[strand_choice]):
            seq += strands[strand_choice][:(beginning + read_length) % len(strands[strand_choice])]
        reads.append(
            SeqRecord(
                seq=seq,
                # id=sequence_record.id + '|read={}_length={}'.format(i, read_length),
                id=sequence_record.id,
                description=''
            )
        )

    print('Simulated {} reads'.format(num_of_reads))
    return reads


parser = argparse.ArgumentParser()
parser.add_argument('paths', type=str, nargs='*',
                    help='Path to a read file')
parser.add_argument('-d', type=float, default=0.03,
                    help='Ratio of heterozygosity from reference sequence')
parser.add_argument('-l', type=int, default=300000,
                    help='Length of generated sequence (if no read file supplied)')
parser.add_argument('-r', type=int, default=200,
                    help='Length of simulated reads')
parser.add_argument('-c', type=int, default=15,
                    help='Coverage of reads')
parser.add_argument('-C', action='store_true',
                    help='Circular genome')
parser.add_argument('-o', '--output', help='Path of mixed read output file', type=str, default='mixed_reads.fasta')

args = parser.parse_args()


loaded_sequences = []
if not args.paths:
    loaded_sequences.append(
        SeqRecord(
            seq=Seq(''.join(np.random.choice(['A', 'C', 'G', 'T'], args.l))),
            id='random_sequence_length={}'.format(args.l),
            description=''
        )
    )
else:
    for path in args.paths:
        loaded_sequences.append(read_sequence_from_file(path))

if len(loaded_sequences) == 1:
    loaded_sequences.append(mutate_sequence(loaded_sequences[0], args.d))

mixed_reads = []
for seq_record in loaded_sequences:
    mixed_reads += simulate_perfect_reads(seq_record, read_length=args.r, coverage=args.c, circular=args.C)

random.shuffle(mixed_reads)
write_sequences_to_file(mixed_reads, args.output, data_format='fasta')
