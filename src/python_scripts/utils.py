import random

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet

from constants import *


def mix_reads(paths, output_path):
    raw_read_sequences = []
    for path in paths:
        for seq_record in SeqIO.parse(path, 'fastq'):
            raw_read_sequences.append(str(seq_record.seq))

    random.shuffle(raw_read_sequences)

    f = open(output_path, 'w')
    for seq in raw_read_sequences:
        f.write('>read_descriptor\n')
        f.write(seq)
        f.write('\n')

    f.close()


def mutate_sequence(sequence, difference_rate):
    mutate_indices = np.nonzero(np.random.binomial(1, difference_rate, len(sequence)))[0]
    mutations = np.random.randint(1, 4, size=len(mutate_indices))
    mutated = [x for x in sequence]
    for i, offset in zip(mutate_indices, mutations):
        mutated[i] = bases[(base_positions[mutated[i]] + offset) % 4]
    return ''.join(mutated)


def simulate_perfect_reads(sequence, read_length=200, coverage=30):
    num_of_reads = int(coverage * len(sequence) / read_length)
    read_beginnings = np.random.randint(0, len(sequence) - read_length, size=num_of_reads)
    complementary = ''.join([complement[x] for x in sequence[::-1]])
    strand_choices = np.random.binomial(1, 0.5, num_of_reads)
    strands = [sequence, complementary]

    reads = []
    for strand_choice, beginning in zip(strand_choices, read_beginnings):
        reads.append(strands[strand_choice][beginning:beginning + read_length])

    return reads


def read_sequence_from_file(path, data_format='fasta'):
    for seq_record in SeqIO.parse(path, data_format):
        yield str(seq_record.seq)


def write_sequences_to_file(sequences, path, data_format='fasta'):
    if data_format in ('fasta', 'plain'):
        with open(path, 'w') as f:
            for i, s in enumerate(sequences):
                if data_format == 'fasta':
                    f.write('>sequence_{}\n'.format(i))
                f.write(s + '\n')
        return

    records = []
    for i, s in enumerate(sequences):
        records.append(SeqRecord(Seq(s, DNAAlphabet), id='sequence_{}'.format(i)))

    with open(path, 'w') as f:
        SeqIO.write(records, f, data_format)


def random_fasta_file(path, num_sequences=10, seq_length=100):
    with open(path, 'w') as f:
        records = []
        for i in range(num_sequences):
            records.append(
                SeqRecord(
                    Seq(''.join(np.random.choice(['A', 'C', 'G', 'T'], seq_length))),
                    id='sequence_{}'.format(i),
                    description=''
                )
            )
        SeqIO.write(records, f, 'fasta')


def compress_kmer(sequence):
    result = 0
    for c in sequence:
        result <<= 2
        result |= base_positions[c]
    return result
