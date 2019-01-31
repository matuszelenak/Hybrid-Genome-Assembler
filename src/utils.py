import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet

from constants import *


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
