from matplotlib import pyplot as plt

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from constants import *


def histogram(path, tool="jellyfish"):
    data = []
    for seq_record in SeqIO.parse(path, 'fasta'):
        data.append(int(seq_record.id))

    data = [x for x in data if 1 < x < 100]

    plt.figure()
    plt.title('Occurence histogram ({})'.format(tool))
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.hist(data, color='blue', bins=120)
    plt.show()


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
