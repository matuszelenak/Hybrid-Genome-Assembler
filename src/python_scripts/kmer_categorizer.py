import argparse
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import compress_kmer
from constants import *


def kmer_occurence_count(path, k, out_path=None):
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

    out_path = out_path or (base_path + '_python_counts.fa')

    with open(out_path, 'w') as f:
        ordered = [SeqRecord(Seq(sequence), id=str(count), description='') for sequence, count in kmer_counts.items()]
        ordered.sort(key=lambda x: str(x.seq))
        SeqIO.write(ordered, f, 'fasta')

    return kmer_counts


def kmer_positions(path, k, kmers, lower, upper, out_path=None):
    kmers = set([sequence for sequence, count in kmers.items() if lower <= count <= upper])

    base_path, extension = os.path.splitext(path)
    if extension in ('.fasta', '.fa'):
        extension = 'fasta'
    else:
        extension = 'fastq'

    out_path = out_path or (base_path + '_python_analysis')
    with open(out_path, 'w') as f:
        for seq_record in SeqIO.parse(path, extension):
            sequence = str(seq_record.seq)
            f.write(seq_record.id + ':[')
            for i in range(k, len(sequence) + 1):
                forward_kmer = sequence[i - k: i]
                backward_kmer = ''.join([complement[x] for x in forward_kmer])[::-1]
                signature = min(forward_kmer, backward_kmer)
                if signature in kmers:
                    f.write('({},{}),'.format(compress_kmer(signature), i - k))
            f.write(']\n')

    return out_path
