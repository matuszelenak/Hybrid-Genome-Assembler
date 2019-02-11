import argparse
import os
import random

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('paths', type=str, nargs='*',
                    help='Path to a read file')
parser.add_argument('-o', '--output', help='Path of mixed read output file', type=str, default='mixed_reads.fasta')

args = parser.parse_args()

reads = []
for path in args.paths:
    _, extension = os.path.splitext(path)
    if extension in ('.fa', 'fasta'):
        data_format = 'fasta'
    else:
        data_format = 'fastq'
    for seq in SeqIO.parse(path, data_format):
        seq.id = os.path.splitext(os.path.split(path)[-1])[0]
        seq.description = ''
        reads.append(seq)

random.shuffle(reads)
with open(args.output, 'w') as f:
    for s in reads:
        f.write(">" + str(s.id) + '\n')
        f.write(str(s.seq) + '\n')
