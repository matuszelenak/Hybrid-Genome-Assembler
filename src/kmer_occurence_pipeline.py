import subprocess

from Bio import SeqIO
import argparse
import random
from matplotlib import pyplot as plt


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


def plot_histogram(occurence, k):
    plt.figure()
    plt.title('Occurence histogram for k = {}'.format(k))
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.hist(occurence, color='blue', bins=120)
    plt.show()


parser = argparse.ArgumentParser()
parser.add_argument('paths', type=str, nargs='+',
                    help='Path to a read file')
parser.add_argument('-o', '--output', help='Path of mixed read output file', type=str, default='mixed_reads')

args = parser.parse_args()

# mix_reads(args.paths, args.output)

base_command = ['./cmake-build-debug/src', '-i', '../data/ecoli/mixed', '--cov', '30']
for k in [10, 11, 12]:
    o = subprocess.check_output(base_command + ['--k', str(k)])
    occurence = [int(x) for x in o.strip().split()]
    plot_histogram(occurence, k)
