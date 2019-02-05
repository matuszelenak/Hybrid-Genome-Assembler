from .utils import *

import random
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('paths', type=str, nargs='+',
                    help='Path to a read file')
parser.add_argument('-d', type=float, default=0.03,
                    help='Ratio of heterozygosity')
parser.add_argument('-o', '--output', help='Path of mixed read output file', type=str, default='mixed_reads')

args = parser.parse_args()


reference_sequence = ""

if args.paths == ['random']:
    reference_sequence = ''.join(np.random.choice(['A', 'C', 'G', 'T'], 30000))
else:
    for path in args.paths[:1]:
        reference_sequence = ''.join([x for x in read_sequence_from_file(path)])
    reference_sequence = ''.join([x for x in reference_sequence if x in ('A', 'C', 'G', 'T')])

mutated = mutate_sequence(reference_sequence, args.d)

write_sequences_to_file([reference_sequence], args.output + '_reference', data_format='fasta')
write_sequences_to_file([mutated], args.output + '_mutated', data_format='fasta')
reference_reads = simulate_perfect_reads(reference_sequence, read_length=200)
mutated_reads = simulate_perfect_reads(mutated, read_length=200)

write_sequences_to_file(reference_reads, args.output + '_reads_reference')
write_sequences_to_file(mutated_reads, args.output + '_reads_mutated')

mixed_reads = reference_reads + mutated_reads
random.shuffle(mixed_reads)
write_sequences_to_file(mixed_reads, args.output, data_format='fasta')
