#!/usr/bin/env python

import argparse
import os
import subprocess
from typing import List, Tuple

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna as DNA


SIMLORD_BIN = '/home/whiskas/miniconda3/envs/simlord/bin/simlord'
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data')


def get_file_type(filename):
    _, ext = os.path.splitext(filename)
    if ext.lower() in ('.fasta', '.fa'):
        return 'fasta'
    if ext.lower() in ('.fq', '.fastq'):
        return 'fastq'
    return None


def sequence_complement(sequence: str):
    c = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get
    return ''.join(map(c, reversed(sequence)))


def indices_to_strand(indices: List[int]) -> str:
    return ''.join(DNA.letters[index] for index in indices)


def strand_to_indices(strand: str) -> List[int]:
    return [DNA.letters.index(base) for base in strand]


class ReadGeneratorBackend:
    @staticmethod
    def get_num_of_reads(seq_length: int, coverage: int, read_length: int) -> int:
        return int(coverage * seq_length / read_length)

    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True, err: float = 0) -> int:
        raise NotImplementedError


class ArtIlluminaBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True, err: float = 0) -> int:
        temporary_file = f'{record.id}_tmp.fasta'
        SeqIO.write([record], temporary_file, 'fasta')

        command = ['art_illumina', '-ss', 'HS25', '-l', '150', '-f', str(coverage), '-na', '-i', temporary_file, '-o', output_prefix]

        print(command)
        subprocess.Popen(command).wait()

        os.remove(temporary_file)

        return sum(1 for _ in SeqIO.parse(f'{output_prefix}.fq', 'fastq'))


class CustomReadGenerator(ReadGeneratorBackend):
    @staticmethod
    def apply_error(sequence: str, err: float) -> str:
        if err == 0:
            return sequence

        sequence_indices = strand_to_indices(sequence)
        substituted_index_shifts = np.random.binomial(1, err, len(sequence)) * np.random.randint(1, 4, size=len(sequence))
        substituted_sequence_indices = (sequence_indices + substituted_index_shifts) % len(DNA.letters)
        return indices_to_strand(substituted_sequence_indices)

    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True, err: float = 0) -> int:
        num_of_reads = ReadGeneratorBackend.get_num_of_reads(len(record), coverage, read_length)
        read_beginnings = np.random.randint(0, len(record) - read_length, size=num_of_reads)

        strands = [str(record.seq), sequence_complement(str(record.seq))]
        strand_choices = np.random.binomial(1, 0.5, len(read_beginnings))

        SeqIO.write(
            (
                SeqRecord(
                    Seq(CustomReadGenerator.apply_error(strands[strand_choice][beginning: beginning + read_length], err)),
                    id=f'{record.id}-{read_id}',
                    description='',
                    letter_annotations={'phred_quality': [41 for _ in range(read_length)]}
                )
                for read_id, (beginning, strand_choice) in enumerate(zip(read_beginnings, strand_choices))
            ),
            f'{output_prefix}.fq',
            'fastq'
        )

        return num_of_reads


class SimLordBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True, err: float = 0) -> int:
        temporary_file = f'{record.id}_tmp.fasta'
        SeqIO.write([record], temporary_file, 'fasta')

        cmd = [SIMLORD_BIN, '-rr', temporary_file, '--no-sam', '--without-ns', '-c', str(coverage), output_prefix]
        if read_length != -1:
            cmd.extend(['-fl', str(read_length)])

        subprocess.Popen(cmd).wait()

        os.rename(f'{output_prefix}.fastq', f'{output_prefix}.fq')

        os.remove(temporary_file)

        return sum(1 for _ in SeqIO.parse(f'{output_prefix}.fq', 'fastq'))


class NanosimHBackend(ReadGeneratorBackend):
    @staticmethod
    def reads_from_sequence(record: SeqRecord, output_prefix: str, coverage: int, read_length: int, complement: bool = True, err: float = 0) -> int:
        temporary_file = f'{record.id}_tmp.fasta'
        SeqIO.write([record], temporary_file, 'fasta')

        avg_read_length = 7777
        required_reads = round((len(record) / avg_read_length) * coverage)

        cmd = ['nanosim-h', temporary_file, '-n', str(required_reads), '--circular', '-o', output_prefix]

        subprocess.Popen(cmd).wait()

        os.remove(f'{output_prefix}.log')
        os.remove(f'{output_prefix}.errors.txt')
        os.remove(temporary_file)

        return sum(1 for _ in SeqIO.parse(f'{output_prefix}.fa', 'fasta'))


backends = (
    ('custom', CustomReadGenerator),
    ('art', ArtIlluminaBackend),
    ('simlord', SimLordBackend),
    ('nanosimh', NanosimHBackend)
)
get_backend = dict(backends).get


def generate_mutated_sequence_pair(genome_size: int, difference_rate: float) -> Tuple[SeqRecord, SeqRecord]:
    original_sequence_indices = np.random.choice(len(DNA.letters), genome_size)
    # Choose which positions to mutate and shift those positions by 1 to 3 slots
    mutation_index_shift = np.random.binomial(1, difference_rate, genome_size) * np.random.randint(1, 4, size=genome_size)

    mutated_sequence_indices = (original_sequence_indices + mutation_index_shift) % len(DNA.letters)

    original_sequence = indices_to_strand(original_sequence_indices)
    original_record = SeqRecord(
        Seq(original_sequence, DNA),
        id=f'Artificial_size_{genome_size}',
    )
    mutated_sequence = indices_to_strand(mutated_sequence_indices)
    mutated_record = SeqRecord(
        Seq(mutated_sequence, DNA),
        id=f'Artificial_mutated_size_{genome_size}'
    )

    SeqIO.write(original_record, f'../data/sequence_data/artificial.fasta', 'fasta')
    SeqIO.write(mutated_record, f'../data/sequence_data/artificial_mutated.fasta', 'fasta')

    return original_record, mutated_record


parser = argparse.ArgumentParser()
# Reference genomes parameters
parser.add_argument('-i', dest='input_files', type=str, nargs='+', help='Paths to files with reference genomes', required=False)
parser.add_argument('-n', dest='genome_size', nargs='?', type=int, help='Size of the two reference genomes')
parser.add_argument('-d', dest='difference_rate', nargs='?', type=float, help='Difference rate between the two sequences')

# Read parameters
parser.add_argument('-b', dest='backend', default='Custom', choices=[b[0] for b in backends], help='Read generation backend: Either Art Illumina or Custom')
parser.add_argument('-c', dest='coverage', nargs='?', default=30, type=int, help='Coverage in reads for one sequence')
parser.add_argument('-r', dest='read_length', nargs='?', default=-1, type=int, help='Read length. Set to -1 to keep platform-specific default')
parser.add_argument('-e', dest='error_rate', nargs='?', default=0.0, type=float, help='Error rate (for custom read backend)')

args = parser.parse_args()

if not args.input_files and not all([args.genome_size, args.difference_rate]):
    raise ValueError('Either specify input files or parameters of the genome for mutation')

if args.input_files:
    sequences = [SeqIO.read(filename, format=get_file_type(filename)) for filename in args.input_files]
else:
    original, mutated = generate_mutated_sequence_pair(args.genome_size, args.difference_rate)
    sequences = [original, mutated]

backend: ReadGeneratorBackend = get_backend(args.backend)
for seq in sequences:
    prefix = os.path.join(DATA_DIR, 'reads', f'{seq.id}_{args.backend}_{args.read_length}b_{args.coverage}x')
    read_count = backend.reads_from_sequence(
        seq,
        prefix,
        args.coverage,
        args.read_length,
        err=args.error_rate
    )
