import subprocess

from Bio import SeqIO

from .utils import random_fasta_file
from .kmer_categorizer import kmer_occurence_count, kmer_positions


def test_pipeline():
    fa_path = './test_data/in.fasta'
    random_fasta_file(fa_path, num_sequences=10, seq_length=60)
    c = kmer_occurence_count(fa_path, 7, out_path='./test_data/python_count.fa')

    count_cmd = './hga count {f} --k {k} --out {out}'.format(
        k=7,
        f=fa_path,
        out='./test_data/hga_count.fa'
    )
    subprocess.run(count_cmd.split())

    python_counts = [x for x in SeqIO.parse('./test_data/python_count.fa', 'fasta')]
    python_counts.sort(key=lambda x: str(x.seq))
    hga_counts = [x for x in SeqIO.parse('./test_data/hga_count.fa', 'fasta')]
    hga_counts.sort(key=lambda x: str(x.seq))

    assert [str(x.seq) for x in python_counts] == [str(x.seq) for x in hga_counts]
    assert [x.id for x in python_counts] == [x.id for x in hga_counts]

    kmer_positions(fa_path, 7, c, 2, 4, out_path='./test_data/python_analysis')

    analysis_cmd = './hga analyze {f} --out {out} --counts {counts} --lower {lower} --upper {upper} --k {k}'.format(
        k=7,
        f=fa_path,
        out='./test_data/hga_analysis',
        counts='./test_data/hga_count.fa',
        lower=2,
        upper=4
    )
    subprocess.run(analysis_cmd.split())

    with open('./test_data/hga_analysis') as f:
        hga_data = [l.strip() for l in f.readlines()]
    with open('./test_data/python_analysis') as f:
        python_data = [l.strip() for l in f.readlines()]

    hga_data.sort(key=lambda x: x.split(':')[0])
    python_data.sort(key=lambda x: x.split(':')[0])

    assert hga_data == python_data


compile_cmd = "g++ ../kmer_categorizer.cpp ../SequenceReader.cpp ../KmerIterator.cpp ../KmerAnalysisWriter.cpp" \
              " -lboost_program_options -std=c++17 -lpthread -lstdc++fs -o hga -O3 -march=native "
compile_cmd = compile_cmd.split()
subprocess.run(compile_cmd)

test_pipeline()
