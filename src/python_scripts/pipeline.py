from Bio import SeqIO
from matplotlib import pyplot as plt
import subprocess
import argparse
import os


def histogram(path, tool="jellyfish"):
    data = []
    for seq_record in SeqIO.parse(path, 'fasta'):
        data.append(int(seq_record.id))

    plt.figure()
    plt.title('Occurence histogram ({})'.format(tool))
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.hist(data, color='blue', bins=120)
    plt.show()


def jellyfish_histogram(path, k):
    base_path, _ = os.path.splitext(path)
    mer_file = base_path + "_jellyfish_{kmer}mer_count.jf".format(kmer=k)
    count_cmd = 'jellyfish count -m {k} -s 100M -t 10 -C {f} -o {o}'.format(
        k=k,
        f=path,
        o=mer_file
    )
    print(count_cmd)
    subprocess.run(count_cmd.split())
    hist_file = base_path + "_jellyfish_{kmer}mer_hist.fa".format(kmer=k)
    histo_cmd = 'jellyfish dump {f} -o {o}'.format(
        f=mer_file,
        o=hist_file
    )
    print(histo_cmd)
    subprocess.check_output(histo_cmd.split())

    histogram(hist_file, 'jellyfish')


def hga_histogram(path, k):
    base_path, _ = os.path.splitext(path)
    compile_cmd = "g++ ../kmer_categorizer.cpp ../SequenceReader.cpp ../KmerIterator.cpp" \
                  " -lboost_program_options -std=c++17 -lstdc++fs -o hga"
    compile_cmd = compile_cmd.split()
    subprocess.run(compile_cmd)

    hist_file = base_path + "_hga_{kmer}mer_hist.fa".format(kmer=k)
    count_cmd = './hga count {f} --k {k} --out {out}'.format(
        k=k,
        f=path,
        out=hist_file
    )
    print(count_cmd)
    subprocess.run(count_cmd.split())

    histogram(hist_file, 'HGA')


parser = argparse.ArgumentParser()
parser.add_argument('path', type=str,
                    help='Path to a read file')
parser.add_argument('-k', type=int, help='Length of k-mer')

args = parser.parse_args()

if not os.path.exists(args.path):
    raise FileNotFoundError("{} does not exist".format(args.path))

jellyfish_histogram(args.path, args.k)
hga_histogram(args.path, args.k)
