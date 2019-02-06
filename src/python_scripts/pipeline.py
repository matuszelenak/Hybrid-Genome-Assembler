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


def jellyfish_path(path, k):
    base_path, _ = os.path.splitext(path)
    mer_file = base_path + "_jellyfish_{kmer}mer_count.jf".format(kmer=k)
    if not os.path.exists(mer_file):
        count_cmd = 'jellyfish count -m {k} -s 100M -t 10 -C {f} -o {o}'.format(
            k=k,
            f=path,
            o=mer_file
        )
        print(count_cmd)
        subprocess.run(count_cmd.split())

    hist_file = base_path + "_jellyfish_{kmer}mer_hist.fa".format(kmer=k)
    if not os.path.exists(hist_file):
        histo_cmd = 'jellyfish dump {f} -o {o} -L 5 -U 100'.format(
            f=mer_file,
            o=hist_file
        )
        print(histo_cmd)
        subprocess.check_output(histo_cmd.split())

    histogram(hist_file, 'jellyfish')

    compile_cmd = "g++ ../kmer_categorizer.cpp ../SequenceReader.cpp ../KmerIterator.cpp" \
                  " -lboost_program_options -std=c++17 -lstdc++fs -o hga -O3 -march=native"
    compile_cmd = compile_cmd.split()
    print(compile_cmd)
    subprocess.run(compile_cmd)

    analysis_file = base_path + "_hga_{kmer}mer_analysis".format(kmer=k)
    if not os.path.exists(analysis_file):
        print("Enter the lower and upper bound for kmer filtering")
        lower, upper = [int(x) for x in input().split()]
        analysis_cmd = './hga analyze {f} --out {out} --counts {counts} --lower {lower} --upper {upper} --k {k}'.format(
            f=path,
            out=analysis_file,
            counts=hist_file,
            lower=lower,
            upper=upper,
            k=k
        )
        print(analysis_cmd)
        subprocess.run(analysis_cmd.split())


def hga_path(path, k):
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
parser.add_argument('-tool', type=str, default='jellyfish', help='Tool for k-mer counting')

args = parser.parse_args()

if not os.path.exists(args.path):
    raise FileNotFoundError("{} does not exist".format(args.path))

jellyfish_path(args.path, args.k)
