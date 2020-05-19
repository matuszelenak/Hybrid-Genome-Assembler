# Hybrid-Genome-Assembler
FMFI UK Master's degree thesis code and sample data

To compile the project in the src folder you need
- cmake > 3.0
- c++ compiler with C++20 support
- libeigen2-dev
- libfmt-dev
- robin-map-dev
- Boost libraries

Run in the src folder
```
cmake .
make	
```

Python dependencies:
Requirements:
- numpy
- biopython, Nanosim-H for read generation
- matplotlib, networkx for plotting

## Generating reads

### Generating from a random sequence

Run ```python read_generator.py -n 500000 -d 0.03 -b art -c 30``` to generate Art Illumina reads of a random 500000 bases long sequence and its 3% mutated pair with 30x coverage. See program help for more options

### Generating from reference files

Run ```python read_generator.py -i sequence_1 sequence_2 ... sequence_n -b nanosimh -c 50``` to generate reads using Nanosim-H from provided sequences in FASTA format with 50x coverage

## Discriminative kmer computation

Run ```jf_occurrences read_file_1 ... read_file_n -k 19``` to compute a histogram of 19-mers from read files. Omit the -k argument to have the program determine the k value itself. Every file is treated as reads from a separate haplotype when calculating the histogram.

After plotting of histogram, you are queried with entering the range of occurrences from which to export kmers, as well as a sampling rate (in range 0 to 1).

## Running read categorization

Run ```./categorization read_file_1 read_file_2 ... read_file_n --kmers file_with_kmers``` to categorize reads using a discriminative kmer set in `file_with_kmers`. See help for more options.

### Replication of pipeline used in the thesis:

Run the following commands

```
python scripts/read_generator.py -i ../data/sequence/ecoli-mg1655.fa ../data/sequences/ecoli-UTI89.fa -b art -c 30
python scripts/read_generator.py -i ../data/sequence/ecoli-mg1655.fa ../data/sequences/ecoli-UTI89.fa -b nanosimh -c 75
./jf_occurrences ecoli_MG1655_art_30x.fq ecoli_UTI89_art_30x.fq -k 19
```
Select a range of 10 to 25, with sampling rate 1

Run 
```
./categorization ecoli_MG1655_nanosimh_75x.fq ecoli_UTI89_nanosimh_75x.fq --kmers 19-mers_10_25_100%.txt -o components
```
In the components folder you will find the reads assigned in components