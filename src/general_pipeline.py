from Bio import SeqIO
import numpy as np
from matplotlib import pyplot as plt
import subprocess

k = 13
f = '../data/nanopore/b12.fastq'
subprocess.check_output(['jellyfish', 'count', '-m', str(k), '-s', '100M', '-t', '10', '-C', f, '-o', '../data/nanopore/{}_mer_count.jf'.format(k)])
subprocess.check_output(['jellyfish', 'histo', '../data/nanopore/{}_mer_count.jf'.format(k), '-o', '../data/nanopore/hist_values'])

with open('../data/nanopore/hist_values') as f:
    lines = [x.strip() for x in f.readlines()]
    pairs = [[int(e) for e in l.split()] for l in lines]

pairs = [p for p in pairs if 2 < p[0] < 300]

x = [p[0] for p in pairs]
y = [p[1] for p in pairs]

plt.figure()
plt.title('Occurence histogram')
plt.xlabel('Coverage')
plt.ylabel('Frequency')
plt.plot(x, y)

plt.show()
