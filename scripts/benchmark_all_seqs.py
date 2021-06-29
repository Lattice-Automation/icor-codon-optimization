import os
from Bio import SeqIO

aa_seqs = 'benchmark_seqs/aa'

numFiles = len(os.listdir(aa_seqs))

for i in range(numFiles):
    