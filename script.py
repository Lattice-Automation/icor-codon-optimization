import os
from Bio import SeqIO

directory = 'benchmark_seqs_aa'
for entry in os.scandir(directory):
    with open() as handle:
    for record in SeqIO.parse(entry, "fasta"):
        record.seq = record.seq()
