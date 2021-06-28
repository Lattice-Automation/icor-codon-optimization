# reformat_seqs.py - fix FASTA formatting for each file.
# After manually importing genes from various sources, the FASTA files were inconsistent.
# This script simply opens each of the files in a directory, and re-writes them using the SeqIO function.
# This script does not change the sequence itself, but helps reformat them by writing each to lines with equal lengths.

import os
from Bio import SeqIO

#directory = 'benchmark_seqs/dna'
directory = 'benchmark_seqs/aa'

numFiles = len(os.listdir(directory))
print("There are %d files in this directory." % numFiles)

for entry in os.scandir(directory):
    record = SeqIO.read(entry, "fasta")
    record.sequence = record.seq
    SeqIO.write(record, entry, "fasta")