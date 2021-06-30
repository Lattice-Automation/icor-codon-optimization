# reformat_seqs.py - fix FASTA formatting for each file.
# After manually importing genes from various sources, the FASTA files were inconsistent.
# This script simply opens each of the files in a directory, and re-writes them using the SeqIO function.
# This script does not change the sequence itself, but helps reformat them by writing each to lines with equal lengths.

import os
from Bio import SeqIO

aa_directory = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\aa"
dna_directory = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\dna"

print("There are %d files in the amino acid directory." % len(os.listdir(aa_directory)))
print("There are %d files in the DNA directory." % len(os.listdir(dna_directory)))

for entry in os.scandir(aa_directory):
    record = SeqIO.read(entry, "fasta")
    SeqIO.write(record, entry, "fasta")

for entry in os.scandir(dna_directory):
    record = SeqIO.read(entry, "fasta")
    num = len(record.seq) % 3
    print(num)
    if (num == 0):
        SeqIO.write(record, entry, "fasta")
    else:
        record.seq = record.seq[:-num]
        SeqIO.write(record, entry, "fasta")