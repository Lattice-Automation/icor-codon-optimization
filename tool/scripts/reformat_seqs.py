# reformat_seqs.py - fix FASTA formatting for each file.
# After manually importing genes from various sources, the FASTA files were inconsistent.
# This script simply opens each of the files in a directory, and re-writes them using the SeqIO function.
# This script does not change the sequence itself, but helps reformat them by writing each to lines with equal lengths.

# Import all necessary modules here
import os
from Bio import SeqIO
import random

# Change this to the directory where your files are stored.
aa_directory = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\aa"
dna_directory = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\dna"

# Print the number of files in each directory.
print("There are %d files in the amino acid directory." % len(os.listdir(aa_directory)))
print("There are %d files in the DNA directory." % len(os.listdir(dna_directory)))

# Iterate over each file in the directory.
for entry in os.scandir(aa_directory):
    record = SeqIO.read(entry, "fasta")
    SeqIO.write(record, entry, "fasta")
    #Although this does not change the actual sequence, it will reformat it with a fixed spacing (makes seqs more legible).

# Iterate over each file in the directory.
for entry in os.scandir(dna_directory):
    record = SeqIO.read(entry, "fasta")

    #Just in case, replace ambigious codons with the corresponding IUPAC ones:
    record.seq = record.seq.replace('K',random.choice(['G','T'])).replace('M',random.choice(['A','C'])).replace('N',random.choice(['A','C','G','T'])).replace('R',random.choice(['A','G'])).replace('W',random.choice(['A','T'])).replace('Y',random.choice(['C','T']))
    
    #if there are sequences that are not divisible by three, then truncate them:
    num = len(record.seq) % 3
    print("Warning: truncated" + entry.name + num)
    #warning: if sequences are being truncated, they are likely not formatted correctly.
    #all CDS should be divisible by three because they are all in frame.