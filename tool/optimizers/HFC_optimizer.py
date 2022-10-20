# import standard modules; re is used for regex on ln 59
import os
from Bio import SeqIO
from Bio.Seq import Seq
import random
import numpy as np
import re

# frequencies are from ecoli_codon_frequencies.xlsx file in summaries dir
# create dict with value being a tuple with the codons and their probabilities/frequncies
frequency = {
    "A": "GCG",
    "R": "CGT",
    "N": "AAC",
    "D": "GAT",
    "C": "TGC",
    "*": "TAA",
    "Q": "CAG",
    "E": "GAA",
    "G": "GGC",
    "H": "CAT",
    "I": "ATT",
    "L": "CTG",
    "K": "AAA",
    "M": "ATG",
    "F": "TTT",
    "P": "CCG",
    "S": "AGC",
    "T": "ACC",
    "W": "TGG",
    "Y": "TAT",
    "V": "GTG"
}

# Amino acid sequence dir to optimize:
# hardcoded path
aa_dir = os.path.join(os.getcwd(), 'benchmark_sequences', 'aa')

# Output dir to store optimized seqs:
# hardcoded path
out_dir = os.path.join(os.getcwd(), 'benchmark_sequences', 'HFC')

for entry in os.scandir(aa_dir):
    name = entry.name.replace("_aa.fasta", "_dna")

    # Replace ambiguities with amino acids from IUPAC guidelines: https://www.bioinformatics.org/sms/iupac.html
    record = SeqIO.read(entry, "fasta")
    seq = record.seq.replace("B", random.choice(["D", "N"])).replace(
        "Z", random.choice(["E", "Q"]))
    seq_arr = []
    for aa in seq:
        # append to the array a random choice of codon using the probabilities given (p)
        seq_arr.append(frequency[aa])

    record.seq = Seq(re.sub('[^GATC]', "", str("".join(seq_arr)).upper()))
    complete_name = os.path.join(out_dir, name)
    SeqIO.write(record, complete_name + ".fasta", "fasta")