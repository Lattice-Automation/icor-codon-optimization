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
    "A": (["GCG", "GCA", "GCT", "GCC"], [0.34, 0.22, 0.17, 0.27]),
    "R": (["AGG", "AGA", "CGG", "CGA", "CGT", "CGC"], [0.03, 0.05, 0.1, 0.07, 0.37, 0.38]),
    "N": (["AAT", "AAC"], [0.46, 0.54]),
    "D": (["GAT", "GAC"], [0.63, 0.37]),
    "C": (["TGT", "TGC"], [0.45, 0.55]),
    "*": (["TGA", "TAG", "TAA"], [0.3, 0.08, 0.62]),
    "Q": (["CAG", "CAA"], [0.66, 0.34]),
    "E": (["GAG", "GAA"], [0.32, 0.68]),
    "G": (["GGG", "GGA", "GGT", "GGC"], [0.15, 0.12, 0.34, 0.39]),
    "H": (["CAT", "CAC"], [0.57, 0.43]),
    "I": (["ATA", "ATT", "ATC"], [0.09, 0.5, 0.41]),
    "L": (["TTG", "TTA", "CTG", "CTA", "CTT", "CTC"], [0.13, 0.13, 0.49, 0.04, 0.11, 0.1]),
    "K": (["AAG", "AAG"], [0.25, 0.75]),
    "M": (["ATG"], [1.0]),
    "F": (["TTT", "TTC"], [0.57, 0.43]),
    "P": (["CCG", "CCA", "CCT", "CCC"], [0.51, 0.2, 0.17, 0.12]),
    "S": (["AGT", "AGC", "TCG", "TCA", "TCT", "TCC"], [0.15, 0.26, 0.15, 0.13, 0.16, 0.15]),
    "T": (["ACG", "ACA", "ACT", "ACC"], [0.26, 0.15, 0.18, 0.42]),
    "W": (["TGG"], [1.0]),
    "Y": (["TAT", "TAC"], [0.58, 0.42]),
    "V": (["GTG", "GTA", "GTT", "GTC"], [0.36, 0.16, 0.27, 0.21])
}

# Amino acid sequence dir to optimize:
# hardcoded path
aa_dir = os.path.join(os.getcwd(), 'benchmark_sequences', 'aa')

# Output dir to store optimized seqs:
# hardcoded path
out_dir = os.path.join(os.getcwd(), 'benchmark_sequences', 'naive')


# Normalize probabilities for frequency if sum is not exactly 1.
def fix_p(p):
    if p.sum() != 1.0:
        p = p*(1./p.sum())
    return p


for entry in os.scandir(aa_dir):
    name = entry.replace("_aa.fasta", "_dna")

    # Replace ambiguities with amino acids from IUPAC guidelines: https://www.bioinformatics.org/sms/iupac.html
    record = SeqIO.read(entry, "fasta")
    seq = record.seq.replace("B", random.choice(["D", "N"])).replace(
        "Z", random.choice(["E", "Q"]))
    seq_arr = []
    for aa in seq:
        # append to the array a random choice of codon using the probabilities given (p)
        seq_arr.append(np.random.choice(
            frequency[aa][0], p=fix_p(np.asarray(frequency[aa][1]))))

    record.seq = Seq(re.sub('[^GATC]', "", str("".join(seq_arr)).upper()))
    complete_name = os.path.join(out_dir, name)
    SeqIO.write(record, complete_name + ".fasta", "fasta")
