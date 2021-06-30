# Naively optimizes a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory.

import random
import os
from Bio import SeqIO
import re
from Bio.Seq import Seq

# Amino acid sequence dir to optimize:
aa_dir = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\aa"

# Output dir to store optimized seqs:
out_dir = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\naive"

# Amino acid to codon table, outputs arr of codons:
def aa2codons(seq : str) -> list:
    _aas = {
        "A": ["GCT GCC GCA GCG"],
        "R": ["CGT CGC CGA CGG AGA AGG"],
        "N": ["AAT AAC"],
        "D": ["GAT GAC"],
        "C": ["TGT TGC"],
        "Q": ["CAA CAG"],
        "E": ["GAA GAG"],
        "G": ["GGT GGC GGA GGG"],
        "H": ["CAT CAC"],
        "I": ["ATT ATC ATA"],
        "L": ["TTA TTG CTT CTC CTA CTG"],
        "K": ["AAA AAG"],
        "M": ["ATG ATG"],
        "F": ["TTT TTC"],
        "P": ["CCT CCC CCA CCG"],
        "S": ["TCT TCC TCA TCG AGT AGC"],
        "T": ["ACT ACC ACA ACG"],
        "W": ["TGG TGG"],
        "Y": ["TAT TAC"],
        "V": ["GTT GTC GTA GTG"],
        "B": ["GAT GAC AAT AAC"],
        "Z": ["GAA GAG CAA CAG"],
        "*": ["TAA TAG TGA"],
    }
    return [_aas[i] for i in seq]

for entry in os.scandir(aa_dir):
    name = entry.name[0:-9] + "_dna"
    record = SeqIO.read(entry,'fasta')
    arr = []

    for i in record.seq:
        arr.append(random.choice(aa2codons(i)[0][0].split()))

    record.seq = Seq(re.sub('[^GATC]',"",str("".join(arr)).upper()))
    complete_name = os.path.join(out_dir, name)
    SeqIO.write(record, complete_name, "fasta")