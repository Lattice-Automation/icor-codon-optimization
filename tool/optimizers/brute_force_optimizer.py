'''
Generates a sequence of codons and then iterates through the sequence, constantly adjusting the current codon to maximize CAI.
Goal of this is to find a combination of codons to maximize CAI (achieve 1.0 CAI).
'''

# Import modules
import os
from Bio import SeqIO
from Bio.Seq import Seq
import random
import numpy as np
import math
import re

# Set input AA sequence directory and output for writing brute sequences
aa_dir = r"..\..\benchmark_sequences\aa"
out_dir = r"..\..\benchmark_sequences\brute"

# Define weights for each codon
weights = [0,1,0.647058823500000,0.500000000000000,0.794117647100000,0.0789473684200000,0.131578947400000,0.263157894700000,0.184210526300000,0.973684210500000,1,0.851851851900000,1,1,0.587301587300000,0.818181818200000,1,0.483870967700000,0.129032258100000,1,1,0.515151515200000,0.470588235300000,1,0.384615384600000,0.307692307700000,0.871794871800000,1,1,0.754385964900000,0.180000000000000,1,0.820000000000000,0.265306122400000,0.265306122400000,1,0.0816326530600000,0.224489795900000,0.204081632700000,0.333333333300000,1,1,1,0.754385964900000,1,0.392156862700000,0.333333333300000,0.235294117600000,0.576923076900000,1,0.576923076900000,0.500000000000000,0.615384615400000,0.576923076900000,0.619047619000000,0.357142857100000,0.428571428600000,1,1,1,0.724137931000000,1,0.444444444400000,0.750000000000000,0.583333333300000]

# Create a list of all the codons and match their corresponding weights
def seq2cai(codonarray):
    output = []
    switcher = {
        'GCG': 1,
        'GCA': 2,
        'GCT': 3,
        'GCC': 4,
        'AGG': 5,
        'AGA': 6,
        'CGG': 7,
        'CGA': 8,
        'CGT': 9,
        'CGC': 10,
        'AAT': 11,
        'AAC': 12,
        'GAT': 13,
        'GAC': 14,
        'TGT': 15,
        'TGC': 16,
        'TGA': 17,
        'TAG': 18,
        'TAA': 19,
        'CAG': 20,
        'CAA': 21,
        'GAG': 22,
        'GAA': 23,
        'GGG': 24,
        'GGA': 25,
        'GGT': 26,
        'GGC': 27,
        'CAT': 28,
        'CAC': 29,
        'ATA': 30,
        'ATT': 31,
        'ATC': 32,
        'TTG': 33,
        'TTA': 34,
        'CTG': 35,
        'CTA': 36,
        'CTT': 37,
        'CTC': 38,
        'AAG': 39,
        'AAA': 40,
        'ATG': 41,
        'TTT': 42,
        'TTC': 43,
        'CCG': 44,
        'CCA': 45,
        'CCT': 46,
        'CCC': 47,
        'AGT': 48,
        'AGC': 49,
        'TCG': 50,
        'TCA': 51,
        'TCT': 52,
        'TCC': 53,
        'ACG': 54,
        'ACA': 55,
        'ACT': 56,
        'ACC': 57,
        'TGG': 58,
        'TAT': 59,
        'TAC': 60,
        'GTG': 61,
        'GTA': 62,
        'GTT': 63,
        'GTC': 64,
    }
    for codon in codonarray:
        output.append(weights[switcher.get(codon,0)])
    length = 1 / len(codonarray)
    return pow(math.prod(output), length)

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

# Converts an amino acid to a random corresponding codon:
for entry in os.scandir(aa_dir):
    # Read in the amino acid sequence:
    name = entry.replace("_aa.fasta", "_dna")
    record = SeqIO.read(entry,'fasta')
    
    masterlist = []
    bestcai = 0
    curcai = 0
    TOTAL_ITERATIONS = 100000

    for curr_iteration in range(0, TOTAL_ITERATIONS):
        codonarr = []
        # Convert amino acid to codons:
        for i in record.seq:
            #Randomly choose a codon from the list of codons for the amino acid:
            codonarr.append(random.choice(aa2codons(i)[0][0].split()))
        masterlist.append(codonarr)
        # With our new codon array, calculate the CAI:
        cai = seq2cai(codonarr)
        if (cai > curcai):
            bestcai = curr_iteration
            curcai = cai
            print('new best cai ' + str(cai))
        curr_iteration += 1
        print(curr_iteration)

    # Write the codon array to a file:
    record.seq = Seq(re.sub('[^GATC]',"",str("".join(masterlist[bestcai])).upper()))
    complete_name = os.path.join(out_dir, name)
    SeqIO.write(record, complete_name + ".fasta", "fasta")