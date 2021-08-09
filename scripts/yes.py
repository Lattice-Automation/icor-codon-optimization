
# Import modules
import os
from Bio import SeqIO
from Bio.Seq import Seq
import random
import numpy as np
import math
import re

# Set input AA sequence directory and output for writing brute sequences
aa_dir = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\aa"

for entry in os.scandir(aa_dir):
    record = SeqIO.read(entry,"fasta")
    print(len(record.seq))