'''
When testing the Genscript OptimumGene tool, they output seqs in a csv format.
This script accepts csv like so:

column 1 | column 2
file_name | seq
...

'''

import csv
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

out_dir = r"C:\Users\risha\Desktop\icor-codon-optimization\benchmark_sequences\genscript"

with open('optimum_seqs.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))
    for i in data:
        file_name = i[0]
        record = SeqRecord(
            Seq(i[1]),
            id=file_name,
            name="file_name",
            description="blank",
        )
        print(record)
        SeqIO.write(record,os.path.join(out_dir, file_name),'fasta')