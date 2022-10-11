'''
When testing the Genscript OptimumGene tool, they output seqs in a csv format.
This script accepts csv like so:

column 1 | column 2
file_name | seq
...

The script will covnert the CSV to sequences that will be written into an output directory which can be specified below.
'''

#import modules
import csv
from Bio import SeqIO
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# output directory to write sequences
out_dir = r"..\..\benchmark_sequences\genscript"

# iterate through the csv file and write sequences to the output directory
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
        SeqIO.write(record, os.path.join(out_dir, file_name), 'fasta')
