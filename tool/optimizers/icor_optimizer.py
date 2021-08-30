# Define variables (must change!)
model_path = r"C:\Users\risha\Desktop\icor-codon-optimization\tool\models\icor.onnx"

# Import packages
from Bio.Seq import Seq
from Bio.Seq import translate
import sys
import onnxruntime as rt
import numpy as np

type = input("Welcome to ICOR! Are you optimizing an amino acid sequence (enter in 'aa' below) or a dna/codon sequence (enter in 'dna' below)?\n\n").strip().upper()
input_seq = input(
    "Enter the coding sequence only.\nEnter in 'demo' to use demo sequence.\n\n").strip().upper()

# Load demo sequence (AKT1 amino acid seq)
if type == 'AA':
    if input_seq == 'DEMO':
        input_seq = "MSDVAIVKEGWLHKRGEYIKTWRPRYFLLKNDGTFIGYKERPQDVDQREAPLNNFSVAQCQLMKTERPRPNTFIIRCLQWTTVIERTFHVETPEEREEWTTAIQTVADGLKKQEEEEMDFRSGSPSDNSGAEEMEVSLAKPKHRVTMNEFEYLKLLGKGTFGKVILVKEKATGRYYAMKILKKEVIVAKDEVAHTLTENRVLQNSRHPFLTALKYSFQTHDRLCFVMEYANGGELFFHLSRERVFSEDRARFYGAEIVSALDYLHSEKNVVYRDLKLENLMLDKDGHIKITDFGLCKEGIKDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHEKLFELILMEEIRFPRTLGPEAKSLLSGLLKKDPKQRLGGGSEDAKEIMQHRFFAGIVWQHVYEKKLSPPFKPQVTSETDTRYFDEEFTAQMITITPPDQDDSMECVDSERRPHFPQFSYSASGTA*"
    if not input_seq.startswith('M') or not input_seq.endswith('*'):
        sys.exit('Invalid amino acid sequence detected.\nThe sequence must start with M and end with * because ICOR only optimizes the codon-sequence region!\nPlease try again.\nRead more: http://www.hgvs.org/mutnomen/references.html#aalist')
elif type == 'DNA':
    if input_seq == 'DEMO':
        input_seq = "ATGAGCGACGTGGCTATTGTGAAGGAGGGTTGGCTGCACAAACGAGGGGAGTACATCAAGACCTGGCGGCCACGCTACTTCCTCCTCAAGAATGATGGCACCTTCATTGGCTACAAGGAGCGGCCGCAGGATGTGGACCAACGTGAGGCTCCCCTCAACAACTTCTCTGTGGCGCAGTGCCAGCTGATGAAGACGGAGCGGCCCCGGCCCAACACCTTCATCATCCGCTGCCTGCAGTGGACCACTGTCATCGAACGCACCTTCCATGTGGAGACTCCTGAGGAGCGGGAGGAGTGGACAACCGCCATCCAGACTGTGGCTGACGGCCTCAAGAAGCAGGAGGAGGAGGAGATGGACTTCCGGTCGGGCTCACCCAGTGACAACTCAGGGGCTGAAGAGATGGAGGTGTCCCTGGCCAAGCCCAAGCACCGCGTGACCATGAACGAGTTTGAGTACCTGAAGCTGCTGGGCAAGGGCACTTTCGGCAAGGTGATCCTGGTGAAGGAGAAGGCCACAGGCCGCTACTACGCCATGAAGATCCTCAAGAAGGAAGTCATCGTGGCCAAGGACGAGGTGGCCCACACACTCACCGAGAACCGCGTCCTGCAGAACTCCAGGCACCCCTTCCTCACAGCCCTGAAGTACTCTTTCCAGACCCACGACCGCCTCTGCTTTGTCATGGAGTACGCCAACGGGGGCGAGCTGTTCTTCCACCTGTCCCGGGAGCGTGTGTTCTCCGAGGACCGGGCCCGCTTCTATGGCGCTGAGATTGTGTCAGCCCTGGACTACCTGCACTCGGAGAAGAACGTGGTGTACCGGGACCTCAAGCTGGAGAACCTCATGCTGGACAAGGACGGGCACATTAAGATCACAGACTTCGGGCTGTGCAAGGAGGGGATCAAGGACGGTGCCACCATGAAGACCTTTTGCGGCACACCTGAGTACCTGGCCCCCGAGGTGCTGGAGGACAATGACTACGGCCGTGCAGTGGACTGGTGGGGGCTGGGCGTGGTCATGTACGAGATGATGTGCGGTCGCCTGCCCTTCTACAACCAGGACCATGAGAAGCTTTTTGAGCTCATCCTCATGGAGGAGATCCGCTTCCCGCGCACGCTTGGTCCCGAGGCCAAGTCCTTGCTTTCAGGGCTGCTCAAGAAGGACCCCAAGCAGAGGCTTGGCGGGGGCTCCGAGGACGCCAAGGAGATCATGCAGCATCGCTTCTTTGCCGGTATCGTGTGGCAGCACGTGTACGAGAAGAAGCTCAGCCCACCCTTCAAGCCCCAGGTCACGTCGGAGACTGACACCAGGTATTTTGATGAGGAGTTCACGGCCCAGATGATCACCATCACACCACCTGACCAAGATGACAGCATGGAGTGTGTGGACAGCGAGCGCAGGCCCCACTTCCCCCAGTTCTCCTACTCGGCCAGCGGCACGGCCTGA"
    if 'U' in input_seq:
        sys.exit('Invalid DNA sequence detected.\nThe sequence must be in DNA form. A "U" was found in your sequence.\nPlease try again.\nRead more: http://www.hgvs.org/mutnomen/references.html#aalist')
    if not input_seq.startswith('ATG'):
        sys.exit('Invalid DNA sequence detected.\nThe sequence must start with ATG because ICOR only optimizes the codon-sequence region! Please try again.\nRead more: http://www.hgvs.org/mutnomen/references.html#aalist')
    if not input_seq.endswith('TAA') and not input_seq.endswith('TGA') and not input_seq.endswith('TAG'):
        sys.exit('Invalid DNA sequence detected.\nThe sequence must end with a stop codon of TAA, TGA, or TAG because ICOR only optimizes the codon-sequence region! Please try again.\nRead more: http://www.hgvs.org/mutnomen/references.html#aalist')
    
    # ICOR accepts the amino acid sequence, so we translate the DNA sequence to amino acid sequence:
    input_seq = Seq(input_seq)
    input_seq = input_seq.translate()

print(input_seq)
# Define categorical labels from when model was trained.
labels = ['AAA', 'AAC','AAG','AAT','ACA','ACG','ACT','AGC','ATA','ATC','ATG','ATT','CAA','CAC','CAG','CCG','CCT','CTA','CTC','CTG','CTT','GAA','GAT','GCA','GCC','GCG','GCT','GGA','GGC','GTC','GTG','GTT','TAA','TAT','TCA','TCG','TCT','TGG','TGT','TTA','TTC','TTG','TTT','ACC','CAT','CCA','CGG','CGT','GAC','GAG','GGT','AGT','GGG','GTA','TGC','CCC','CGA','CGC','TAC','TAG','TCC','AGA','AGG','TGA']

# Define aa to integer table
def aa2int(seq: str) -> list:
    _aa2int = {
        'A': 1,
        'R': 2,
        'N': 3,
        'D': 4,
        'C': 5,
        'Q': 6,
        'E': 7,
        'G': 8,
        'H': 9,
        'I': 10,
        'L': 11,
        'K': 12,
        'M': 13,
        'F': 14,
        'P': 15,
        'S': 16,
        'T': 17,
        'W': 18,
        'Y': 19,
        'V': 20,
        'B': 21,
        'Z': 22,
        'X': 23,
        '*': 24,
        '-': 25,
        '?': 26
    }
    return [_aa2int[i] for i in seq]


# Create empty array to fill
oh_array = np.zeros(shape=(26, len(input_seq)))

# Load placements from aa2int
aa_placement = aa2int(input_seq)

# One-hot encode the amino acid sequence:
i = 0
while i < len(aa_placement):
    oh_array[aa_placement[i], i] = 1
    i += 1

oh_array = [oh_array]
x = np.array(np.transpose(oh_array))

y = x.astype(np.float32)

y = np.reshape(y, (y.shape[0], 1, 26))

# Start ICOR session using model.
sess = rt.InferenceSession(model_path)
input_name = sess.get_inputs()[0].name

# Get prediction:
pred_onx = sess.run(None, {input_name: y})

# Get the index of the highest probability from softmax output:
pred_indices = []
for pred in pred_onx[0]:
    pred_indices.append(np.argmax(pred))

str = ""
for index in pred_indices:
    str += labels[index]
print('==== OUTPUT ====\n' + str)

output = input(
    "Would you like to write this into a file? (Y or N)\n\n").strip().upper()

if (output == 'Y'):
    with open('output.txt', 'w') as f:
        f.write(str)
    print('\nOutput written to output.txt')
else:
    print('\nNo output written. Done!')
