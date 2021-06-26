# ICOR Codon Optimization Tool
---
- [Benchmarking](#Benchmarking)
    - [Sequences](#Benchmark%20Sequences)
---
## Benchmarking

### Benchmark Sequences
`benchmark_seqs_dna` contains `x` DNA sequences in the FASTA format for benchmarking. `benchmark_seqs_aa` contains `x` AA sequences that correspond to the DNA sequence folder. Below is a table depicting the importance of each file.

Descriptions marked with "**" indicate genes that have been expressed in *E. coli* in past studies. The codon optimization tool has been trained on *E. coli genes*. Therefore, it is best applied to improve expression for genes that would be expressed by using *E. coli* as a cell factory. *E. coli* is the ideal cell host for recombinant expression because of its efficiency and low-cost.

|        File         | Description |
|       :---:         | ----------- |
| FALVAC-1  | ** FALVAC-1 is a vaccine against Plasmodium Falciparum. This is a good benchmark for a synthetic gene as it is a real candidate vaccine that is recombinantly produced ([Reference Paper](https://doi.org/10.1016%2Fj.pep.2003.11.006)).|
| PEA   | ** *Pseudomonas aeruginosa* exotoxin A (PEA) is an important pathogenic factor [1](https://doi.org/10.1073/pnas.85.9.2939). It retains high immunogenicity even after detoxification, enabling its use as vaccine adjuvants and vaccine carriers. This is a good benchmark as it can be produced for a low-cost at a large-scale in *E. coli*. Further, the ([Reference Paper](https://www.sciencedirect.com/science/article/pii/S1046592810000501?via%3Dihub#bib2)) finds that codon optimization enhances expression of PEA in *E. coli* -- thus, if the tool presented in this research can achieve similar/better results than the paper and/or other approaches, it will be considered an improvement. |
| msox   | ** Monomeric sarcosine oxidase (Msox) is an important diagnostic enzyme and can catalyze the oxidation of N-methyl with the FAD cofactor. ([Study](https://www.sciencedirect.com/science/article/pii/S0168165615301905?via%3Dihub)) found that the expression of recombinant SOX "has been extensively studied." Thus, the sox gene serves as a good benchmark because the extent to which codon optimization affects it has been studied. |
| soxB   | ** Sarcosine oxidase subunit B (soxB) is an important diagnostic enzyme and can catalyze the oxidation of N-methyl with the FAD cofactor. ([Study](https://www.sciencedirect.com/science/article/pii/S0168165615301905?via%3Dihub)) found that the expression of recombinant SOX "has been extensively studied." Thus, the sox gene serves as a good benchmark because the extent to which codon optimization affects it has been studied. |
| PRL-3/PTP4A3  | ** Protein tyrosine phosphatase 4A3 (PTP4A3/PRL-3) is a well-documented protein that can be produced recombinantly in *E. coli.* This protein has been identified as a potential target to treat some cancers. |
| hPDF  | ** Human peptide deformylase (hPDF) is a target for cancer therapeutics. However, its expression is not very efficient in *E. coli.* This serves as a good target for benchmarks as past [studies](https://pubmed.ncbi.nlm.nih.gov/19825416/) have noted the valuable potential of codon optimization of this gene. |
| PF3D7  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| PA  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MmpL3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MmpL3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MmpL3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MmpL3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |