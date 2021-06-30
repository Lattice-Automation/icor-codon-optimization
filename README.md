# ICOR Codon Optimization Tool
---
- [Benchmarking](#Benchmarking)
- [Scripts](#Scripts)
---
## Benchmarking

### Benchmark Sequences
`benchmark_sequences` is a folder that contains sequences for benchmarking purposes:
- `aa` which consists of 40 amino acid sequences in the FASTA format.
- `dna` which consists of 40 DNA sequences in the FASTA format.

Descriptions marked with "**" indicate genes that have been expressed in *E. coli* in past studies. The codon optimization tool has been trained on *E. coli genes*. Therefore, it is best applied to improve expression for genes that would be expressed by using *E. coli* as a cell factory. *E. coli* is the ideal cell host for recombinant expression because of its efficiency and low-cost.

|        File         | Description |
|       :---:         | ----------- |
| FALVAC-1  | ** FALVAC-1 is a vaccine against Plasmodium Falciparum. This is a good benchmark for a synthetic gene as it is a real candidate vaccine that is recombinantly produced ([Reference Paper](https://doi.org/10.1016%2Fj.pep.2003.11.006)).|
| PEA   | ** *Pseudomonas aeruginosa* exotoxin A (PEA) is an important pathogenic factor [1](https://doi.org/10.1073/pnas.85.9.2939). It retains high immunogenicity even after detoxification, enabling its use as vaccine adjuvants and vaccine carriers. This is a good benchmark as it can be produced for a low-cost at a large-scale in *E. coli*. Further, the ([Reference Paper](https://www.sciencedirect.com/science/article/pii/S1046592810000501?via%3Dihub#bib2)) finds that codon optimization enhances expression of PEA in *E. coli* -- thus, if the tool presented in this research can achieve similar/better results than the paper and/or other approaches, it will be considered an improvement. |
| msox   | ** Monomeric sarcosine oxidase (Msox) is an important diagnostic enzyme and can catalyze the oxidation of N-methyl with the FAD cofactor. ([Study](https://www.sciencedirect.com/science/article/pii/S0168165615301905?via%3Dihub)) found that the expression of recombinant SOX "has been extensively studied." Thus, the sox gene serves as a good benchmark because the extent to which codon optimization affects it has been studied. |
| soxB   | ** Sarcosine oxidase subunit B (soxB) is an important diagnostic enzyme and can catalyze the oxidation of N-methyl with the FAD cofactor. ([Study](https://www.sciencedirect.com/science/article/pii/S0168165615301905?via%3Dihub)) found that the expression of recombinant SOX "has been extensively studied." Thus, the sox gene serves as a good benchmark because the extent to which codon optimization affects it has been studied. |
| PRL-3/PTP4A3  | ** Protein tyrosine phosphatase 4A3 (PTP4A3/PRL-3) is a well-documented protein that can be produced recombinantly in *E. coli.* This protein has been identified as a potential target to treat some cancers. |
| hPDF  | ** Human peptide deformylase (hPDF) is a target for cancer therapeutics. However, its expression is not very efficient in *E. coli.* This serves as a good target for benchmarks as past [studies](https://pubmed.ncbi.nlm.nih.gov/19825416/) have noted the valuable potential of codon optimization of this gene. |
| PF3D7  | ** Circumsporozoite protein is used as a surface antigen on the sporozite of a malaria parasite. It has been studied in papers that attempt to quantify gene expression. |
| PA  | ** Polymerase acidic protein (PA) plays a role in viral RNA transcription and replication. It is from the Influenza A virus. It has been studied in ([codon optimization papers](https://www.nature.com/articles/s41598-020-74091-z)).  |
| PIM1  | ** PIM-1 kinase can be found in humans for signal transduction and the development of lymphoid malignancies. Studies have measured the expression of this gene in *E. coli* before, enabling suitable comparisons to be made [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2970903/) [2](https://pubmed.ncbi.nlm.nih.gov/16508102/). |
| MAPKAPK5  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| LCK  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| FGFR4  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| FLT1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| AKT1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| RPS6KB1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| BRAF1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| PAK1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| PLK1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MAPK1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| CSNK1A1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| NPR1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| GSK3B  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| MmpL3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| EMG1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| UBTF  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| PDCD11  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| NOC2L  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| CEBPZ  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| CDK1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| BIRC5  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| SMARCD1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| KIF11  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| NGFR  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| TAS2R10  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| OPRM1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| LEMD3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| CLN3  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| LAMP1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| TAP1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |
| CAV1  | ** Mycobacterial membrane protein Large 3 (MmpL3) is a  |

## Scripts
The following is a description of the purpose for each script in the repository.

`reformat_seqs.py`
> Iterate through each file in a directory and reformat the sequence uniformly.

`naively_optimize.py`
> Naively optimizes a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory.

`benchmark_seqs.ipynb`
> An interactive notebook that helps benchmark a directory containing FASTA sequences across the following metrics:
- Codon Adaptation Index (CAI)
- GC Content
- CFD (known un-optimized gene that reduces efficiency)
- Negative CIS elements
- Negative repeat elements