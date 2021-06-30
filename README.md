# ICOR Codon Optimization Tool
---
- [Benchmarking](#Benchmarking)
- [Scripts](#Scripts)
- [Summaries](#Summaries)
- [Resources](#Resources)
---
## Benchmarking

### Benchmark Summaries
Benchmark summaries & overviews can be found in the `summaries` folder.

### Benchmark Sequences
`benchmark_sequences` is a folder that contains sequences for benchmarking purposes:
- `aa` which consists of 40 amino acid sequences in the FASTA format.
- `dna` which consists of 40 DNA sequences in the FASTA format.
- `super_naive`  consists of 40 DNA sequences in the FASTA format optimized by the super_naive script.
- `naive` which consists of 40 DNA sequences in the FASTA format optimized by the naive script.

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

`super_naive_optimizer.py`
> Super naive optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It randomly selects a codon given an amino acid, making it a very naive approach.

`naive_optimizer.py`
> Naive optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It selects codons to match the natural frequency that occurs within E. coli. This is what many tools in the industry use as well. This tool/script is built upon the `ecoli_codon_frequencies.csv` file in the summaries directory.

`run_benchmark.ipynb`
> An interactive notebook that helps benchmark a directory containing FASTA sequences across the following metrics:
- Codon Adaptation Index (CAI)
- GC Content
- CFD (known un-optimized gene that reduces efficiency)
- Negative CIS elements
- Negative repeat elements

## Summaries
The following is a description of the purpose for each summary in the summaries folder.

`benchmark_genes.csv`
> Description of each benchmark gene used. Also is above, in README file.

`codon_map.xlsx`
> Contains the codon map used for the AA2Codons dictionary.

`super_naive_benchmarks.csv`
> Contains the benchmarks for super_naively-created sequences.

`naive_benchmarks.csv`
> Contains the benchmarks for naively-created sequences.

`original_benchmarks.csv`
> Contains the benchmarks for the original sequences.

`ICOR_benchmarks.csv`
> Contains the benchmarks for the ICOR-optimized sequences.

`benchmarks_overview.xlsx`
> Contains an overview of the benchmarks, comparing each of the "tools" for each of the benchmarks. This is the sheet to look at if you would like to be able to see the metrics differences between the tools.

`ecoli_codon_frequencies.xlsx`
> Contains the codon frequencies found in E. coli for each amino acid. The naive tool was built upon these frequencies.

## Resources
- Python 3.9.4
- biopython
- numpy
- [AA -> Codons dict](https://www.mathworks.com/help/bioinfo/ref/aa2nt.html)
- selenium
    - Chrome (chromedriver does not seem to work for chromium, needs to use an actual chrome installation)
- web_driver