<p align="center">
  <img src="/assets/icor-flat-small.png">
</p>

# ICOR: Improving Codon Optimization with Recurrent neural networks
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5173008.svg)](https://doi.org/10.5281/zenodo.5173008)

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

`benchmark_genes.pdf`
> A document that contains all of the benchmarking genes and descriptions of them.

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

`Benchmarking Results & Comparison - ICOR Codon Optimization.pdf`
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
