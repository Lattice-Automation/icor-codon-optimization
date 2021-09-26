<p align="center">
  <img src="/assets/icor-flat-small.png">
</p>

<p align="center">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5173008.svg)](https://doi.org/10.5281/zenodo.5173008)
[![LICENSE](https://img.shields.io/badge/license-none%20(yet)-brightgreen)](https://github.com/Lattice-Automation/icor-codon-optimization/blob/master/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/Lattice-Automation/icor-codon-optimization)

</p>

<h3 align="center"> ICOR: Improving Codon Optimization with Recurrent neural networks <h4>

---
- [About](#About)
- [Assets](#Assets)
- [Benchmark Results](#Benchmark-Results)
- [Benchmark Sequences](#Benchmark-Sequences)
- [ICOR Tool](#Tool)
  - [Models](#Models)
  - [Optimizers](#Optimizers)
  - [Scripts](#Scripts)
  - [Resources](#Resources)
- [Dependencies](#Dependencies)
---

### Quickstart
I really like having a quickstart section that gives me a single command to install prereqs, a single command to run all tests (if any), and a single command to run the application. Something like:

```bash
# Install prereqs
pip install -r requriements.txt # or an install_prereqs.sh script if you have more diverse dependencies

# run tests (if you decided to add tests in the future)
pytest

# run models
python ./tool/optimizers/brute_force_optimizer.py
```

### About
In protein sequences—as there are 61 sense codons but only 20 standard amino acids—most amino acids are encoded by more than one codon. Although such synonymous codons do not alter the encoded amino acid sequence, their selection can dramatically affect the production of the resulting protein. Codon optimization of synthetic DNA sequences for maximum expression is an important segment of heterologous expression. However, existing solutions are primarily based on choosing high-frequency codons only, neglecting the important effects of rare codons. In this paper, we propose a novel recurrent-neural-network (RNN) based codon optimization tool, ICOR, that aims to learn codon usage bias on a genomic dataset of Escherichia coli. We compile a dataset of over 42,000 non-redundant, robust genes that are used for deep learning. The model uses a bidirectional long short-term memory-based architecture, allowing for the sequential information of genes to be learnt. Our tool can predict synonymous codons for synthetic genes towards optimal expression in E. coli. We demonstrate that sequential context achieved via RNN may yield codon selection that is more similar to the host genome, therefore improving protein expression more than frequency-based approaches. On a benchmark set of over 40 select DNA sequences, ICOR tool improved the codon adaptation index by 41.69% compared to the original sequence. Our resulting algorithm is provided as an open-source software package along with the benchmark set of sequences.

### Assets
Assets including images and branding for the ICOR tool, hosted on the [biotools by Lattice Automation](https://tools.latticeautomation.com/) website.

### Benchmark Results
`benchmark_results` is a folder that contains the following summaries, each in the CSV format:
- `brute_benchmarks` which consists of the benchmark results for the brute force optimized sequences.
- `icor_benchmarks` which consists of the benchmark results for the ICOR optimized sequences.
- `naive_benchmarks` which consists of the benchmark results for the naively optimized sequences.
- `original_benchmarks` which consists of the benchmark results for the original, unoptimized sequences.
- `super_naive_benchmarks` which consists of the benchmark results for the super naively optimized sequences.
- `genscript_benchmarks` which consists of the benchmark results for the [Genscript Gensmart™](https://www.genscript.com/gensmart-free-gene-codon-optimization.html) optimized sequences.

### Benchmark Sequences
`benchmark_sequences` is a folder that contains sequences for benchmarking purposes, each in the FASTA format:
- `aa` which consists of the 40 original amino acid sequences.
- `all_original` which consists of all 40 amino acid and DNA sequences compiled into one file.
- `brute` which consists of 40 DNA sequences optimized by the brute force optimizer.
- `dna` which consists of the 40 original DNA sequences.
- `icor` which consists of 40 DNA sequences optimized by the ICOR optimizer.
- `naive` which consists of 40 DNA sequences optimized by the naive optimizer.
- `super_naive`  consists of 40 DNA sequences optimized by the super naive optimizer.
- `genscript` consists of 40 DNA sequences optimized by the [Genscript Gensmart™](https://www.genscript.com/gensmart-free-gene-codon-optimization.html) tool.

### Tool
The ICOR tool has been divided into four directories: models, optimizers, resources, and scripts. In the `/tool/optimizers` directory sits the `icor_optimizer.py` file: an interactive script to optimize a sequence utilizing the trained ICOR model.

Supporting files were used to train, evaluate, and test the ICOR model. Descriptions for these can be found below:

#### Models
The models directory contains the trained ICOR model in the [ONNX](https://onnx.ai) (open-neural-network-exchange) format. Below is a preview of the model architecture:

<img src="/assets/icor-small-visualization.png">
The ICOR model was trained in the MATLAB environment. For more details on model architecture, please review our manuscript file in the base of the repository. Upon submission, this will be changed to a DOI/biorxiv link.

`benchmark_genes.pdf`
> A document that contains all of the benchmarking genes and descriptions of them.

#### Optimizers
`brute_force_optimizer.py`
> Brute force optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It generates 10,000 sequences and chooses the one with the highest CAI.

`icor_optimizer.py`
> ICOR optimizer outputs a text file given a sequence input of amino acids or DNA. It is an interactive Python command-line script. It runs an inference through the ICOR model.
  
`naive_optimizer.py`
> Naive optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It selects codons to match the natural frequency that occurs within E. coli. This is what many tools in the industry use as well. This tool/script is built upon the `ecoli_codon_frequencies.csv` file in the summaries directory.
  
`super_naive_optimizer.py`
> Super naive optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It randomly selects a codon given an amino acid, making it a very naive approach.

#### Scripts
The following is a description of the purpose for each script in the repository.

`convert_to_cds.py`  
> Takes an input of DNA sequences and fetches their CDS only from the NCBI nuccore database. Rewrites files with CDS.

`csv_to_seqs.py`
> Takes an input of a CSV from the GenScript Gensmart tool and writes them into files containing the sequences in the FASTA format.
  
`reformat_seqs.py`
> Iterate through each file in a directory and reformat the sequence uniformly.

`run_benchmark.ipynb`
> An interactive notebook that helps benchmark a directory containing FASTA sequences across the following metrics:
- Codon Adaptation Index (CAI)
- GC Content
- CFD (known un-optimized gene that reduces efficiency)
- Negative CIS elements
- Negative repeat elements
  
#### Resources
The following is a description of the purpose for each resource in the resources folder.

`Benchmarking Results & Comparison - ICOR Codon Optimization.pdf`
> Contains an overview of the benchmarks, comparing each of the "tools" for each of the benchmarks. This is the sheet to look at if you would like to be able to see the metrics differences between the tools.
  
`benchmark_genes.pdf`
> A table for all of the benchmark genes used for validation.
  
`codon_map.xlsx`
> Contains the codon map used for the AA2Codons dictionary.
  
`ecoli_codon_frequencies.csv`
> Contains the codon frequency weights for each codon/amino acid used in the E. coli genomes. The naive tool was built upon these frequencies.

### Dependencies
- Python 3.9.4
  - biopython
  - numpy
  - web_driver
  - onnxruntime
  - re
  - selenium
    - Chrome (chromedriver does not seem to work for chromium, needs to use an actual chrome installation)
- [AA -> Codons dict](https://www.mathworks.com/help/bioinfo/ref/aa2nt.html)
