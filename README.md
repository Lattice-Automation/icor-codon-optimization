<p align="center">
  <img src="/assets/icor-flat-small.png">
</p>

<p align="center">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5529209.svg)](https://doi.org/10.5281/zenodo.5529209)
[![LICENSE](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/Lattice-Automation/icor-codon-optimization/blob/master/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/Lattice-Automation/icor-codon-optimization)

</p>

<h3 align="center"> ICOR: Improving Codon Optimization with Recurrent neural networks <h3>

---
- [About](#About)
- [Quickstart](#Quickstart)
- [Assets](#Assets)
- [Benchmark Results](#Benchmark-Results)
- [Benchmark Sequences](#Benchmark-Sequences)
- [ICOR Tool](#Tool)
  - [Models](#Models)
  - [Optimizers](#Optimizers)
  - [Scripts](#Scripts)
  - [Resources](#Resources)
- [Dependencies](#Dependencies)
- [Training](#Training)
---

### About
In protein sequences—as there are 61 sense codons but only 20 standard amino acids—most amino acids are encoded by more than one codon. Although such synonymous codons do not alter the encoded amino acid sequence, their selection can dramatically affect the production of the resulting protein. Codon optimization of synthetic DNA sequences for maximum expression is an important segment of heterologous expression. However, existing solutions are primarily based on choosing high-frequency codons only, neglecting the important effects of rare codons. In this paper, we propose a novel recurrent-neural-network (RNN) based codon optimization tool, ICOR, that aims to learn codon usage bias on a genomic dataset of Escherichia coli. We compile a dataset of over 42,000 non-redundant, robust genes that are used for deep learning. The model uses a bidirectional long short-term memory-based architecture, allowing for the sequential information of genes to be learnt. Our tool can predict synonymous codons for synthetic genes towards optimal expression in E. coli. We demonstrate that sequential context achieved via RNN may yield codon selection that is more similar to the host genome, therefore improving protein expression more than frequency-based approaches. On a benchmark set of over 40 select DNA sequences, ICOR tool improved the codon adaptation index by 41.69% compared to the original sequence. Our resulting algorithm is provided as an open-source software package along with the benchmark set of sequences.

### Quickstart
Quickstart to run install prerequisites and then launch the ICOR optimization script. Please use the following commands to use our software.

```bash
# Install package
git clone https://github.com/Lattice-Automation/icor-codon-optimization.git

# Install prereqs
pip3 install -r requirements.txt

# Run ICOR optimizer
python3 ./tool/optimizers/icor_optimizer.py
```
Note: ICOR does not support Apple M1/M2 Silicon. This is because Microsoft's onnxruntime does not have support for M1 yet. ICOR should work fine for other MacOSX versions, and of course, Linux & Windows.

Please ensure that scripts are ran from the base folder `icor-codon-optimization` to take advantage of the relative pathing that we have implemented throughout the codebase.

The ICOR software package was developed on Windows with Python v3.9.4 and Linux Ubuntu with Python v3.8.13.
Python v3 is the basic requirement, however, newer versions of Python are suggested for compatibility.

### Assets
Assets including images and branding for the ICOR tool, hosted on the [biotools by Lattice Automation](https://tools.latticeautomation.com/) website.

### Benchmark Results
`benchmark_results` is a folder that contains the following summaries, each in the CSV format:
- `ERC_benchmarks` which consists of the benchmark results for the extended random choice (ERC) optimized sequences.
- `icor_benchmarks` which consists of the benchmark results for the ICOR optimized sequences.
- `BFC_benchmarks` which consists of the benchmark results for the background frequency choice (BFC) optimized sequences.
- `original_benchmarks` which consists of the benchmark results for the original, unoptimized sequences.
- `HFC_benchmarks` which consists of the benchmark results for the highest frequency choice (HFC) optimized sequences.
- `URC_benchmarks` which consists of the benchmark results for the uniform random choice (URC) optimized sequences.
- `genscript_benchmarks` which consists of the benchmark results for the [Genscript Gensmart™](https://www.genscript.com/gensmart-free-gene-codon-optimization.html) optimized sequences.

### Benchmark Sequences
`benchmark_sequences` is a folder that contains sequences for benchmarking purposes, each in the FASTA format:
- `aa` which consists of the 40 original amino acid sequences.
- `all_original` which consists of all 40 amino acid and DNA sequences compiled into one file.
- `ERC` which consists of 40 DNA sequences optimized by the ERC optimizer.
- `dna` which consists of the 40 original DNA sequences.
- `icor` which consists of 40 DNA sequences optimized by the ICOR optimizer.
- `BFC` which consists of 40 DNA sequences optimized by the BFC optimizer.
- `URC` which consists of 40 DNA sequences optimized by the URC optimizer.
- `HFC` which consists of 40 DNA sequences optimzied by the HFC optimizer.
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
`ERC_optimizer.py`
> Extended random choice optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It generates 100,000 sequences using the URC algorithm and chooses the one with the highest CAI.

`icor_optimizer.py`
> ICOR optimizer outputs a text file given a sequence input of amino acids or DNA. It is an interactive Python command-line script. It runs an inference through the ICOR model.
  
`BFC_optimizer.py`
> Background frequency choice optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It selects codons to match the natural frequency that occurs within E. coli. This is what many tools in the industry use as well. This tool/script is built upon the `ecoli_codon_frequencies.csv` file in the summaries directory.
  
`URC_optimizer.py`
> Uniform random choice optimizer creates a directory containing amino acid sequences in the FASTA format and saves these "optimized" / "generated" DNA sequences in a directory. It randomly selects a codon given an amino acid, making it a very naive approach.

`HFC_optimizer.py`
> Chooses the highest frequency codon in E. coli only, having just one codon for every amino acid. This approach should have a CAI of 1.0. Creates a directory containing amino acid sequences in the FASTA format and saves the sequences into this directory. 

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

### Training
Please refer to the training directory for scripts/code regarding the training of the ICOR neural network.
The full dataset as described in our paper is also given in this directory.

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
