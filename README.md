# SIMPApy

Normalized Single Sample Integrated Multi-Omics Pathway Analysis for Python

## Description

`SIMPApy` is a Python package for performing Gene Set Enrichment Analysis (GSEA) on multi-omics data and integrating the results. It supports RNA sequencing, DNA methylation, and copy number variation data types.

## Installation

```bash
pip install SIMPApy
```

## Features

- Run GSEA on single-omics data for single samples normalized to a control population (SOPA).
- Calculate normalized single sample gene rankings for different omics data types (RNA-seq, DNA methylation, CNV).
- Integrate GSEA results from multiple omics platforms in control-normalized single samples (SIMPA).
- Calculate Multiomics Pathway Enrichment Score (MPES) for analysis of differentially activated pathways.

## Usage

### Calculate Rankings for Different Omics Types

```python
import SIMPApy as sp
import pandas as pd

# Load example data found in data directory
rna_data = pd.read_csv("rna.csv", index_col=0)
cnv_data = pd.read_csv("cn.csv", index_col=0)
dna_data = pd.read_csv("dna.csv", index_col=0) # M-values
meth_data = pd.read_csv("meth.csv", index_col=0) # beta values

hallmark = "path/to/h.all.v2023.1.Hs.symbols.gmt" # hallmarks gene set

# Calculate rankings for RNA-seq data
rna_rankings = sp.calculate_ranking(
    df=rna_data,
    omic="RNA",
    alpha=0.05
)

# Calculate rankings for CNV data
cnv_data = pd.read_csv("cnv_data.csv", index_col=0)
cnv_rankings = sp.calculate_ranking(
    df=cnv_data,
    omic="CNV"
)
```

### Integrate Multi-Omics GSEA Results

```python
import SIMPApy as sp

# Run SIMPA for a single sample
integrated_results = sp.simpa(
    sample_id="sample1",
    rna_dir="path/to/rna/gsea/results",
    cnv_dir="path/to/cnv/gsea/results",
    dna_dir="path/to/dna/gsea/results",
    output_dir="path/to/output"
)

# Process multiple samples
sample_ids = ["sample1", "sample2", "sample3"]
sp.run_simpa_batch(
    sample_ids=sample_ids,
    rna_dir="path/to/rna/gsea/results",
    cnv_dir="path/to/cnv/gsea/results",
    dna_dir="path/to/dna/gsea/results",
    output_dir="path/to/output"
)
```

## Requirements

- Python ≥ 3.8
- gseapy==1.1.3
- numpy==1.23.5
- scipy==1.14.0
- pandas==2.2.2
- pydeseq2==0.4.10
- matplotlib==3.9.2
- plotly==5.24.1
- scikit-learn==1.5.1
- seaborn==0.13.2
- statsmodels==0.14.1
- ipywidgets==8.1.5
- pillow==10.4.0


## License

Apache-2.0