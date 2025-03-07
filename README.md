# SIMPApy

Normalized Single Sample Integrated Multi-Omics Pathway Analysis for Python

## Description

`SIMPApy` is a Python package for performing Gene Set Enrichment Analysis (GSEA) on multi-omics data and integrating the results. It supports RNA sequencing, DNA methylation, and copy number variation data types.

## Installation


To install SIMPApy, create a new virtual environment making sure to install the requirements found at the SIMPApy Git repository.
Afterwards, use:
```bash
pip install SIMPApy 
```

## Features

- Run GSEA on single-omics data for single samples normalized to a control population (SOPA).
- Calculate normalized single sample gene rankings for different omics data types (RNA-seq, DNA methylation, CNV).
- Integrate GSEA results from multiple omics platforms in control-normalized single samples (SIMPA).
- Calculate Multiomics Pathway Enrichment Score (MPES) for analysis of differentially activated pathways.

## Usage

### Description of SOPA and SIMPA methodology
To use SOPA, we need raw data for a single -omic, and as follows:
- The dataframes should be:
    1. For RNAseq data: TPM, or normalized counts. FPKM, RPKM may work, but the tool was not validated on them.
    2. For CNV data: Copy numbers. Each gene must have its individual copy numbers.
    3. For DNAm data: Beta values. The values must be mapped to gene names rather than CpG sites.
- The dataframes must contain sample divided into 2 groups, cases and controls. For the case group, each sample must start with the name 'tm', while for the control group, each sample must start with 'tw'.
- It is possible to use more than 2 groups. For example:
        
        Assume 3 groups (group1, group2, group3). Analyses could be done between group1 and group2, and then group3 and group2, or group3 and group1, according to preference.

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


## Single sample gene ranking
### RNAseq and DNAm M-values
```python
rnaranks = sp.calculate_ranking(rna,omic='rna', alpha=0.05)
# assuming the following:
# df: pandas DataFrame with gene expression data (genes as rows, samples as columns).
```
After running the above functions, the following code should be used to retrieve the dataset:
```python
# make a dataframe of rnaranks where index is gene names and columns are samples and only include "weighted" column
rnaranks_df = pd.concat({k: v['weighted'] for k, v in rnaranks.items()}, axis=1)
```
### CNV data
Here, make sure to have copy number data. If GISTIC2 data is present (often through log(Copy_Number +1)), then conversion is usually done through 2**(GISTIC2_value -1)
```python
cnvs = sp.calculate_ranking(cn, omic = 'cnv'):
# cn: Pandas DataFrame with gene-level copy numbers. 
# Rows are genes and columns are samples ('tm' for cases, 'tw' for controls).
```
Afterwards, retrieve the dataframe:
```python
# retrieve the final dataframe with weights
cnranks_df = pd.concat({k: v['adjusted_weight'] for k, v in cnranks.items()}, axis=1)
```
## Running SOPA
To run SOPA, we need 2 available files:

- A single sample gene ranking dataset
- a number of gene sets as a .gmt file
```python
# to run SOPA on a single sample
single_sample_results = sp.sopa(ranking, gene_set, minisz, seeder) # assuming is a dataframe of 1 column of ranking values, and rows (indices) as genes
# minisz is the minimal gene set matching size
# seeder is the random seed to maintain result consistency
```
if multiple samples are available in the single sample dataframe, each column must be a sample with ranking values as rows. Gene names must be ENSEMBL IDs. We could run SOPA after defining the function above and running the code:
```python
single_samples_output = sp.sopa_population(ranking_dataframe, gene_set_gmt_file, folder_to_contain_outputs_for_single_sample_enrichment_analysis)
```
## SIMPA
To run SIMPA, we need to have RNAseq, CNV, and DNA methylation SOPA results in 3 different folders.
First, we get file names. This is to be done once, as the output of SOPA is similar across all 3 omics
```python
# get file list
file_list= glob.glob(os.path.join(dir, '*_gsea_results.csv')) # dir is to be replaced with r"X:\sopa_results_folder_Location\".
# sample ids
sample_ids = [os.path.basename(f).split('_')[0] for f in file_list]
```
Then, we define the paths to RNA, CNVs, and DNAm SOPA results:
```python
# This function requires 3 data directories to run, one for each omic type
# rna_dir = 'path/to/rna/sopa/results'
# cnv_dir = 'path/to/cnv/sopa/results'
# dna_dir = 'path/to/dna/sopa/results'

# run SIMPA
simpa_res = sp.run_simpa_batch(rna_dir, cnv_dir, dna_dir, output_directory)

```
### Visualization of SIMPA results
Prior to visualization using 3D box, we need 4 - 5 datasets, and as follows:
- Mandatory:
    - SIMPA results dataframe from running run_simpa_batch
    - RNAseq raw data (preferably TPM values, as algorithm clips values at 1,000 TPM)
    - CNVs raw data (copy numbers)
    - DNAm raw data (beta values)
- Optional:
    - Clinical information dataset

Then, we retrieve 2 datasets, one for cases (TMAs here) and one for controls (TWAs in this example):

```python
tmas, twas = sp.process_multiomics_data(simpa_res, # SIMPA results
                                            rna, # TPM values
                                             cn, # copy numbers
                                             dna, # beta values
                                             pop_info) # for samples clinical information
```
Afterwards, we could use:

```python
sp.create_interactive_plot(twas, # loading dataframe 
                            "TWA") # header title
```

The could would show the figure below, which can be rotated and interacted with

![.](https://github.com/hasanalsharoh/SIMPApy/blob/main/images/3D_interactive_plot.png)

In the figure, each dot is a single gene in a single sample. On hovering each datapoint, the following data is shown:
            
            - Gene: Gene name of hovered datapoint
            
            - Sample: Sample name datapoint belongs to
            
            - Term: Term the datapoint enriches (in general, the Term/pathway would match the dropdown menu's Term)
            
            - Cancer Type: clinical information (if available)
            
            - AJCC Stage: clinical (if available) 
            
            - DNAm: DNA methylation beta value
            
            - TPM: TPM value (clipped to 1,000)
            
            - CNV: Copy number value
            
            - Enrichment: NES or MPES
            
            - Tagged: 1 if leading edge (impactful on enrichment), 0 if not. This is determined by the GSEA algorithm
            
            - Multiomics FDR: BH-corrected p-value for the Term/pathway the datapoint belongs to.

Figure settings allow the following:
    
    - Term (dropdown menu): Filtering based on pathway/term
    
    - Omic Type (dropdown menu): Changes color based on enrichment of the selected omic (options are RNA, CNV, DNA, MPES)
    
    - Cancer Type (dropdown menu): Filtering based on clinical information
    
    - AJCC Stage (dropdown menu): Filtering based on clinical information
    
    - Filter by (buttons): This setting is used prior to search in the Search box. If sample is selected, the search box filters all datapoints to match the sample searched for. If gene selected, the search box will filter the datapoints to include only genes matching the gene name searched for.
    
    - Search (search box): After the selection of Sample/Gene from the *Filter by* option, the search box will retrieve matching variables (not case sensitive). After writing, click *Search*.
    
    - Normalize Untagged (checkbox): This setting allows for normalizing untagged datapoints (Tagged = 0), by controlling their NES values to 0. This option is mainly for single -omics, as on MPES, all -omics must be tagged for a single datapoint for it to not be normalized.
    
    - Reset View (button): Restores the view to unfiltered.
    
    - Export (300DPI): This options exports the current view of the 3D box.
            
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