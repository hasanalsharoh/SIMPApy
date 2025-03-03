import pandas as pd
import numpy as np
import sys
import os
from SIMPApy.preprocess import process_multiomics_data

def test_process_multiomics_data():
    """
    Test the process_multiomics_data function with mock data
    """
    print("Creating mock data...")
    
    # Create sample gene list
    genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
    
    # Create sample samples list
    tm_samples = ['tm_sample1', 'tm_sample2']
    tw_samples = ['tw_sample1', 'tw_sample2']
    samples = tm_samples + tw_samples
    
    # Create mock sigimgsea dataframe
    sigimgsea_data = []
    terms = ['Pathway1', 'Pathway2']
    
    for sample in samples:
        for term in terms:
            # Create different matched and tag genes for each sample/term
            if 'tm' in sample:
                matched_genes = 'GENE1;GENE3;GENE5'
                rna_lead = 'GENE1'
                cnv_lead = 'GENE3'
                dna_lead = 'GENE5'
            else:
                matched_genes = 'GENE2;GENE4'
                rna_lead = 'GENE2'
                cnv_lead = 'GENE4'
                dna_lead = ''
            
            sigimgsea_data.append({
                'sample_name': sample,
                'Term': term,
                'mpes': 1.5,
                'fdr_bh': 0.01,
                'rna_lead_genes': rna_lead,
                'cnv_lead_genes': cnv_lead,
                'dna_lead_genes': dna_lead,
                'rna_nes': 2.5,
                'cnv_nes': 1.8,
                'dna_nes': -2.1,
                'matched_genes': matched_genes
            })
    
    sigimgsea = pd.DataFrame(sigimgsea_data)
    
    # Create mock expression data (tp)
    tp = pd.DataFrame(np.random.uniform(0, 1500, size=(len(genes), len(samples))),
                     index=genes, columns=samples)
    
    # Create mock copy number data (cn)
    cn = pd.DataFrame(np.random.uniform(-2, 2, size=(len(genes), len(samples))),
                     index=genes, columns=samples)
    
    # Create mock methylation data (meth)
    meth = pd.DataFrame(np.random.uniform(0, 1, size=(len(genes), len(samples))),
                     index=genes, columns=samples)
    
    # Create mock population info
    pop_info = pd.DataFrame({
        'key': samples,
        'cancer_type': ['BRCA', 'LUAD', 'BRCA', 'LUAD'],
        'ajcc_pathologic_stage': ['Stage I', 'Stage II', 'Stage III', 'Stage IV']
    })
    
    print("Processing mock data...")
    
    # Run the function
    tmas, twas = process_multiomics_data(sigimgsea, tp, cn, meth, pop_info)
    
    # Basic validation of output
    print("\nValidation Results:")
    print(f"TMAS DataFrame shape: {tmas.shape}")
    print(f"TWAS DataFrame shape: {twas.shape}")
    
    # Check for required columns in output
    required_columns = [
        'sample_name', 'Term', 'gene_name', 'fdr', 'mpes',
        'rna_value', 'cnv_value', 'dna_value',
        'rna_nes', 'cnv_nes', 'dna_nes',
        'rna_tag', 'cnv_tag', 'dna_tag',
        'cancer_type', 'ajcc_pathologic_stage'
    ]
    
    missing_columns = [col for col in required_columns if col not in tmas.columns]
    if missing_columns:
        print(f"ERROR: Missing columns in output: {missing_columns}")
    else:
        print("All required columns present in output")
    
    # Check for tag gene identification
    tm_tag_gene_counts = tmas[['rna_tag', 'cnv_tag', 'dna_tag']].sum()
    tw_tag_gene_counts = twas[['rna_tag', 'cnv_tag', 'dna_tag']].sum()
    
    print("\nTag gene counts in TMAS:")
    print(tm_tag_gene_counts)
    print("\nTag gene counts in TWAS:")
    print(tw_tag_gene_counts)
    
    # Check that the TPM values are clipped at 1000
    tpm_max = tmas['rna_value'].max()
    print(f"\nMaximum TPM value after clipping: {tpm_max}")
    if tpm_max > 1000:
        print("ERROR: TPM values were not clipped properly")
    else:
        print("TPM values clipped correctly")
    
    # Check sample separation
    tm_samples_in_tmas = tmas['sample_name'].unique()
    tw_samples_in_twas = twas['sample_name'].unique()
    
    print(f"\nTMAS samples: {tm_samples_in_tmas}")
    print(f"TWAS samples: {tw_samples_in_twas}")
    
    # Return results for additional inspection if needed
    return tmas, twas

if __name__ == "__main__":
    print("Testing SIMPApy preprocess.process_multiomics_data function")
    tmas, twas = test_process_multiomics_data()
    print("\nTest completed. Inspect the returned DataFrames for further validation.")