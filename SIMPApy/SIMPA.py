"""
Integration module for normalized single sample integrated multiomics pathway analysis.

This module provides functions to integrate SOPA results from multiple
omics data types (RNA-seq, CNV, and DNA methylation) to identify
consistently enriched pathways.
"""

import os
import pandas as pd
import numpy as np
import glob
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from typing import Dict, List, Union, Optional, Tuple, Any


def calculate_wcos_mpes(row: pd.Series) -> pd.Series:
    """
    Calculates Weighted Combined Omics Score (WCOS) and Multiomics Pathway Enrichment Score (MPES) 
    for a single pathway across multiple omics platforms.
    
    Args:
        row: A pandas Series containing GSEA results for a pathway across omics types.
        
    Returns:
        A pandas Series with WCOS values for each omic type and the MPES score.
    """
    wcos_values = []
    for omic in ['rna', 'cnv', 'dna']:
        fdr = row[f'{omic}_fdr']
        nes = row[f'{omic}_nes']

        # Handle cases where pval (and thus FDR) might be exactly 0 or 1
        fdr = max(1e-16, min(fdr, 1 - 1e-16))

        # Correctly count leading-edge genes from the lead_genes string
        leading_edge_genes = len(row[f'{omic}_lead_genes'].split(';')) if isinstance(row[f'{omic}_lead_genes'], str) else 0
        
        matched_genes = len(row['matched_genes'].split(';')) if isinstance(row['matched_genes'], str) else 0

        l = leading_edge_genes / matched_genes if matched_genes > 0 else 0
        wcos = (1 - fdr) * nes * np.log(1+l)
        wcos_values.append(wcos)

    # Normalize WCOS values
    wcos_mean = np.mean(wcos_values)
    wcos_std = np.std(wcos_values)

    # Calculate MPES
    mpes = ((np.sum(wcos_values)) - wcos_mean) / np.sqrt((wcos_std**2)/3)

    return pd.Series({
        'rna_wcos': wcos_values[0],
        'cnv_wcos': wcos_values[1],
        'dna_wcos': wcos_values[2],
        'mpes': mpes
    })


def sort_gene_list(gene_str: str) -> str:
    """
    Sort a semicolon-separated string of genes alphabetically.
    
    Args:
        gene_str: String of semicolon-separated gene names.
        
    Returns:
        Sorted string of semicolon-separated gene names.
    """
    if pd.isna(gene_str) or not isinstance(gene_str, str):
        return gene_str
    genes = gene_str.split(';')
    return ';'.join(sorted(genes))


def simpa(
    sample_id: str, 
    rna_dir: str, 
    cnv_dir: str, 
    dna_dir: str,
    output_dir: Optional[str] = None
) -> Optional[pd.DataFrame]:
    """
    Integrates GSEA results from RNA, CNV, and DNA methylation for a single sample.
    
    Args:
        sample_id: Sample identifier used in filenames.
        rna_dir: Directory containing RNA-seq GSEA results.
        cnv_dir: Directory containing CNV GSEA results.
        dna_dir: Directory containing DNA methylation GSEA results.
        output_dir: Optional directory to save integrated results. If None, results are not saved.
        
    Returns:
        A pandas DataFrame with integrated GSEA results or None if an error occurs.
    """
    try:
        # Read files and ensure Term is not the index
        rna = pd.read_csv(os.path.join(rna_dir, f'{sample_id}_gsea_results.csv'))
        cnv = pd.read_csv(os.path.join(cnv_dir, f'{sample_id}_gsea_results.csv'))
        dna = pd.read_csv(os.path.join(dna_dir, f'{sample_id}_gsea_results.csv'))
    except FileNotFoundError:
        print(f"GSEA results file not found for sample: {sample_id}")
        return None
    except pd.errors.EmptyDataError:
        print(f"GSEA results file is empty for sample: {sample_id}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading GSEA results for sample {sample_id}: {e}")
        return None

    # Sort matched_genes alphabetically in each dataframe
    for df in [rna, cnv, dna]:
        df['matched_genes'] = df['matched_genes'].apply(sort_gene_list)

    # Rename columns with omic prefixes, excluding 'Term' and 'matched_genes'
    for df, prefix in [(rna, 'rna'), (cnv, 'cnv'), (dna, 'dna')]:
        df.columns = [f'{prefix}_{col}' if col not in ['Term', 'matched_genes'] else col 
                     for col in df.columns]
        # Fix the (leading edge genes) column name
        if f'{prefix}_lead_genes' in df.columns:
            df.rename(columns={f'{prefix}_lead_genes': f'{prefix}_lead_genes'}, inplace=True)

    # Replace 0 and 1 values with 1e-16 and 1-1e-16 in pval columns
    for df, prefix in [(rna, 'rna'), (cnv, 'cnv'), (dna, 'dna')]:
        df[f'{prefix}_pval'] = df[f'{prefix}_pval'].replace({0.000000: 1e-16, 1.000000: 1 - 1e-16})

    # Merge dataframes on Term
    combined = rna.merge(cnv, on=['Term', 'matched_genes'], how='inner')
    combined = combined.merge(dna, on=['Term', 'matched_genes'], how='inner')

    # Calculate z-scores
    for prefix in ['rna', 'cnv', 'dna']:
        combined[f'{prefix}_z'] = norm.ppf(1 - combined[f'{prefix}_pval'])
        # Multiply z-scores by the sign of NES to get the correct direction
        combined[f'{prefix}_z'] = combined[f'{prefix}_z'] * np.sign(combined[f'{prefix}_nes'])

    # Calculate combined z-score
    z_scores = combined[['rna_z', 'cnv_z', 'dna_z']]
    combined['combined_z'] = z_scores.sum(axis=1) / np.sqrt(3)

    # Calculate combined p-value
    combined['combined_pval'] = norm.sf(abs(combined['combined_z'])) * 2  # Two-tailed test

    # Multiple testing correction with proper error handling
    try:
        # Remove any infinite or NA values
        valid_pvals = combined['combined_pval'].replace([np.inf, -np.inf], np.nan).dropna()
        
        if len(valid_pvals) > 0:
            # Perform correction only on valid p-values
            reject, pvals_corrected, _, _ = multipletests(valid_pvals, method='fdr_bh')
            
            # Initialize FDR column with NaN
            combined['fdr_bh'] = np.nan
            
            # Update FDR values only for rows with valid p-values
            combined.loc[valid_pvals.index, 'fdr_bh'] = pvals_corrected
        else:
            # If no valid p-values, set all FDR to NaN
            combined['fdr_bh'] = np.nan
            print(f"Warning: No valid p-values for multiple testing correction in sample {sample_id}")
    except Exception as e:
        print(f"Warning: Error in multiple testing correction for sample {sample_id}: {e}")
        combined['fdr_bh'] = np.nan

    # Sort by FDR (putting NaN values at the end)
    combined = combined.sort_values(by='fdr_bh', na_position='last')

    # Apply WCOS and MPES calculations
    wcos_mpes_results = combined.apply(calculate_wcos_mpes, axis=1)
    combined = pd.concat([combined, wcos_mpes_results], axis=1)

    # Keep relevant columns
    try:
        combined = combined[['Term', 'combined_pval', 'combined_z', 'fdr_bh', 'matched_genes',
                          'rna_lead_genes', 'cnv_lead_genes', 'dna_lead_genes',
                          'rna_nes', 'cnv_nes', 'dna_nes', 
                          'rna_wcos', 'cnv_wcos', 'dna_wcos', 'mpes']]
    except KeyError as e:
        # Handle case where columns might have different names
        print(f"Warning: Some expected columns are missing: {e}")
        # Keep all columns if we can't find the expected ones
        pass

    # Save results if output_dir is provided
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f'{sample_id}_integrated_gsea_results.csv')
        combined.to_csv(output_file)
        print(f"Saved integrated results for {sample_id} to {output_file}")

    return combined


def run_simpa_batch(
    sample_ids: List[str],
    rna_dir: str,
    cnv_dir: str,
    dna_dir: str,
    output_dir: str
) -> None:
    """
    Run SIMPA integration for multiple samples.
    
    Args:
        sample_ids: List of sample identifiers.
        rna_dir: Directory containing RNA-seq GSEA results.
        cnv_dir: Directory containing CNV GSEA results.
        dna_dir: Directory containing DNA methylation GSEA results.
        output_dir: Directory to save integrated results.
        
    Returns:
        None. Results are saved to files in the output directory.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for sample_id in sample_ids:
        combined_df = simpa(sample_id, rna_dir, cnv_dir, dna_dir)
        
        if combined_df is not None:
            output_file = os.path.join(output_dir, f'{sample_id}_integrated_gsea_results.csv')
            combined_df.to_csv(output_file)

        # Clear memory
        del combined_df
        if 'rna' in locals():
            del rna
        if 'cnv' in locals():
            del cnv
        if 'dna' in locals():
            del dna

    print('Integration done! Results saved in:', output_dir)

def load_simpa(directory):
    """
    Loads and processes SIMPA results from a directory of CSV files.

    Args:
        directory: The path to the directory containing the SIMPA results files.

    Returns:
        A pandas DataFrame with columns: sample_name, term, fdr, pval.
    """

    all_results = []

    # Use glob to find all CSV files matching the pattern
    file_pattern = os.path.join(directory, "tm*_integrated_gsea_results.csv")
    file_paths = glob.glob(file_pattern)
    
    # also get tw files
    file_pattern = os.path.join(directory, "tw*_integrated_gsea_results.csv")
    file_paths.extend(glob.glob(file_pattern))

    for file_path in file_paths:
        # Extract sample name from filename
        file_name = os.path.basename(file_path)
        sample_name = file_name.split("_integrated_gsea_results")[0]  # Extract tm(n) or tw(n)

        # Load the CSV file into a DataFrame
        try:
            df = pd.read_csv(file_path)
        except pd.errors.ParserError:
            print(f"Error: Could not parse {file_path} as a CSV file. Skipping.")
            continue

        # Select relevant columns and add sample name
        df = df[['Term','combined_pval', 'combined_z', 'fdr_bh','matched_genes',
                          'rna_lead_genes', 'cnv_lead_genes', 'dna_lead_genes', 'rna_nes', 'cnv_nes', 'dna_nes', 'mpes']]
        df['sample_name'] = sample_name

        # Append to the list of results
        all_results.append(df)

    # Concatenate all results into a single DataFrame
    final_df = pd.concat(all_results, ignore_index=True)

    return final_df