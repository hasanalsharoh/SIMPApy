"""
Core functions for the SIMPApy package.

This module contains the main functions for running SOPA on ranked gene data.
"""

import gseapy as gp
import pandas as pd
import numpy as np
import os
import glob
from typing import Dict, List, Union, Optional, Tuple
import logging # Import the logging library

def _sopa(
    ranking: pd.Series,
    gene_set: Union[Dict, str],
    minisz: int = 3,
    seeder: int = 7,
    threads: int = 8,
    permutation_num: int = 1000,
    **kwargs
) -> pd.DataFrame:
    """
    Run SOPA on a ranked gene list.
    
    Args:
        ranking: A pandas Series with gene names as index and ranking values.
        gene_set: Gene set database in GMT format or a dictionary.
        minisz: Minimum size of gene sets to consider. Default is 3.
        seeder: Random seed for reproducibility. Default is 7.
        threads: Number of threads to use. Default is 8.
        permutation_num: Number of permutations for calculating FDR. Default is 1000.
        **kwargs: Additional arguments passed to gp.prerank().
        
    Returns:
        A pandas DataFrame with SOPA GSEA results sorted by FDR.
    """
    # Pass all arguments to prerank, including any additional arguments
    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=gene_set,
        min_size=minisz,
        seed=seeder,
        threads=threads,
        **kwargs
    )
    
    out = []
    for term in list(pre_res.results):
        out.append([
            term,
            pre_res.results[term]['fdr'],
            pre_res.results[term]['es'],
            pre_res.results[term]['nes'],
            pre_res.results[term]['pval'],
            pre_res.results[term]['matched_genes'],
            pre_res.results[term]['gene %'],
            pre_res.results[term]['lead_genes'],
            pre_res.results[term]['tag %']
        ])
    
    out_df = pd.DataFrame(
        out,
        columns=['Term', 'fdr', 'es', 'nes', 'pval', 'matched_genes', 'gene %', 'lead_genes', 'tag %']
    ).sort_values('fdr').reset_index(drop=True)
    
    return out_df


def sopa(
    ranks: pd.DataFrame,
    gene_set: Union[Dict, str],
    output_dir: str,
    minisz: int = 3,
    seeder: int = 7,
    **kwargs
) -> None:
    """
    Run SOPA on all samples in the ranking dataframe and save results as CSV files.
    
    Args:
        ranks: DataFrame with genes as index and samples as columns.
        gene_set: Gene set database in GMT format or a dictionary.
        output_dir: Directory where results will be saved.
        minisz: Minimum size of gene sets to consider. Default is 3.
        seeder: Random seed for reproducibility. Default is 7.
        **kwargs: Additional arguments passed to sopa().
        
    Returns:
        None. Results are saved to files in the output directory.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate through each sample in the ranks DataFrame
    for col in ranks.columns:
        # Sort rankings and handle infinities
        ranking = ranks[col].replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
                
        # Run sopa with all the parameters
        gsea_result = _sopa(
            ranking=ranking,
            gene_set=gene_set,
            minisz=minisz,
            seeder=seeder,
            **kwargs
        )

        # Construct the output file path (use os.path.join for cross-platform compatibility)
        output_file = os.path.join(output_dir, f"{col}_gsea_results.csv")

        # Save the GSEA results to a CSV file
        gsea_result.to_csv(output_file, sep=',')

        # Clean up
        del gsea_result, ranking

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def load_sopa(directory: str) -> pd.DataFrame:
    """
    Loads and processes sopa results from a directory of CSV files.

    Args:
        directory: The path to the directory containing the SOPA results files.

    Returns:
        A pandas DataFrame with columns: sample_name, term, fdr, pval.
        Returns an empty DataFrame if no valid files are found or processed.
    """

    all_results = []
    required_columns = ['Term', 'fdr', 'es', 'nes', 'matched_genes', 'gene %', 'tag %', 'lead_genes']

    # Use glob to find all CSV files matching the patterns
    file_pattern_tm = os.path.join(directory, "tm*_gsea_results.csv")
    file_pattern_tw = os.path.join(directory, "tw*_gsea_results.csv")
    file_paths = glob.glob(file_pattern_tm) + glob.glob(file_pattern_tw) # Combine lists directly

    if not file_paths:
        logging.warning(f"No files matching 'tm*_gsea_results.csv' or 'tw*_gsea_results.csv' found in directory: {directory}")
        return pd.DataFrame(columns=['sample_name'] + required_columns) # Return empty DataFrame with expected columns

    for file_path in file_paths:
        # Extract sample name from filename
        file_name = os.path.basename(file_path)
        sample_name = file_name.split("_gsea_results")[0]  # Extract tm(n) or tw(n)

        try:
            # Load the CSV file into a DataFrame
            df = pd.read_csv(file_path)

            # Check for and select required columns
            try:
                df_selected = df[required_columns].copy() # Use .copy() to avoid SettingWithCopyWarning
                df_selected['sample_name'] = sample_name
                all_results.append(df_selected)
            except KeyError as e:
                # Handle missing columns after successful parsing
                missing_cols = list(set(required_columns) - set(df.columns))
                logging.error(f"Error: File {file_path} is missing required columns: {missing_cols}. Skipping.")
                continue # Skip this file

        except pd.errors.ParserError:
            # Handle files that cannot be parsed as CSV
            logging.error(f"Error: Could not parse {file_path} as a CSV file. Skipping.")
            continue
        except FileNotFoundError:
            # Handle case where file disappears between glob and read_csv (unlikely but possible)
            logging.error(f"Error: File {file_path} not found. Skipping.")
            continue
        except Exception as e:
             # Catch any other unexpected errors during file processing
             logging.error(f"An unexpected error occurred processing {file_path}: {e}. Skipping.")
             continue


    # Concatenate all results into a single DataFrame
    if not all_results:
        logging.warning(f"No valid SOPA result files were processed successfully in directory: {directory}")
        return pd.DataFrame(columns=['sample_name'] + required_columns) # Return empty DataFrame if no files were valid

    final_df = pd.concat(all_results, ignore_index=True)

    return final_df
