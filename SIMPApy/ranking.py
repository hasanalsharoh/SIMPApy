"""
Gene ranking functions for different omics data types.

This module contains functions to calculate rankings for RNA-seq, DNA methylation,
and copy number variation data.
"""

import pandas as pd
import numpy as np
from scipy.stats import norm
from typing import Dict, List, Union, Optional, Tuple


def _calculate_msd(df: pd.DataFrame, alpha: float = 0.05) -> pd.Series:
    """
    Calculates the Minimum Significant Difference (MSD) for each gene in the dataframe.

    Args:
        df: pandas DataFrame with gene expression data.
        alpha: Significance level for the Z-score. Default is 0.05.

    Returns:
        pandas Series with MSD for each gene.
    """
    # Separate TWA (control) group columns
    twa_cols = [col for col in df.columns if col.startswith('tw')]
    twa_df = df[twa_cols]

    # Calculate the mean expression for each gene across the TWA group
    gene_means = twa_df.mean(axis=1)

    # Calculate the Sum of Squares Within (SSW) for each gene
    ssw = ((twa_df.subtract(gene_means, axis=0))**2).sum(axis=1)

    # Calculate the standard error (SE)
    n = len(twa_cols)  # Number of samples in the TWA group
    se = np.sqrt(ssw / (n - 1))

    # Calculate the Z-score for the given alpha level
    z_alpha = norm.ppf(1 - alpha/2)

    # Calculate MSD
    msd = z_alpha * se

    return msd


def calculate_ranking(
    df: pd.DataFrame, 
    omic: str = "RNA", 
    alpha: float = 0.05
) -> Dict[str, pd.DataFrame]:
    """
    Calculate rankings for different types of omics data.
    
    Args:
        df: pandas DataFrame with omics data. Rows are genes/features, columns are samples.
        omic: Type of omics data. Must be "RNA", "DNAm", or "CNV". Default is "RNA".
        alpha: Significance level for RNA and DNAm rankings. Default is 0.05.
        
    Returns:
        A dictionary of DataFrames, where each key is a sample name containing a 
        DataFrame with gene rankings.
    """
    if omic.upper() in ["RNA", "DNAM"]:
        # Calculate MSD first for RNA and DNAm
        msd = _calculate_msd(df, alpha)
        
        # Separate TWA (control) group columns
        twa_cols = [col for col in df.columns if col.startswith('tw')]
        twa_df = df[twa_cols]

        # Calculate the mean expression for each gene across the TWA group
        gene_means = twa_df.mean(axis=1)

        # Dictionary to store the ranked DataFrames
        ranked_dfs = {}

        # Iterate over each sample (column) including TWA samples
        for sample in df.columns:
            # Calculate the difference (D_(x,s))
            d_xs = df[sample] - gene_means

            # Adjust sign of MSD based on D_(x,s)
            msd_signed = msd * np.sign(d_xs)

            # Calculate weighted score
            weighted_score = d_xs / msd

            # Create a DataFrame for the current sample
            sample_df = pd.DataFrame({
                'D_xs': d_xs, 
                'MSD': msd_signed, 
                'weighted': weighted_score, 
                'Significant': abs(d_xs) > msd
            })

            # Rank genes based on D_(x,s)
            sample_df['Rank'] = sample_df['D_xs'].rank(ascending=False)

            # Store the DataFrame in the dictionary
            ranked_dfs[sample] = sample_df

            del sample_df

        return ranked_dfs
    
    elif omic.upper() == "CNV":

        control_data = df.filter(regex='^tw')
        N = len(control_data.columns)
        epsilon = 0.01 # Small constant to prevent division by zero

        # Pre-compute all necessary stats for the control group
        control_counts_df = control_data.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
        mu_controls = control_data.mean(axis=1)
        sigma_controls = control_data.std(axis=1)

        ranked_dfs = {}

        # 2. Loop through all samples
        for sample_name in df.columns:
            sample_series = df[sample_name]
            scores = []

            # Loop through each gene in the current sample
            for gene, cn_value in sample_series.items():
                
                if cn_value != 2:
                    
                    # Look up k: number of controls with the same CN value
                    k = control_counts_df.loc[gene, cn_value] if cn_value in control_counts_df.columns else 0
                    
                    # Construct 2x2 table cells for a stable Odds Ratio calculation
                    a, b = 1.5, 0.5
                    c, d = k + 0.5, (N - k) + 0.5
                    
                    # Calculate the corrected odds ratio
                    or_corrected = (a * d) / (b * c)
                    
                    # Handle edge case for log transform if OR is somehow non-positive
                    if or_corrected <= 0:
                        or_corrected = epsilon
                        
                    # Look up the pre-computed standard deviation for the gene
                    sigma_for_gene = sigma_controls.loc[gene]
                    
                    # Calculate the final score using the enhanced formula
                    score = (np.sign(cn_value - 2) * np.log10(or_corrected)) / (sigma_for_gene + epsilon)

                else: # cn_value == 2
                    
                    # Look up pre-computed mean and std dev for the gene
                    mu_for_gene = mu_controls.loc[gene]
                    sigma_for_gene = sigma_controls.loc[gene]
                    
                    # Calculate the Z-score relative to the control mean to capture nuance
                    score = (2 - mu_for_gene) / (sigma_for_gene + epsilon)
                
                scores.append(score)

            df_sample = pd.DataFrame(
                {'adjusted_weight': scores},
                index=df.index
            )
            ranked_dfs[sample_name] = df_sample
            
        return ranked_dfs
    
    else:
        raise ValueError("Omic type must be 'RNA', 'DNAm', or 'CNV'")