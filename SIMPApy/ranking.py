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

            # Delete the DataFrame to free up memory (optional)
            del sample_df

        return ranked_dfs
    
    elif omic.upper() == "CNV":
        # Calculate baseline importance (Bg) for each gene
        control_data = df.filter(regex='^tw')  # Control data
        control_std = control_data.std(axis=1)
        control_max = control_data.max(axis=1)
        Bg = control_std / control_max

        # Create a dictionary to store results
        adjusted_weights_dict = {}

        # Loop through all samples (cases and controls)
        for sample_name in df.columns:
            x_s = df[sample_name]  # Copy numbers for the current sample

            # Calculate non-linear weight (w(x_s,g))
            w_x_s_g = 2 ** np.abs(x_s - 2) * np.sign(x_s - 2)

            # Calculate adjusted weight (w_adjusted(x_s,g))
            adjusted_weight = w_x_s_g.where(w_x_s_g != 0, Bg)

            # Create a DataFrame for the current sample
            df_sample = pd.DataFrame({
                'adjusted_weight': adjusted_weight
            })

            adjusted_weights_dict[sample_name] = df_sample

        return adjusted_weights_dict
    
    else:
        raise ValueError("Omic type must be 'RNA', 'DNAm', or 'CNV'")