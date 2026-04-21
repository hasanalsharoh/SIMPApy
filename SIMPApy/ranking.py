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


def _calculate_msd_robust(
    df: pd.DataFrame,
    alpha: float = 0.05,
    asymmetric: bool = True,
    kappa: float = 1.4826,
) -> pd.DataFrame:
    """
    Robust Minimum Significant Deviation via median and (optionally directional) MAD.

    Returns a DataFrame indexed by gene with columns:
      - 'center'  : gene-wise median of the TWA group
      - 'msd_up'  : MSD applied when D_{x,s} >= 0
      - 'msd_dn'  : MSD applied when D_{x,s} <  0

    When asymmetric=False, msd_up == msd_dn (classical symmetric robust MSD).

    Fallback hierarchy for zero-scale genes:
      1. IQR / 1.349  (Gaussian-consistent, nonzero unless >=25% ties at median)
      2. global median of non-zero robust scales across genes
    """
    twa_cols = [c for c in df.columns if c.startswith('tw')]
    if len(twa_cols) < 3:
        raise ValueError("Robust MSD requires >= 3 TWA samples.")
    twa = df[twa_cols]

    # Robust center and deviations
    center = twa.median(axis=1)
    devs = twa.subtract(center, axis=0)

    # Symmetric MAD (used as a sanity fallback inside asymmetric branch)
    mad_sym = devs.abs().median(axis=1)
    sigma_sym = kappa * mad_sym

    # IQR-based and global fallbacks
    iqr = twa.quantile(0.75, axis=1) - twa.quantile(0.25, axis=1)
    sigma_iqr = iqr / 1.349
    pool = sigma_sym[sigma_sym > 0]
    global_sigma = float(pool.median()) if not pool.empty else 1e-6
    global_sigma = max(global_sigma, 1e-6)

    def _apply_fallbacks(sigma: pd.Series) -> pd.Series:
        s = sigma.copy().astype(float)
        bad = (s <= 0) | ~np.isfinite(s)
        s.loc[bad] = sigma_iqr.loc[bad]
        bad = (s <= 0) | ~np.isfinite(s)
        s.loc[bad] = global_sigma
        return s

    if asymmetric:
        pos = devs.where(devs > 0)            # NaN where dev <= 0
        neg = (-devs).where(devs < 0)         # positive magnitudes, NaN elsewhere
        mad_up = pos.median(axis=1)
        mad_dn = neg.median(axis=1)
        sigma_up = kappa * mad_up
        sigma_dn = kappa * mad_dn
        # If a tail is empty/zero, fall back to symmetric MAD *before* IQR
        sigma_up = sigma_up.where(sigma_up > 0, sigma_sym)
        sigma_dn = sigma_dn.where(sigma_dn > 0, sigma_sym)
    else:
        sigma_up = sigma_sym.copy()
        sigma_dn = sigma_sym.copy()

    sigma_up = _apply_fallbacks(sigma_up)
    sigma_dn = _apply_fallbacks(sigma_dn)

    z = norm.ppf(1 - alpha / 2)
    return pd.DataFrame({
        'center': center,
        'msd_up': z * sigma_up,
        'msd_dn': z * sigma_dn,
    })

def calculate_ranking(
    df: pd.DataFrame,
    omic: str = "RNA",
    alpha: float = 0.05,
    robust: bool = True,
    asymmetric: bool = True,
):
    """
    Parameters
    omic : str
        Type of omic data: "RNA", "DNAm", or "CNV".
    alpha : float
        Significance level for ranking (used in MSD calculation).
    robust : bool
        If True, use median + MAD-based MSD (non-parametric).
    asymmetric : bool
        If True (and robust=True), use directional (double) MAD
        to handle skewed baseline distributions.
        (asymetric and robust are highly recommended together and are the default)
        Using assymetric and robust together is preferred for non parametric data, and is the more conservative approach.
        Having both false is the classical parametric approach which may be more powerful if assumptions are met, but is more sensitive to outliers and non-normality.
    """
    if omic.upper() in ["RNA", "DNAM"]:
        twa_cols = [c for c in df.columns if c.startswith('tw')]

        if robust:
            tbl = _calculate_msd_robust(df, alpha=alpha, asymmetric=asymmetric)
            center = tbl['center']
            msd_up, msd_dn = tbl['msd_up'], tbl['msd_dn']
        else:
            center = df[twa_cols].mean(axis=1)
            msd_par = _calculate_msd(df, alpha)
            msd_up = msd_par
            msd_dn = msd_par

        ranked_dfs = {}
        for sample in df.columns:
            d_xs = df[sample] - center
            # Directional scale: msd_up for positive deviations, msd_dn for negative
            scale = np.where(d_xs >= 0, msd_up.values, msd_dn.values)
            # Avoid zero division from any edge case
            scale = np.where(scale > 0, scale, np.nan)
            weighted = d_xs.values / scale
            msd_signed = np.where(d_xs >= 0, msd_up.values, -msd_dn.values)

            ranked_dfs[sample] = pd.DataFrame({
                'D_xs': d_xs.values,
                'MSD': msd_signed,
                'weighted': weighted,
                'Significant': np.abs(d_xs.values) > np.where(d_xs >= 0, msd_up.values, msd_dn.values),
                'Rank': pd.Series(d_xs.values, index=df.index).rank(ascending=False).values,
            }, index=df.index)
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