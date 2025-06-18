"""
Analysis module for the SIMPApy package.

This module provides functions for statistical analysis and visualization of predefined groups.
"""
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

def group_diffs(
    data: pd.DataFrame,
    pathway_col: str,
    value_col: str,
    group_col: str,
    group1_prefix: str = 'tm',
    group2_prefix: str = 'tw',
    adj_method: str = 'fdr_bh',
    dropna: bool = True
) -> pd.DataFrame:
    """
    Performs Mann-Whitney U test for group differences for each pathway.

    Args:
        data (pd.DataFrame): DataFrame containing the pathway analysis results.
                             Must contain columns for pathway names, values, and sample groups.
        pathway_col (str): The name of the column containing pathway/term names.
        value_col (str): The name of the column with the numerical values to compare (e.g., 'mpes', 'combined_z').
        group_col (str): The name of the column that defines the sample groups (e.g., 'sample_name').
        group1_prefix (str): The prefix or identifier for the first group within the group_col.
        group2_prefix (str): The prefix or identifier for the second group within the group_col.
        adj_method (str): Method for p-value correction using statsmodels.multitest.
                          Defaults to 'fdr_bh'.
        dropna (bool): If True, drops NaN values from group values before testing. 
                       Note: mannwhitneyu may not handle NaNs, so setting to False could cause errors.
                       Defaults to True.

    Returns:
        pd.DataFrame: A DataFrame with group differences for each pathway, including
                      mean difference, p-value, adjusted p-value, and -log10(p_adj).
    """
    pathway_results = []
    unique_pathways = data[pathway_col].unique()

    group1_samples = data[data[group_col].str.startswith(group1_prefix)]
    group2_samples = data[data[group_col].str.startswith(group2_prefix)]

    for pathway in unique_pathways:
        pathway_data = data[data[pathway_col] == pathway]
        
        group1_values = pathway_data[pathway_data[group_col].isin(group1_samples[group_col])][value_col]
        group2_values = pathway_data[pathway_data[group_col].isin(group2_samples[group_col])][value_col]

        if dropna:
            group1_values = group1_values.dropna()
            group2_values = group2_values.dropna()

        if len(group1_values) > 1 and len(group2_values) > 1:
            mean_diff = group1_values.mean() - group2_values.mean()
            stat, p_value = mannwhitneyu(group1_values, group2_values, alternative='two-sided')
            
            pathway_results.append({
                'pathway': pathway,
                'mean_diff': mean_diff,
                'p_value': p_value,
                'statistic': stat
            })

    if not pathway_results:
        return pd.DataFrame()

    results_df = pd.DataFrame(pathway_results)
    
    # P-value adjustment
    if 'p_value' in results_df.columns:
        p_values = results_df['p_value'].dropna()
        if not p_values.empty:
            _, p_adj, _, _ = multipletests(p_values, method=adj_method)
            results_df.loc[p_values.index, 'p_adj'] = p_adj
            results_df['neg_log10_p_adj'] = -np.log10(results_df['p_adj'])

    return results_df.sort_values('p_value').reset_index(drop=True)


def plot_volcano(
    data: pd.DataFrame,
    x_col: str = 'mean_diff',
    y_col: str = 'neg_log10_p_adj',
    p_thresh: float = 0.05,
    title: str = "Group Differences Volcano Plot",
    xlabel: str = "Mean Difference",
    ylabel: str = "-log10(Adjusted P-value)",
    save_path: str = None
) -> None:
    """
    Generates and optionally saves a volcano plot from group difference results.

    Args:
        data (pd.DataFrame): DataFrame from group_diffs function.
        x_col (str): Column for the x-axis (mean difference).
        y_col (str): Column for the y-axis (-log10 adjusted p-value).
        p_thresh (float): Significance threshold for p-value.
        title (str): The title of the plot.
        xlabel (str): The label for the x-axis.
        ylabel (str): The label for the y-axis.
        save_path (str): If provided, the plot is saved to this path.
    """
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=data, x=x_col, y=y_col, color='blue')

    neg_log10_p_threshold = -np.log10(p_thresh)
    plt.axhline(y=neg_log10_p_threshold, color='gray', linestyle='--')

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)

    if save_path:
        plt.savefig(save_path, dpi=300)
    
    plt.show()

def calculate_correlation(data: pd.DataFrame, x_col: str, y_col: str, group_col: str = None) -> pd.DataFrame:
    """
    Calculates Spearman correlation between two columns, optionally grouped by a third column.

    Args:
        data: DataFrame containing the data.
        x_col: The name of the first column for correlation.
        y_col: The name of the second column for correlation.
        group_col: If provided, correlations are calculated for each group defined by this column.

    Returns:
        A DataFrame with the correlation results.
    """
    correlation_results = []

    if group_col:
        for group_name, group_df in data.groupby(group_col):
            corr, p_val = spearmanr(group_df[x_col], group_df[y_col])
            correlation_results.append({
                'group': group_name,
                'correlation': corr,
                'p_value': p_val
            })
    else:
        corr, p_val = spearmanr(data[x_col], data[y_col])
        correlation_results.append({
            'group': 'overall',
            'correlation': corr,
            'p_value': p_val
        })

    return pd.DataFrame(correlation_results)


def plot_correlation_scatterplot(
    data: pd.DataFrame, 
    x_col: str, 
    y_col: str, 
    hue_col: str = None, 
    palette: dict = None,
    title: str = "Correlation Scatter Plot",
    xlabel: str = "X Value",
    ylabel: str = "Y Value",
    save_path: str = None
) -> None:
    """
    Generates and optionally saves a scatter plot with a regression line.

    Args:
        data: DataFrame containing the data.
        x_col: The name of the column for the x-axis.
        y_col: The name of the column for the y-axis.
        hue_col: If provided, points are colored based on this column.
        palette: A dictionary to map hue values to colors.
        title: The title of the plot.
        xlabel: The label for the x-axis.
        ylabel: The label for the y-axis.
        save_path: If provided, the plot is saved to this path.
    """
    plt.figure(figsize=(8, 8))
    sns.lmplot(
        data=data,
        x=x_col,
        y=y_col,
        hue=hue_col,
        palette=palette,
        scatter_kws={'alpha':0.5}
    )
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)

    if save_path:
        plt.savefig(save_path, dpi=300)
        
    plt.show()