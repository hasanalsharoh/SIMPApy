import numpy as np
import pandas as pd


def _extract_tag_genes(df, group):
    """
    Extracts tag genes and NES values for a given group ("tm" or "tw"),
    unnests the tag gene lists into separate rows, and adds a column indicating
    whether the tag gene is present in the matched_genes list.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing enrichment results
    group : str
        Group identifier ("tm" or "tw")
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with tag genes information
    """
    matched_genes_list = []
    for index, row in df.iterrows():
        sample = row['sample_name']
        term = row['Term']
        mpes = row['mpes']
        fdr = row['fdr_bh']

        for omic in ['rna', 'cnv', 'dna']:
            tag_genes = row[f'{omic}_lead_genes']
            nes = row[f'{omic}_nes']
            matched_genes_str = row['matched_genes']
            if isinstance(matched_genes_str, str):  # Check if matched_genes is a string
                matched_genes = matched_genes_str.split(';')

                for gene in matched_genes:
                    # Determine the corresponding tag_gene if available
                    tag_gene = np.nan
                    if isinstance(tag_genes, str):
                        tag_genes_list = tag_genes.split(';')
                        if gene in tag_genes_list:
                            tag_gene = gene

                    matched_genes_list.append({
                        'sample_name': sample,
                        'Term': term,
                        'omic': omic,
                        'matched_gene': gene,
                        'nes': nes,
                        'fdr': fdr,
                        'mpes': mpes,
                        'group': group,
                        'tag_gene': tag_gene  # Add tag_gene information
                    })

    return pd.DataFrame(matched_genes_list)


def _create_aggregated_dataframes(melted_raw_data):
    """
    Creates aggregated dataframes (tmas and twas) with omics values, NES, FDR, and tag information.
    
    Parameters
    ----------
    melted_raw_data : pandas.DataFrame
        Melted DataFrame with raw omics data
        
    Returns
    -------
    tuple
        Two DataFrames (tmas, twas) with aggregated data
    """
    # Separate into 'tm' and 'tw' groups
    tm_df = melted_raw_data[melted_raw_data['group'] == 'tm']
    tw_df = melted_raw_data[melted_raw_data['group'] == 'tw']

    def aggregate_group(df):
        """Aggregates data for a single group ('tm' or 'tw')."""
        # Create an empty list to store the aggregated data
        aggregated_data = []

        # Get unique combinations of sample_name and Term
        unique_combinations = df[['sample_name', 'Term']].drop_duplicates()

        for index, row in unique_combinations.iterrows():
            sample_name = row['sample_name']
            term = row['Term']

            # Filter data for the current sample_name and Term
            filtered_df = df[(df['sample_name'] == sample_name) & (df['Term'] == term)]
            # --- Get the FDR value for this sample_name and Term (FDR is not omic-specific) ---
            fdr_value = filtered_df['fdr'].iloc[0] if not filtered_df.empty else np.nan  # All rows have same fdr

            # Get unique gene names for this combination
            unique_genes = filtered_df['gene_name'].unique()

            for gene in unique_genes:
                gene_data = {'sample_name': sample_name, 'Term': term, 'gene_name': gene, 'fdr': fdr_value}

                for omic in ['rna', 'cnv', 'dna']:
                    omic_data = filtered_df[(filtered_df['omic'] == omic) & (filtered_df['gene_name'] == gene)]

                    if not omic_data.empty:
                        gene_data[f'{omic}_value'] = omic_data['raw_value'].iloc[0]
                        gene_data[f'{omic}_nes'] = omic_data['nes'].iloc[0]
                        gene_data[f'{omic}_tag'] = omic_data['tag_gene'].iloc[0] if not pd.isna(omic_data['tag_gene'].iloc[0]) else ""
                        if 'mpes' not in gene_data:  # Take mpes from the first non-empty omic data (assuming it's consistent)
                          gene_data['mpes'] = omic_data['mpes'].iloc[0]
                    else:
                        gene_data[f'{omic}_value'] = np.nan
                        gene_data[f'{omic}_nes'] = np.nan
                        gene_data[f'{omic}_tag'] = ""

                aggregated_data.append(gene_data)

        return pd.DataFrame(aggregated_data)

    # Aggregate for each group
    tmas = aggregate_group(tm_df)
    twas = aggregate_group(tw_df)

    return tmas, twas


def process_multiomics_data(sigimgsea, tp, cn, meth, pop_info=None):
    """
    Process multi-omics data for visualization, combining enrichment results
    with raw data matrices.
    
    Parameters
    ----------
    sigimgsea : pandas.DataFrame
        DataFrame containing enrichment results
    tp : pandas.DataFrame
        Gene expression (TPM) data with genes as rows and samples as columns
    cn : pandas.DataFrame
        Copy number variation data with genes as rows and samples as columns
    meth : pandas.DataFrame
        DNA methylation data with genes as rows and samples as columns
    pop_info : pandas.DataFrame, optional
        Sample metadata with cancer types and clinical information
        
    Returns
    -------
    tuple
        Two DataFrames (tmas, twas) with processed data ready for visualization
    """
    # Separate sigimgsea data based on sample name prefix
    tmas_sig = sigimgsea[sigimgsea['sample_name'].str.contains('tm')]
    twas_sig = sigimgsea[sigimgsea['sample_name'].str.contains('tw')]
    
    # Extract tag genes for each group
    tmas_genes = _extract_tag_genes(tmas_sig, 'tm')
    twas_genes = _extract_tag_genes(twas_sig, 'tw')
    tag_genes = pd.concat([tmas_genes, twas_genes])
    
    # Get unique genes and samples from tag genes
    genes = list(tag_genes['matched_gene'].unique())
    samples = list(tag_genes['sample_name'].unique())
    
    # Filter and align raw data matrices
    tp_filtered = tp[samples].loc[genes]
    cn_filtered = cn[samples].loc[genes]
    meth_filtered = meth[samples].loc[genes]
    
    # Reset indices for each DataFrame
    tp_reset = tp_filtered.reset_index()
    cn_reset = cn_filtered.reset_index()
    meth_reset = meth_filtered.reset_index().rename(columns={'index': 'gene_name'})
    
    # Add omic column to each DataFrame
    tp_reset['omic'] = 'rna'
    cn_reset['omic'] = 'cnv'
    meth_reset['omic'] = 'dna'
    
    # Merge the DataFrames based on 'gene_name' and sample columns
    merged_raw_data = pd.concat([tp_reset, cn_reset, meth_reset], ignore_index=True)
    
    # Melt the merged DataFrame
    melted_raw_data = merged_raw_data.melt(
        id_vars=['gene_name', 'omic'],
        var_name='sample_name',
        value_name='raw_value'
    )
    
    # Set group based on sample name prefix
    tag_genes['group'] = np.where(tag_genes['sample_name'].str.startswith('tm'), 'tm', 'tw')
    melted_raw_data['group'] = np.where(melted_raw_data['sample_name'].str.startswith('tm'), 'tm', 'tw')
    
    # Merge with tag_genes to associate raw values with enrichment results
    melted_raw_data = pd.merge(
        melted_raw_data,
        tag_genes,
        left_on=['sample_name', 'gene_name', 'omic', 'group'],
        right_on=['sample_name', 'matched_gene', 'omic', 'group'],
        how='left'
    )
    
    # Create aggregated dataframes
    tmas, twas = _create_aggregated_dataframes(melted_raw_data)
    
    # Convert tag indicators to binary values
    for df in [tmas, twas]:
        df['rna_tag'] = np.where(df['rna_tag'] == '', 0, 1)
        df['cnv_tag'] = np.where(df['cnv_tag'] == '', 0, 1)
        df['dna_tag'] = np.where(df['dna_tag'] == '', 0, 1)
        
        # Clip TPM values in rna_value columns to 1000
        df['rna_value'] = df['rna_value'].clip(upper=1000)
    
    # Add population info if provided
    if pop_info is not None:
        # Try to add cancer type and staging info, handling potential missing columns gracefully
        required_cols = ['key']
        optional_cols = ['cancer_type', 'ajcc_pathologic_stage']
        
        # Determine which optional columns are available
        available_cols = ['key'] + [col for col in optional_cols if col in pop_info.columns]
        
        # Merge with available population info
        tmas = pd.merge(tmas, pop_info[available_cols], left_on='sample_name', right_on='key', how='left')
        twas = pd.merge(twas, pop_info[available_cols], left_on='sample_name', right_on='key', how='left')
        
        # Drop key column if it was added during merge
        if 'key' in tmas.columns:
            tmas.drop(columns=['key'], inplace=True)
        if 'key' in twas.columns:
            twas.drop(columns=['key'], inplace=True)
    
    # Ensure cancer_type and ajcc_pathologic_stage columns exist even if pop_info was not provided
    if 'cancer_type' not in tmas.columns:
        tmas['cancer_type'] = "Not Available"
        twas['cancer_type'] = "Not Available"
    
    if 'ajcc_pathologic_stage' not in tmas.columns:
        tmas['ajcc_pathologic_stage'] = "Not Available"
        twas['ajcc_pathologic_stage'] = "Not Available"
    
    return tmas, twas