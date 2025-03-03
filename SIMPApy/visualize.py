import numpy as np
import plotly.graph_objects as go
import pandas as pd
import traceback
from ipywidgets import Dropdown, HBox, VBox, Text, Button, Checkbox, RadioButtons, Output
from IPython.display import display


def _create_traces(data, selected_term=None, selected_omic='rna', selected_sample=None, selected_gene=None,
                  selected_cancer_type=None, normalize_untagged=False, selected_ajcc_stage=None):
    """
    Create plotly trace for the 3D scatter plot visualization.
    
    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame with processed multi-omics data
    selected_term : str, optional
        Term to filter data by
    selected_omic : str, optional
        Omic type to use for coloring ('rna', 'cnv', 'dna', or 'mpes')
    selected_sample : str, optional
        Sample to filter data by
    selected_gene : str, optional
        Gene to filter data by
    selected_cancer_type : str, optional
        Cancer type to filter data by
    normalize_untagged : bool, optional
        Whether to normalize values for untagged genes
    selected_ajcc_stage : str, optional
        AJCC pathologic stage to filter data by
        
    Returns
    -------
    tuple
        Plotly trace object and filtered DataFrame
    """
    filtered_data = data.copy()

    if selected_term:
        filtered_data = filtered_data[filtered_data['Term'] == selected_term]
    if selected_sample:
        filtered_data = filtered_data[filtered_data['sample_name'] == selected_sample]
    if selected_gene:
        filtered_data = filtered_data[filtered_data['gene_name'] == selected_gene]
    if selected_cancer_type and selected_cancer_type != "All":
        filtered_data = filtered_data[filtered_data['cancer_type'] == selected_cancer_type]
    # --- AJCC Stage Filtering ---
    if selected_ajcc_stage and selected_ajcc_stage != "All Stages":
        filtered_data = filtered_data[filtered_data['ajcc_pathologic_stage'] == selected_ajcc_stage]

    if selected_omic == 'mpes':
        nes_col = 'mpes'
        filtered_data['all_tagged'] = (filtered_data['rna_tag'] != 0) & (filtered_data['cnv_tag'] != 0) & (filtered_data['dna_tag'] != 0)
        filtered_data['all_tagged'] = filtered_data['all_tagged'].astype(int)
        tag_col = 'all_tagged'
    else:
        nes_col = f'{selected_omic}_nes'
        tag_col = f'{selected_omic}_tag'

    if normalize_untagged:
        filtered_data.loc[filtered_data[tag_col] == 0, nes_col] = 0

    filtered_data['fdr'] = filtered_data['fdr'].replace([np.inf, -np.inf], np.nan).round(5)
    filtered_data['ajcc_pathologic_stage'] = filtered_data['ajcc_pathologic_stage'].fillna("Not Available")

    trace = go.Scatter3d(
        x=filtered_data['dna_value'],
        y=filtered_data['rna_value'],
        z=filtered_data['cnv_value'],
        mode='markers',
        marker=dict(
            size=8,
            opacity=0.8,
            color=filtered_data[nes_col],
            colorscale="Bluered",
            colorbar=dict(title=f"{selected_omic.upper()}"),
            showscale=True
        ),
        text=filtered_data.apply(lambda x: f"{x['gene_name']}_{x['sample_name']}", axis=1),
        hovertemplate=(
            "<b>Gene:</b> %{customdata[0]}<br>"
            "<b>Sample:</b> %{customdata[1]}<br>"
            "<b>Term:</b> %{customdata[2]}<br>"
            "<b>Cancer Type:</b> %{customdata[3]}<br>"
            "<b>AJCC Stage:</b> %{customdata[7]}<br>"            
            "<b>DNAm:</b> %{x:.3f}<br>"
            "<b>TPM:</b> %{y:.1f}<br>"
            "<b>CNV:</b> %{z:.1f}<br>"
            "<b>Enrichment:</b> %{customdata[4]:.3f}<br>"
            "<b>Tagged:</b> %{customdata[5]}<br>"
            "<b>Multiomics FDR:</b> %{customdata[6]:.5f}<br>"
            "<extra></extra>"
        ),
        customdata=filtered_data[['gene_name', 'sample_name', 'Term', 'cancer_type', nes_col, tag_col, 'fdr', 'ajcc_pathologic_stage']].values,
        name=selected_omic.upper()
    )

    return trace, filtered_data


def create_interactive_plot(data, title_suffix=""):
    """
    Create an interactive 3D scatter plot for visualizing multi-omics data.
    
    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame with processed multi-omics data
    title_suffix : str, optional
        Suffix to add to the plot title
        
    Returns
    -------
    None
        Displays the interactive plot in the notebook
    """
    class PlotState:
        def __init__(self):
            self.current_term = None
            self.current_omic = 'rna'
            self.current_cancer_type = "All"
            self.current_sample = None
            self.current_gene = None
            self.normalize_untagged = False
            self.filter_type = 'Sample'
            self.color_min = None
            self.color_max = None
            self.current_ajcc_stage = "All Stages"  # Initialize AJCC Stage

    state = PlotState()

    unique_terms = sorted(data['Term'].unique())
    initial_term = unique_terms[0]
    state.current_term = initial_term

    unique_cancer_types = ["All"] + sorted(data['cancer_type'].unique())

    # --- Create unique AJCC Stages list (including "Not Available") ---
    # Handle NaN before sorting
    data['ajcc_pathologic_stage'] = data['ajcc_pathologic_stage'].fillna("Not Available")
    unique_ajcc_stages = ["All Stages"] + sorted(data['ajcc_pathologic_stage'].unique(), 
                                                key=lambda x: (x == "Not Available", x))

    nes_columns_for_range = ['rna_nes', 'cnv_nes', 'dna_nes']
    state.color_min = min(data[nes_columns_for_range].min())
    state.color_max = max(data[nes_columns_for_range].max())

    fig = go.FigureWidget()

    initial_trace, _ = _create_traces(
        data,
        selected_term=initial_term,
        selected_omic='rna',
        selected_ajcc_stage="All Stages"  # Initial AJCC Stage
    )
    fig.add_trace(initial_trace)

    omic_dropdown = Dropdown(
        options=['rna', 'cnv', 'dna', 'mpes'],
        value='rna',
        description='Omic Type:',
        disabled=False,
    )

    term_dropdown = Dropdown(
        options=unique_terms,
        value=initial_term,
        description='Term:',
        disabled=False,
    )

    cancer_type_dropdown = Dropdown(
        options=unique_cancer_types,
        value="All",
        description='Cancer Type:',
        disabled=False,
    )

    # --- AJCC Stage Dropdown ---
    ajcc_stage_dropdown = Dropdown(
        options=unique_ajcc_stages,
        value="All Stages",
        description='AJCC Stage:',
        disabled=False,
    )

    filter_type = RadioButtons(
        options=['Sample', 'Gene'],
        value='Sample',
        description='Filter by:',
        disabled=False
    )

    search_box = Text(
        value='',
        placeholder='Type sample/gene name...',
        description='Search:',
        disabled=False
    )

    search_button = Button(description="Search")

    normalize_checkbox = Checkbox(
        value=False,
        description='Normalize Untagged',
        disabled=False
    )

    reset_button = Button(description="Reset View")
    export_button = Button(description="Export (300 DPI)")
    output_widget = Output()

    def export_high_res(b):
        with output_widget:
            try:
                print("Starting export process...")
                width_px = 4 * 300
                height_px = 4 * 300
                print("Creating static figure...")
                static_fig = go.Figure(data=fig.data, layout=fig.layout)
                print("Updating layout...")
                static_fig.update_layout(
                    width=width_px,
                    height=height_px
                )
                print("Saving file...")
                import os
                # Create a more flexible file path approach
                desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
                file_path = os.path.join(desktop_path, f"{title_suffix}.png")
                static_fig.write_image(file_path, scale=1, width=width_px, height=height_px)
                print(f"Successfully exported high-resolution image as '{file_path}'")
            except Exception as e:
                print(f"Error during export: {str(e)}")
                print("\nFull traceback:")
                traceback.print_exc()
                print("\nPlease ensure you have installed required packages:")
                print("pip install kaleido plotly -U")

    def update_plot():
        new_trace, filtered_data = _create_traces(
            data,
            selected_term=state.current_term,
            selected_omic=state.current_omic,
            selected_sample=state.current_sample,
            selected_gene=state.current_gene,
            selected_cancer_type=state.current_cancer_type,
            normalize_untagged=state.normalize_untagged,
            selected_ajcc_stage=state.current_ajcc_stage  # Pass AJCC Stage
        )

        if state.current_omic == 'mpes':
            new_trace.marker.update(
                colorscale="Bluered",
                cmin=-5,
                cmax=5,
                showscale=True,
                colorbar=dict(title=f"{state.current_omic.upper()}")
            )
        else:
            new_trace.marker.update(
                colorscale="Bluered",
                cmin=state.color_min,
                cmax=state.color_max,
                showscale=True,
                colorbar=dict(title=f"{state.current_omic.upper()}")
            )

        with fig.batch_update():
            fig.data[0].x = new_trace.x
            fig.data[0].y = new_trace.y
            fig.data[0].z = new_trace.z
            fig.data[0].marker = new_trace.marker
            fig.data[0].customdata = new_trace.customdata
            fig.data[0].hovertemplate = new_trace.hovertemplate

    def on_term_change(change):
        state.current_term = change['new']
        state.current_sample = None
        state.current_gene = None
        update_plot()

    def on_omic_change(change):
        state.current_omic = change['new']
        update_plot()

    def on_cancer_type_change(change):
        state.current_cancer_type = change['new']
        state.current_sample = None
        state.current_gene = None
        update_plot()

    # --- AJCC Stage Change Handler ---
    def on_ajcc_stage_change(change):
        state.current_ajcc_stage = change['new']
        update_plot()

    def on_normalize_change(change):
        state.normalize_untagged = change['new']
        update_plot()

    def on_filter_type_change(change):
        state.filter_type = change['new']
        state.current_sample = None
        state.current_gene = None
        search_box.value = ''
        update_plot()

    def on_search_click(b):
        search_term = search_box.value.strip()
        if search_term:
            if state.filter_type == 'Sample':
                samples = data['sample_name'].unique()
                matching_samples = [s for s in samples if search_term.lower() in s.lower()]
                if matching_samples:
                    state.current_sample = matching_samples[0]
                    state.current_gene = None
            else:
                genes = data['gene_name'].unique()
                matching_genes = [g for g in genes if search_term.lower() in g.lower()]
                if matching_genes:
                    state.current_gene = matching_genes[0]
                    state.current_sample = None
            update_plot()

    def handle_click(trace, points, selector):
        if points.point_inds:
            selected_data = trace.customdata[points.point_inds[0]]
            if state.filter_type == 'Gene':
                state.current_gene = selected_data[0]
                state.current_sample = None
            else:
                state.current_sample = selected_data[1]
                state.current_gene = None
            update_plot()

    def reset_plot(b):
        state.current_sample = None
        state.current_gene = None
        search_box.value = ''
        if state.filter_type != 'Sample':
            state.filter_type = 'Sample'
            filter_type.value = 'Sample'
        state.current_ajcc_stage = "All Stages"  # Reset AJCC stage
        ajcc_stage_dropdown.value = "All Stages"
        update_plot()

    term_dropdown.observe(on_term_change, names='value')
    omic_dropdown.observe(on_omic_change, names='value')
    cancer_type_dropdown.observe(on_cancer_type_change, names='value')
    # --- Observe AJCC Stage Dropdown ---
    ajcc_stage_dropdown.observe(on_ajcc_stage_change, names='value')
    normalize_checkbox.observe(on_normalize_change, names='value')
    filter_type.observe(on_filter_type_change, names='value')
    search_button.on_click(on_search_click)
    reset_button.on_click(reset_plot)
    fig.data[0].on_click(handle_click)
    export_button.on_click(export_high_res)

    fig.update_layout(
        title=f"{title_suffix}",
        title_x=0.5,
        title_y=0.95,
        scene=dict(
            xaxis=dict(
                title='DNAm',
                backgroundcolor="white",
                gridcolor="grey"
            ),
            yaxis=dict(
                title='TPM',
                backgroundcolor="white",
                gridcolor="grey"
            ),
            zaxis=dict(
                title='CNV',
                backgroundcolor="white",
                gridcolor="grey"
            ),
            bgcolor="white"
        ),
        width=800,
        height=700,
    )
    controls_row1 = HBox([term_dropdown, omic_dropdown, cancer_type_dropdown, ajcc_stage_dropdown])  # Added to row 1
    controls_row2 = HBox([filter_type, search_box, search_button, normalize_checkbox, reset_button, export_button])
    display(VBox([controls_row1, controls_row2, fig, output_widget]))