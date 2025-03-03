import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np
from SIMPApy.visualize import create_interactive_plot, _create_traces

import plotly.graph_objects as go

# Import the function from parent directory


def create_test_dataframe():
    """Create a sample dataframe for testing visualization functions"""
    data = {
        'Term': ['Term1', 'Term1', 'Term2', 'Term2'],
        'sample_name': ['Sample1', 'Sample2', 'Sample1', 'Sample2'],
        'gene_name': ['Gene1', 'Gene2', 'Gene1', 'Gene2'],
        'cancer_type': ['Type1', 'Type1', 'Type2', 'Type2'],
        'ajcc_pathologic_stage': ['Stage1', 'Stage2', 'Stage1', np.nan],
        'dna_value': [0.1, 0.2, 0.3, 0.4],
        'rna_value': [10.5, 20.5, 30.5, 40.5],
        'cnv_value': [1.1, 1.2, 1.3, 1.4],
        'rna_nes': [2.0, -2.0, 3.0, -3.0],
        'cnv_nes': [1.5, -1.5, 2.5, -2.5],
        'dna_nes': [1.0, -1.0, 2.0, -2.0],
        'rna_tag': [1, 0, 1, 0],
        'cnv_tag': [1, 0, 1, 0],
        'dna_tag': [1, 0, 1, 0],
        'fdr': [0.01, 0.05, 0.02, np.inf],
        'mpes': [0.5, -0.5, 1.5, -1.5]
    }
    return pd.DataFrame(data)


class TestVisualizeModule(unittest.TestCase):
    
    def setUp(self):
        self.test_data = create_test_dataframe()
    
    def test_create_traces_basic(self):
        """Test that _create_traces produces the expected trace object"""
        trace, filtered_data = _create_traces(self.test_data)
        
        # Check return types
        self.assertIsInstance(trace, go.Scatter3d)
        self.assertIsInstance(filtered_data, pd.DataFrame)
        self.assertEqual(filtered_data.shape[0], self.test_data.shape[0])
    
    def test_create_traces_filtering(self):
        """Test filtering capabilities of _create_traces"""
        # Test filtering by term
        trace, filtered_data = _create_traces(self.test_data, selected_term='Term1')
        self.assertEqual(filtered_data.shape[0], 2)
        self.assertTrue(all(filtered_data['Term'] == 'Term1'))
        
        # Test filtering by AJCC stage
        trace, filtered_data = _create_traces(self.test_data, selected_ajcc_stage='Stage1')
        self.assertEqual(filtered_data.shape[0], 2)
        self.assertTrue(all(filtered_data['ajcc_pathologic_stage'] == 'Stage1'))
    
    @patch('SIMPApy.visualize.display')
    @patch('SIMPApy.visualize.VBox')
    @patch('SIMPApy.visualize.HBox') 
    @patch('SIMPApy.visualize.go.FigureWidget')
    def test_create_interactive_plot_basic(self, mock_figure_widget, mock_hbox, mock_vbox, mock_display):
        """Basic smoke test for create_interactive_plot"""
        # Setup mock figure
        mock_fig = MagicMock()
        mock_data = MagicMock()
        mock_data.on_click = MagicMock()
        mock_fig.data = [mock_data]
        mock_fig.add_trace = MagicMock()
        mock_fig.update_layout = MagicMock()
        mock_figure_widget.return_value = mock_fig
        
        # Mock VBox and HBox to return MagicMocks
        mock_vbox.return_value = MagicMock()
        mock_hbox.return_value = MagicMock()
        
        # Call the function
        create_interactive_plot(self.test_data, title_suffix="Test")
        
        # Verify basic functionality
        mock_figure_widget.assert_called_once()
        mock_fig.add_trace.assert_called_once()
        mock_fig.update_layout.assert_called_once()
        mock_display.assert_called_once()
    
    @patch('SIMPApy.visualize.display')
    @patch('SIMPApy.visualize.VBox')
    @patch('SIMPApy.visualize.HBox')
    @patch('SIMPApy.visualize.go.FigureWidget')
    @patch('SIMPApy.visualize._create_traces')
    def test_create_interactive_plot_widgets_integration(self, mock_create_traces, 
                                                      mock_figure_widget, mock_hbox,
                                                      mock_vbox, mock_display):
        """Test that widgets are properly integrated with trace creation"""
        # Setup mocks
        mock_trace = MagicMock()
        mock_filtered_data = self.test_data.copy()
        mock_create_traces.return_value = (mock_trace, mock_filtered_data)
        
        mock_fig = MagicMock()
        mock_fig.data = [MagicMock()]
        mock_figure_widget.return_value = mock_fig
        
        # Mock VBox and HBox to return MagicMocks instead of actual widgets
        mock_vbox.return_value = MagicMock()
        mock_hbox.return_value = MagicMock()
        
        # Execute function
        create_interactive_plot(self.test_data, "Test Plot")
        
        # Verify _create_traces was called with initial parameters
        mock_create_traces.assert_called_with(
            self.test_data,
            selected_term=self.test_data['Term'].unique()[0],
            selected_omic='rna',
            selected_ajcc_stage="All Stages"
        )

if __name__ == '__main__':
    unittest.main()