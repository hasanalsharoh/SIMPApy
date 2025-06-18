import unittest
import pandas as pd
from unittest.mock import patch

from SIMPApy.analyze import group_diffs, plot_volcano, calculate_correlation, plot_correlation_scatterplot

class TestAnalyzeFunctions(unittest.TestCase):
    
    def setUp(self):
        """Set up a sample DataFrame for testing."""
        self.sample_data = pd.DataFrame({
            'Term': ['pathway1', 'pathway1', 'pathway1', 'pathway1', 'pathway2', 'pathway2', 'pathway2', 'pathway2'],
            'sample_name': ['tm1', 'tm2', 'tw1', 'tw2', 'tm1', 'tm2', 'tw1', 'tw2'],
            'mpes': [1.2, 1.5, -1.1, -1.3, 0.5, 0.6, -0.4, -0.7],
            'value1': [1, 2, 3, 4, 5, 6, 7, 8],
            'value2': [8, 7, 6, 5, 4, 3, 2, 1]
        })

    def test_group_diffs(self):
        """Test the group_diffs function for correctness."""
        diff_results = group_diffs(
            data=self.sample_data,
            pathway_col='Term',
            value_col='mpes',
            group_col='sample_name',
            group1_prefix='tm',
            group2_prefix='tw'
        )
        self.assertIsInstance(diff_results, pd.DataFrame)
        self.assertEqual(len(diff_results), 2)
        self.assertIn('pathway', diff_results.columns)
        self.assertIn('p_adj', diff_results.columns)
        self.assertIn('neg_log10_p_adj', diff_results.columns)
        # For this specific data, tm groups have higher mpes, so mean_diff should be positive
        self.assertTrue((diff_results['mean_diff'] > 0).all())

    @patch('matplotlib.pyplot.show')
    def test_plot_volcano(self, mock_show):
        """Test that the volcano plot function runs without error."""
        diff_results = group_diffs(self.sample_data, 'Term', 'mpes', 'sample_name')
        plot_volcano(diff_results)
        mock_show.assert_called_once()

    def test_calculate_correlation(self):
        """Test the correlation calculation."""
        # Test overall correlation
        corr_df = calculate_correlation(self.sample_data, 'value1', 'value2')
        self.assertEqual(len(corr_df), 1)
        self.assertAlmostEqual(corr_df['correlation'].iloc[0], -1.0)
        
        # Test grouped correlation
        self.sample_data['group'] = ['A', 'A', 'B', 'B', 'A', 'A', 'B', 'B']
        corr_df_grouped = calculate_correlation(self.sample_data, 'value1', 'value2', group_col='group')
        self.assertEqual(len(corr_df_grouped), 2)
        # Correlation within each group should also be -1.0 for this data
        self.assertAlmostEqual(corr_df_grouped.loc[corr_df_grouped['group'] == 'A', 'correlation'].iloc[0], -1.0)
        self.assertAlmostEqual(corr_df_grouped.loc[corr_df_grouped['group'] == 'B', 'correlation'].iloc[0], -1.0)

    @patch('matplotlib.pyplot.show')
    def test_plot_correlation_scatterplot(self, mock_show):
        """Test that the scatter plot function runs without error."""
        plot_correlation_scatterplot(
            data=self.sample_data,
            x_col='value1',
            y_col='value2',
            hue_col='Term'
        )
        # lmplot calls show() on the FacetGrid object, not directly on pyplot
        # So we check if any plot was attempted to be shown
        self.assertGreater(mock_show.call_count, 0)
        
if __name__ == '__main__':
    unittest.main()