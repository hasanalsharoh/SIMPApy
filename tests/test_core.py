import unittest
import pandas as pd
import numpy as np
import os
import tempfile
import shutil
from unittest.mock import patch, MagicMock

from SIMPApy.core import sopa, sopa_population, load_sopa

class TestSopaFunctions(unittest.TestCase):
    
    def setUp(self):
        # Create sample data for testing
        self.genes = ['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7']
        self.ranking = pd.Series([0.8, 0.6, 0.4, 0.2, -0.2, -0.4, -0.6], index=self.genes)
        
        # Sample gene sets
        self.gene_sets = {
            'pathway1': ['gene1', 'gene3', 'gene5'],
            'pathway2': ['gene2', 'gene4', 'gene6'],
            'pathway3': ['gene1', 'gene2', 'gene7']
        }
        
        # Create a temporary directory for output files
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        # Clean up temporary directory
        shutil.rmtree(self.test_dir)
    
    @patch('gseapy.prerank')
    def test_sopa_with_valid_input(self, mock_prerank):
        # Mock the return value of gseapy.prerank
        mock_result = MagicMock()
        mock_result.results = {
            'pathway1': {
                'fdr': 0.01, 'es': 0.5, 'nes': 0.6, 'pval': 0.005,
                'matched_genes': 3, 'gene %': 30, 'lead_genes': 2, 'tag %': 70
            },
            'pathway2': {
                'fdr': 0.05, 'es': 0.3, 'nes': 0.4, 'pval': 0.03,
                'matched_genes': 2, 'gene %': 20, 'lead_genes': 1, 'tag %': 60
            }
        }
        mock_prerank.return_value = mock_result
        
        # Call sopa function
        result = sopa(self.ranking, self.gene_sets)
        
        # Assert that prerank was called with the correct parameters
        mock_prerank.assert_called_once()
        args, kwargs = mock_prerank.call_args
        self.assertEqual(kwargs['min_size'], 3)
        self.assertEqual(kwargs['seed'], 7)
        
        # Assert the result structure
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 2)  # Two pathways
        self.assertIn('Term', result.columns)
        self.assertIn('fdr', result.columns)
        self.assertTrue((result['fdr'] == [0.01, 0.05]).all() or (result['fdr'] == [0.05, 0.01]).all())
    
    @patch('SIMPApy.core.sopa')
    def test_sopa_population(self, mock_sopa):
        # Create a sample ranks DataFrame
        ranks = pd.DataFrame({
            'sample1': [0.8, 0.6, -0.4],
            'sample2': [0.7, -0.5, 0.3]
        }, index=['gene1', 'gene2', 'gene3'])
        
        # Mock return value of sopa
        mock_result = pd.DataFrame({
            'Term': ['pathway1', 'pathway2'],
            'fdr': [0.01, 0.05],
            'es': [0.5, 0.3],
            'nes': [0.6, 0.4],
            'pval': [0.005, 0.03],
            'matched_genes': [3, 2],
            'gene %': [30, 20],
            'lead_genes': [2, 1],
            'tag %': [70, 60]
        })
        mock_sopa.return_value = mock_result
        
        # Call sopa_population
        output_dir = os.path.join(self.test_dir, 'sopa_results')
        sopa_population(ranks, self.gene_sets, output_dir)
        
        # Assert directory was created
        self.assertTrue(os.path.exists(output_dir))
        
        # Assert files were created
        self.assertTrue(os.path.exists(os.path.join(output_dir, "sample1_gsea_results.csv")))
        self.assertTrue(os.path.exists(os.path.join(output_dir, "sample2_gsea_results.csv")))
        
        # Assert sopa was called twice (once for each sample)
        self.assertEqual(mock_sopa.call_count, 2)
    
    def test_load_sopa(self):
        # Create mock result files
        os.makedirs(os.path.join(self.test_dir, 'results'), exist_ok=True)
        
        # Create mock result CSV files
        sample1_df = pd.DataFrame({
            'Term': ['pathway1', 'pathway2'],
            'fdr': [0.01, 0.05],
            'es': [0.5, 0.3],
            'nes': [0.6, 0.4],
            'pval': [0.005, 0.03],
            'matched_genes': [3, 2],
            'gene %': [30, 20],
            'lead_genes': ['gene1', 'gene2'],
            'tag %': [70, 60]
        })
        
        sample2_df = pd.DataFrame({
            'Term': ['pathway1', 'pathway3'],
            'fdr': [0.02, 0.04],
            'es': [0.4, 0.3],
            'nes': [0.5, 0.4],
            'pval': [0.01, 0.02],
            'matched_genes': [2, 3],
            'gene %': [20, 30],
            'lead_genes': ['gene2', 'gene3'],
            'tag %': [60, 70]
        })
        
        # Save mock result files
        sample1_df.to_csv(os.path.join(self.test_dir, 'results', 'tm1_gsea_results.csv'), index=False)
        sample2_df.to_csv(os.path.join(self.test_dir, 'results', 'tw1_gsea_results.csv'), index=False)
        
        # Call load_sopa
        results = load_sopa(os.path.join(self.test_dir, 'results'))
        
        # Assert results structure
        self.assertIsInstance(results, pd.DataFrame)
        self.assertEqual(len(results), 4)  # 2 pathways × 2 samples
        self.assertIn('sample_name', results.columns)
        self.assertIn('Term', results.columns)
        self.assertIn('fdr', results.columns)
        
        # Assert sample names are correctly extracted
        self.assertTrue('tm1' in results['sample_name'].values)
        self.assertTrue('tw1' in results['sample_name'].values)
    
    @unittest.expectedFailure
    def test_load_sopa_with_invalid_file(self):
        # Create directory with an invalid CSV file
        os.makedirs(os.path.join(self.test_dir, 'invalid'), exist_ok=True)
        
        # Create a valid file
        pd.DataFrame({'Term': ['pathway1']}).to_csv(
            os.path.join(self.test_dir, 'invalid', 'tm1_gsea_results.csv'), index=False
        )
        
        # Create an invalid file
        with open(os.path.join(self.test_dir, 'invalid', 'tm2_gsea_results.csv'), 'w') as f:
            f.write('This is not a valid CSV file"""""')
        
        # Test that it handles the error gracefully
        results = load_sopa(os.path.join(self.test_dir, 'invalid'))
        
        # Only the valid file should be processed
        self.assertEqual(len(results), 1)
        self.assertEqual(results['sample_name'].iloc[0], 'tm1')

if __name__ == '__main__':
    unittest.main()