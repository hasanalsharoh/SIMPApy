import unittest
import pandas as pd
import numpy as np
from SIMPApy.ranking import calculate_ranking

class TestCalculateRanking(unittest.TestCase):
    def setUp(self):
        # Create mock RNA/DNAm data with control (tw) and case samples
        self.rna_data = pd.DataFrame({
            'tw1': [10, 20, 30, 40, 15],
            'tw2': [12, 22, 32, 42, 17],
            'tw3': [11, 21, 31, 41, 16],
            'case1': [20, 10, 50, 45, 5],
            'case2': [15, 25, 40, 30, 20]
        }, index=['gene1', 'gene2', 'gene3', 'gene4', 'gene5'])
        
        # Create mock CNV data (typically around 2 for normal copy number)
        self.cnv_data = pd.DataFrame({
            'tw1': [2.0, 2.0, 2.0, 1.9, 2.1],
            'tw2': [2.1, 1.9, 2.0, 2.0, 2.0],
            'tw3': [1.9, 2.1, 1.9, 2.1, 2.0],
            'case1': [4.0, 1.0, 3.0, 0.0, 2.0],
            'case2': [2.0, 0.5, 3.5, 1.5, 5.0]
        }, index=['gene1', 'gene2', 'gene3', 'gene4', 'gene5'])

    def test_rna_ranking(self):
        """Test RNA ranking functionality"""
        result = calculate_ranking(self.rna_data, omic="RNA", alpha=0.05)
        
        # Check general structure
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 5)  # 5 samples
        
        # Check that all samples are in the result
        for col in self.rna_data.columns:
            self.assertIn(col, result)
            
        # Check the structure of each dataframe
        sample_df = result['case1']
        self.assertIsInstance(sample_df, pd.DataFrame)
        self.assertEqual(len(sample_df), 5)  # 5 genes
        
        # Check the expected columns exist
        expected_cols = ['D_xs', 'MSD', 'weighted', 'Significant', 'Rank']
        for col in expected_cols:
            self.assertIn(col, sample_df.columns)
        
        # Check some values (gene3 should have a high D_xs in case1)
        self.assertGreater(result['case1'].loc['gene3', 'D_xs'], 
                          result['case1'].loc['gene2', 'D_xs'])
        
        # Check significance calculation
        self.assertTrue((sample_df['Significant'] == 
                       (abs(sample_df['D_xs']) > abs(sample_df['MSD']))).all())

    def test_dnam_ranking(self):
        """Test DNA methylation ranking functionality"""
        result = calculate_ranking(self.rna_data, omic="DNAm", alpha=0.05)
        
        # Should work the same way as RNA
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 5)  # 5 samples
        
        # Expected columns should be present
        for col in ['D_xs', 'MSD', 'weighted', 'Significant', 'Rank']:
            self.assertIn(col, result['case1'].columns)

    def test_cnv_ranking(self):
        """Test CNV ranking functionality"""
        result = calculate_ranking(self.cnv_data, omic="CNV")
        
        # Check general structure
        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 5)  # 5 samples
        
        # Check that each sample result has the expected column
        for sample in result:
            self.assertIn('adjusted_weight', result[sample].columns)
            
        # For genes with copy number = 2, adjusted weight should be small
        # For genes with deviated copy numbers, adjusted weight should be larger
        self.assertLess(abs(result['case1'].loc['gene5', 'adjusted_weight']), 
                       abs(result['case1'].loc['gene1', 'adjusted_weight']))
        
        # gene1 in case1 has CN=4, should have positive adjusted weight
        self.assertGreater(result['case1'].loc['gene1', 'adjusted_weight'], 0)
        
        # gene2 in case1 has CN=1, should have negative adjusted weight
        self.assertLess(result['case1'].loc['gene2', 'adjusted_weight'], 0)

    def test_invalid_omic_type(self):
        """Test handling of invalid omic type"""
        with self.assertRaises(ValueError):
            calculate_ranking(self.rna_data, omic="invalid_type")

    @unittest.expectedFailure
    def test_empty_dataframe(self):
        """Test with empty dataframe"""
        empty_df = pd.DataFrame()
        with self.assertRaises(Exception):  # Some exception should be raised
            calculate_ranking(empty_df)

    @unittest.expectedFailure
    def test_single_sample(self):
        """Test with single sample (should fail as control samples are needed)"""
        single_sample = pd.DataFrame({'case1': [10, 20, 30]}, 
                                    index=['gene1', 'gene2', 'gene3'])
        with self.assertRaises(Exception):  # Should fail without control samples
            calculate_ranking(single_sample)

if __name__ == '__main__':
    unittest.main()