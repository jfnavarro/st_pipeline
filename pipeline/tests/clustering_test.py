#! /usr/bin/env python
""" 
Unit-test the package clustering
"""

import unittest
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive

class TestClustering(unittest.TestCase):
         
    def test_naive_clustering(self):
        """
        Tests the function numberOfClusters that takes
        as input a list of short reads(molecular barcodes) the minimum distance
        allowed and the min cluster size allowed
        """
        molecular_barcodes = ['AAAAAAAAAA', 'AAAAAAAAAA', 'AAAAAAAAAA',
                              'AAAAAAAAAA', 'AAAAAAAAAA', 'AAAAAAAAAA',
                              'AAAAAAAABB', 'AAAAAAAABB', 'AAAAAAAABB']
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, len("AAAAAAAAAA"), 2)
        self.assertTrue(len(clusters) == 1 and len(clusters[0]) == 9)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, len("AAAAAAAAAA"), 10)
        self.assertTrue(len(clusters) == 0)
      
        molecular_barcodes = ['AAAAAAAAAA', 'AAAAAAAABB', 'AAAAAAAACC', 'AAAAAADDDD', 'AAAAAAAAAA', 
                              'AAAAAAAAAA', 'ZZZZAABBBB', 'AAAAAABBBB', 'AAAAAABBBB']
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, len("AAAAAAAAAA"), 2)
        self.assertTrue(len(clusters) == 2 and len(clusters[0]) == 5 and len(clusters[1]) == 2)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 1, 0, len("AAAAAAAAAA"), 2)
        self.assertTrue(len(clusters) == 2 and len(clusters[0]) == 3 and len(clusters[1]) == 2)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 5, 0, len("AAAAAAAAAA"), 2)
        self.assertTrue(len(clusters) == 1 and len(clusters[0]) == 8)

if __name__ == '__main__':
    unittest.main()        
