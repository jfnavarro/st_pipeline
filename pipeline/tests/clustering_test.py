#! /usr/bin/env python
""" 
Unit-test the package clustering
"""
from operator import itemgetter
import unittest
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive

class TestClustering(unittest.TestCase):
         
    def test_naive_clustering(self):
        """
        Tests the function numberOfClusters that takes
        as input a list of short reads(molecular barcodes) the minimum distance
        allowed and the min cluster size allowed
        """
        molecular_barcodes = [('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'),
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'),
                              ('name','AAAAAAAABB','QQQQQQQQQQ'),
                              ('name','AAAAAAAABB','QQQQQQQQQQ'), 
                              ('name','AAAAAAAABB','QQQQQQQQQQ')]
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, 10, 2)
        # should be one cluster of 10 reads so only a random one is returned
        self.assertTrue(len(clusters) == 1)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, 10, 10)
        # should not return any cluster so it will return the 10 reads
        self.assertTrue(len(clusters) == 9)

        molecular_barcodes = [('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAABB','QQQQQQQQQQ'), 
                              ('name','AAAAAAAACC','QQQQQQQQQQ'),
                              ('name','AAAAAADDDD','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'), 
                              ('name','AAAAAAAAAA','QQQQQQQQQQ'),
                              ('name','ZZZZAABBBB','QQQQQQQQQQ'),
                              ('name','AAAAAABBBB','QQQQQQQQQQ'), 
                              ('name','AAAAAABBBB','QQQQQQQQQQ')]
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 3, 0, 10, 2)
        # should make two cluster, one of 5 reads and one of 2 so the output should have 4 reads
        self.assertTrue(len(clusters) == 4)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 1, 0, 10, 2)
        # should make two clusters one of 3 reads and another one of 2 so the output should have 6 reads
        self.assertTrue(len(clusters) == 6)
        
        clusters = countMolecularBarcodesClustersNaive(molecular_barcodes, 5, 0, 10, 2)
        # should make one cluster of 8 reads so the output should have 2 reads
        self.assertTrue(len(clusters) == 2)

    def test_prefix_tree_clustering(self):
        #TODO
        self.assertTrue(True)
        
    def test_hierarchical_clustering(self):
        #TODO
        self.assertTrue(True)
        
        
if __name__ == '__main__':
    unittest.main()        
