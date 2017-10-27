#! /usr/bin/env python
""" 
Unit-test the package clustering
"""
import unittest
from stpipeline.common.clustering import *
 
class TestClustering(unittest.TestCase):
          
    @classmethod
    def setUpClass(self):
        self.molecular_barcodes1 = ['AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'AAAAAAAABB',
                                    'AAAAAAAABB',
                                    'AAAAAAAABB']
         
        self.molecular_barcodes2 = ['AAAAAAAAAA',
                                    'AAAAAAAABB',
                                    'AAAAAAAACC',
                                    'AAAAAADDDD',
                                    'AAAAAAAAAA',
                                    'AAAAAAAAAA',
                                    'ZZZZAABBBB',
                                    'AAAAAABBBB',
                                    'AAAAAABBBB']
         
        self.molecular_barcodes3 = ['AAAA',
                                    'DAAA',
                                    'ABAA',
                                    'ABBA',
                                    'ABBB',
                                    'BBBB',
                                    'CCCC',
                                    'ACCC']
         
    def test_naive_clustering(self):
        clusters = countUMINaive(self.molecular_barcodes1, 0)
        self.assertTrue(len(clusters) == 2)
        clusters = countUMINaive(self.molecular_barcodes1, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = countUMINaive(self.molecular_barcodes1, 2)
        self.assertTrue(len(clusters) == 1)
        
        clusters = countUMINaive(self.molecular_barcodes2, 0)
        self.assertTrue(len(clusters) == 6)
        clusters = countUMINaive(self.molecular_barcodes2, 1)
        self.assertTrue(len(clusters) == 6)
        clusters = countUMINaive(self.molecular_barcodes2, 2)
        self.assertTrue(len(clusters) == 4)
        
        clusters = countUMINaive(self.molecular_barcodes3, 0)
        self.assertTrue(len(clusters) == 8)
        clusters = countUMINaive(self.molecular_barcodes3, 1)
        self.assertTrue(len(clusters) == 5)
        clusters = countUMINaive(self.molecular_barcodes3, 3)
        self.assertTrue(len(clusters) == 4)
        
    def test_hierarchical_clustering(self):
        clusters = countUMIHierarchical(self.molecular_barcodes1, 0)
        self.assertTrue(len(clusters) == 2)
        clusters = countUMIHierarchical(self.molecular_barcodes1, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = countUMIHierarchical(self.molecular_barcodes1, 2)
        self.assertTrue(len(clusters) == 1)
        
        clusters = countUMIHierarchical(self.molecular_barcodes2, 0)
        self.assertTrue(len(clusters) == 6)
        clusters = countUMIHierarchical(self.molecular_barcodes2, 1)
        self.assertTrue(len(clusters) == 6)
        clusters = countUMIHierarchical(self.molecular_barcodes2, 2)
        self.assertTrue(len(clusters) == 3)
        
        clusters = countUMIHierarchical(self.molecular_barcodes3, 0)
        self.assertTrue(len(clusters) == 8)
        clusters = countUMIHierarchical(self.molecular_barcodes3, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = countUMIHierarchical(self.molecular_barcodes3, 3)
        self.assertTrue(len(clusters) == 1)
        
    def test_dedup_adj(self):
        clusters = dedup_adj(self.molecular_barcodes1, 0)
        self.assertTrue(len(clusters) == 2)
        clusters = dedup_adj(self.molecular_barcodes1, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = dedup_adj(self.molecular_barcodes1, 2)
        self.assertTrue(len(clusters) == 1)
        
        clusters = dedup_adj(self.molecular_barcodes2, 0)
        self.assertTrue(len(clusters) == 6)
        clusters = dedup_adj(self.molecular_barcodes2, 1)
        self.assertTrue(len(clusters) == 6)
        clusters = dedup_adj(self.molecular_barcodes2, 2)
        self.assertTrue(len(clusters) == 4)
        
        clusters = dedup_adj(self.molecular_barcodes3, 0)
        self.assertTrue(len(clusters) == 8)
        clusters = dedup_adj(self.molecular_barcodes3, 1)
        self.assertTrue(len(clusters) == 6)
        clusters = dedup_adj(self.molecular_barcodes3, 3)
        self.assertTrue(len(clusters) == 2)
        
    def test_dedup_dir_adj(self):
        clusters = dedup_dir_adj(self.molecular_barcodes1, 0)
        self.assertTrue(len(clusters) == 2)
        clusters = dedup_dir_adj(self.molecular_barcodes1, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = dedup_dir_adj(self.molecular_barcodes1, 2)
        self.assertTrue(len(clusters) == 1)
        
        clusters = dedup_dir_adj(self.molecular_barcodes2, 0)
        self.assertTrue(len(clusters) == 6)
        clusters = dedup_dir_adj(self.molecular_barcodes2, 1)
        self.assertTrue(len(clusters) == 6)
        clusters = dedup_dir_adj(self.molecular_barcodes2, 2)
        self.assertTrue(len(clusters) == 4)
        
        clusters = dedup_dir_adj(self.molecular_barcodes3, 0)
        self.assertTrue(len(clusters) == 8)
        clusters = dedup_dir_adj(self.molecular_barcodes3, 1)
        self.assertTrue(len(clusters) == 2)
        clusters = dedup_dir_adj(self.molecular_barcodes3, 3)
        self.assertTrue(len(clusters) == 1)
        
    def test_affinity(self):
        # TODO must create test cases
        self.assertTrue(True)
           
if __name__ == '__main__':
    unittest.main()        
