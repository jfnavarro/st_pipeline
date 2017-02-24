# #! /usr/bin/env python
# """ 
# Unit-test the package clustering
# """
# import unittest
# from stpipeline.common.clustering import countMolecularBarcodesClustersNaive, countMolecularBarcodesClustersHierarchical
# 
# class TestClustering(unittest.TestCase):
#          
#     @classmethod
#     def setUpClass(self):
#         self.molecular_barcodes1 = [('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAABB', []),
#                                     ('AAAAAAAABB', []),
#                                     ('AAAAAAAABB', [])]
#         
#         self.molecular_barcodes2 = [('AAAAAAAAAA', []),
#                                     ('AAAAAAAABB', []),
#                                     ('AAAAAAAACC', []),
#                                     ('AAAAAADDDD', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('AAAAAAAAAA', []),
#                                     ('ZZZZAABBBB', []),
#                                     ('AAAAAABBBB', []),
#                                     ('AAAAAABBBB', [])]
#         
#         self.molecular_barcodes3 = [('AAAA', []),
#                                     ('DAAA', []),
#                                     ('ABAA', []),
#                                     ('ABBA', []),
#                                     ('ABBB', []),
#                                     ('BBBB', []),
#                                     ('CCCC', []),
#                                     ('ACCC', [])]
#         
#     def test_naive_clustering(self):
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes1, 3, 2)
#         # should be one cluster of 10 reads so only a random one is returned
#         self.assertTrue(len(clusters) == 1)
#         
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes1, 3, 10)
#         # should not return any cluster so it will return the 9 reads
#         self.assertTrue(len(clusters) == 9)
# 
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes2, 3, 2)
#         # should make four clusters
#         self.assertTrue(len(clusters) == 4)
#         
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes2, 1, 2)
#         # should make siz clusters
#         self.assertTrue(len(clusters) == 6)
#         
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes2, 5, 2)
#         # should make two clusters
#         self.assertTrue(len(clusters) == 2)
#         
#         clusters = countMolecularBarcodesClustersNaive(self.molecular_barcodes3, 1, 2)
#         # should make 5 clusters
#         self.assertTrue(len(clusters) == 5)
#         
#     def test_hierarchical_clustering(self):
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes1, 3, 2)
#         # should be one cluster 
#         self.assertTrue(len(clusters) == 1)
#         
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes1, 3, 10)
#         # should not return any cluster so it will return 9
#         self.assertTrue(len(clusters) == 9)
# 
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes2, 3, 2)
#         # should make four clusters in complete mode or 3 in single mode
#         self.assertTrue(len(clusters) == 3)
#         
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes2, 1, 2)
#         # should make six clusters
#         self.assertTrue(len(clusters) == 6)
#         
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes2, 5, 2)
#         # should make two clusters
#         self.assertTrue(len(clusters) == 1)
#         
#         clusters = countMolecularBarcodesClustersHierarchical(self.molecular_barcodes3, 1, 2)
#         # should make four clusters in complete mode or 2 in single mode
#         self.assertTrue(len(clusters) == 2)
#           
# if __name__ == '__main__':
#     unittest.main()        
