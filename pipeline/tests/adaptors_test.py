#! /usr/bin/env python
""" 
Unit-test the package adaptors
"""

import unittest
from stpipeline.common.adaptors import removeAdaptor

class TestAdaptors(unittest.TestCase):
       
    @classmethod
    def setUpClass(self):
        return
        
    @classmethod
    def tearDownClass(self):
        return
      
    def test_removeAdaptor(self):
        """
        Test that the function removeAdaptor removes correctly
        adaptors when using different parameters.
        It also tests that it fails when it has to
        """
        
        fake_qual = "AAAAAAAAAAAAAAAAAAAA"
        fake_name = "FAKE"
        seq_adaptor_beginning = (fake_name, "TTTTTAAAAAAAAAAAAAAA", fake_qual)
        seq_adaptor_middle = (fake_name, "AAAAAAAAAATTTTTAAAAA", fake_qual)
        seq_adaptor_end = (fake_name, "AAAAAAAAAAAAAAATTTTT", fake_qual)
        
        test_3end = removeAdaptor(seq_adaptor_beginning, "TTTTT", 0, "3")
        self.assertTrue(len(test_3end[1]) == 15 and len(test_3end[2]) == 15)
        test_3end = removeAdaptor(seq_adaptor_middle, "TTTTT", 0, "3")
        self.assertTrue(len(test_3end[1]) == 5 and len(test_3end[2]) == 5)
        test_3end = removeAdaptor(seq_adaptor_end, "TTTTT", 0, "3")
        self.assertTrue(len(test_3end[1]) == 0 and len(test_3end[2]) == 0)
        
        test_5end = removeAdaptor(seq_adaptor_beginning, "TTTTT", 0, "5")
        self.assertTrue(len(test_5end[1]) == 0 and len(test_5end[2]) == 0)
        test_5end = removeAdaptor(seq_adaptor_middle, "TTTTT", 0, "5")
        self.assertTrue(len(test_5end[1]) == 10 and len(test_5end[2]) == 10)
        test_5end = removeAdaptor(seq_adaptor_end, "TTTTT", 0, "5")
        self.assertTrue(len(test_5end[1]) == 15 and len(test_5end[2]) == 15)
        
        #test_discard = removeAdaptor(seq_adaptor_middle, "TTTTT", 0, "discard")
        #self.assertTrue(test_discard is None)
        
        #test_wrong_option = removeAdaptor(seq_adaptor_middle, "TTTTT", 0, "wrong")
        #self.assertTrue(test_wrong_option == seq_adaptor_middle)
        
        #self.assertRaises(ValueError, removeAdaptor, seq_adaptor_middle, "TTTTT", -1, "3")
        #self.assertRaises(ValueError, removeAdaptor, seq_adaptor_middle, "TTTTT", 21, "3")
        #self.assertRaises(ValueError, removeAdaptor, (fake_name,fake_qual), "TTTTT", 0, "3")

if __name__ == '__main__':
    unittest.main()    