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
        seq_adaptor_beginning = "TTTTTAAAAAAAAAAAAAAA"
        seq_adaptor_middle = "AAAAAAAAAATTTTTAAAAA"
        seq_adaptor_end = "AAAAAAAAAAAAAAATTTTT"
        adaptor = "TTTTT"
        
        result_seq, result_qual = removeAdaptor(seq_adaptor_beginning, fake_qual, adaptor)
        self.assertTrue(len(result_seq) == 0 and len(result_qual) == 0)
        result_seq, result_qual = removeAdaptor(seq_adaptor_middle, fake_qual, adaptor)
        self.assertTrue(len(result_seq) == 10 and len(result_qual) == 10)
        result_seq, result_qual = removeAdaptor(seq_adaptor_end, fake_qual, adaptor)
        self.assertTrue(len(result_seq) == 15 and len(result_qual) == 15)

        #self.assertRaises(ValueError, removeAdaptor, seq_adaptor_middle, "TTTTT", -1, "3")
        #self.assertRaises(ValueError, removeAdaptor, seq_adaptor_middle, "TTTTT", 21, "3")
        #self.assertRaises(ValueError, removeAdaptor, (fake_name,fake_qual), "TTTTT", 0, "3")

if __name__ == '__main__':
    unittest.main()    