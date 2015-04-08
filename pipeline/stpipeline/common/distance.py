#!/usr/bin/env python
""" 
This module contains some functions for computing distance between
sequences
"""

def hamming_distance(s1, s2):
    """
    Returns the Hamming distance between equal-length sequences.
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
