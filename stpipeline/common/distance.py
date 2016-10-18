""" 
This module contains some functions for computing distance between
sequences
"""

def hamming_distance(s1, s2):
    """
    Returns the Hamming distance between equal-length sequences.
    :param s1: the first string/sequence
    :param s1: the second string/sequence
    :type s1: str
    :type s2: str
    :return: the number of number of different elements in both sequences
    :rtype: int
    :raises: ValueError
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
