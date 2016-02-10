#!/usr/bin/env python
""" 
This module contains some functions to filter out adaptors in the reads
"""

def removeAdaptor(sequence, quality, adaptor):
    """
    :param sequence the sequence of the read
    :param quality the quality of the read
    :param adaptor is a string containing the adaptor sequence
    Tries to find the given adaptor in the given fastq read (sequence, quality)
    If adaptor is found removes the adaptor and everything after the adaptor's
    first position (this function is meant to be used in reverse reads)
    """
    if len(sequence) < len(adaptor):
        return sequence, quality
    # Find adaptor and trim if found
    pos = sequence.find(adaptor)
    if pos != -1:
        return sequence[:pos], quality[:pos]
    else:
        return sequence, quality
