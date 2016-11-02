""" 
This module contains some functions to find and removes adaptors in fastq reads
"""
import regex

def removeAdaptor(sequence, quality, adaptor, missmatches=2):
    """
    Tries to find the given adaptor in the given fastq read (sequence, quality)
    If adaptor is found removes the adaptor and everything after the adaptor's
    first position.
    :param sequence: the sequence of the read
    :param quality: the quality of the read
    :param adaptor: the adaptor sequence
    :param missmatches: allow for missmatches when finding the adaptor
    :type sequence: str
    :type quality: str
    :type adaptor: str
    :type missmatches: int
    :return: a tuple (sequence,quality) with the adaptor trimmed
    :rtype: tuple
    """
    if len(sequence) < len(adaptor) or len(sequence) != len(quality):
        return sequence, quality
    # Find adaptor and trim from the first position of 
    # the adaptor till the end of the read
    #TODO this is slow, find a faster approach
    candidate = regex.findall(r'(?:%s){s<=%s}' % (adaptor, missmatches), sequence)
    if len(candidate) > 0:
        pos = sequence.find(candidate[0])
        return sequence[:pos], quality[:pos]
    else:
        return sequence, quality
            
    
    
    
