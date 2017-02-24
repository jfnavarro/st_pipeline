""" 
This module contains some functions to find and removes adaptors in fastq reads
"""
import regex

def removeAdaptor(sequence, quality, adaptor, missmatches=2):
    """
    Tries to find the given adaptor sequence in the given fastq read (sequence, quality)
    If adaptor is found removes the adaptor and everything after the adaptor's
    first position and returns the trimmed fastq record.
    :param sequence: the sequence of the read
    :param quality: the quality of the read
    :param adaptor: the adaptor sequence
    :param missmatches: allow number of missmatches when searching for the adaptor
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
    # TODO this is slow, find a faster approach
    candidate = regex.findall(r'(?:%s){s<=%s}' % (adaptor, missmatches), sequence)
    if len(candidate) > 0:
        local_seq = candidate[0]
        # Miss-matches may happen at the start
        # so we account for it
        local_pos = local_seq.find(adaptor[0])
        # We now look for the first base of the matched adaptor
        pos = sequence.find(local_seq[local_pos:])
        # Return the sequence/quality with the adaptor and everything after it trimmed
        return sequence[:pos], quality[:pos]
    else:
        return sequence, quality
            
    
    
    
