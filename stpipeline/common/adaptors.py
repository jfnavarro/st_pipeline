""" 
This module contains some functions to find and removes adaptors in fastq reads
"""

def removeAdaptor(sequence, quality, adaptor):
    """
    Tries to find the given adaptor in the given fastq read (sequence, quality)
    If adaptor is found removes the adaptor and everything after the adaptor's
    first position.
    :param sequence: the sequence of the read
    :param quality: the quality of the read
    :param adaptor: the adaptor sequence
    :type sequence: str
    :type quality: str
    :type adaptor: str
    :return: a tuple (sequence,quality) with the adaptor trimmed
    :rtype: tuple
    """
    if len(sequence) < len(adaptor) or len(sequence) != len(quality):
        return sequence, quality
    # Find adaptor and trim if found
    pos = sequence.find(adaptor)
    if pos != -1:
        return sequence[:pos], quality[:pos]
    else:
        return sequence, quality
