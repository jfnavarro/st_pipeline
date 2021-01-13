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
    if missmatches == 0:
        pos = sequence.find(adaptor)
    else:
        candidates = regex.findall(r'(?:%s){s<=%s}' % (adaptor, missmatches),
                                   sequence, overlapped=False)
        if len(candidates) > 0:
            local_seq = candidates[0]
            # Miss-matches may happen at the start
            # so we account for it
            local_pos = 0
            if adaptor[0] != local_seq[0]:
                local_pos = local_seq.find(adaptor[0])
            # We now look for the first base of the matched adaptor
            pos = sequence.find(local_seq[local_pos:])
        else:
            pos = -1
    # Trim only if pos is correct          
    if pos != -1:
        return sequence[:pos], quality[:pos]
    else:
        return sequence, quality
