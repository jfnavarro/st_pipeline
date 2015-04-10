#!/usr/bin/env python
""" 
This module contains some functions to handle adaptors in the reads
and filters reads with adaptors
"""

def removeAdaptor(read, adaptor, trimming, action="3"):
    """
    :param read is a tuple (name,sequence,quality)
    :param adaptor is a string containing the adaptor sequence
    :param trimming is the bases to trim in the read (>=0 and <= len(sequence))
    :param action defines that to do if adaptor is found
    Tries to find the given adaptor in the given read
    If adaptor is found three actions can be performed 
      - 3 : removes the adaptor and everything before it from the 3 prime end (default)
      - 5 : removes the adaptor and everything after it from the 5 prime end
      - discard both (returns empty read)
    """
    
    if len(read) != 3:
        raise ValueError("Read must be a tuple (name,sequence,quality)")
    
    if trimming < 0 or trimming > len(read[1]):
        raise ValueError("Incorrect value for trimming")
        
    seq = read[1][trimming:]
    seq_trimmed = read[1][0:trimming]
    name = read[0]
    qual = read[2][trimming:]
    qual_trimmed = read[2][0:trimming]
    pos = seq.find(adaptor)
    length = len(adaptor)
    
    if pos != -1 and action == "3":
        return name, seq_trimmed + seq[pos + length:], qual_trimmed + qual[pos + length:]
    elif pos != -1 and action == "5":
        return name, seq_trimmed + seq[:pos], qual_trimmed + qual[:pos]
    elif pos != -1 and action == "discard":
        return None
    else: 
        return read
