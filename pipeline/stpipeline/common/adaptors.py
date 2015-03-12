#!/usr/bin/env python
""" 
This module contains some functions to handle adaptors in the reads
and filters reads with adaptors
"""

def removeAdaptor(read, adaptor, action="5"):
    """
    Tries to find the given adaptor in the given read
    If adaptor is found three actions can be performed 
      - 5 : removes the adaptor and everything before it from the 5 prime end
      - 3 : removes the adaptor and everything after it from the 3 prime end
      - discard both (returns empty read)
    """
    pos = read.find(adaptor)
    length = len(adaptor)
    if pos != -1 and action == "5":
        return read[pos + length:]
    elif pos != -1 and action == "3":
        return read[:pos + length]
    elif pos != -1 and action == "discard":
        return ""
    else: 
        return read
