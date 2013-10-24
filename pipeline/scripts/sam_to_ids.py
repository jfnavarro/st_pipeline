#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" Script to converted demupltiplexed sam file into a tab-delimited file.
"""

import argparse
from main.core import mapping
from main.common.fastq_utils import *
from main.common.utils import *
import sys
import pysam

def main(input, out):
    
    if not os.path.isfile(input) or out == "":
        print "Wrong parameters"
        sys.exit(1)
        
    outfile = safeOpenFile(out,"w")    
    input = pysam.Samfile(input, "r")
    
    for read in input:
        # filtering out unmapped reads
        if not read.is_unmapped:
            sequence = read.seq
            quality = read.qual
            header = read.qname
            print header
            print read.mrnm
            print read.rname
            print read.tid
            print read.tags
        else:
            # not mapped stuff discard
            pass  
            
    input.close()
    outfile.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=str,
                        help='Name of the input file (fastq)')
    parser.add_argument('-o', '--out', type=str,
                        help='Name of merged output file')

    args = parser.parse_args()
    main(args.input, args.out)
