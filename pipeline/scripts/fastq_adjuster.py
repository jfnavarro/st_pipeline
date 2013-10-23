#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" Script to adjust st fastq reads to have the same length.
"""

import argparse
from main.core import mapping
from main.common.fastq_utils import *
from main.common.utils import *
import sys

def main(input, out, fixed_lenght = 300):
    
    if not os.path.isfile(input) or not not os.path.isfile(out):
        print "Wrong parameters"
        sys.exit(1)
        
    inputfile = safeOpenFile(input,"r")
    outfile = safeOpenFile(out,"w")
    out_writer = writefq(outfile)
    
    for record in readfq(inputfile):
        head = record[0]
        seq = record[1]
        qual = record[2]
        seq_lenght = len(seq)
        while(seq_lenght < fixed_lenght):
            seq += "L"
            qual += "L"
            seq_lenght += 1       
        new_record = (head, seq, qual)
        out_writer.send(new_record)
        
    inputfile.close()
    outfile.close()
    out_writer.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=str,
                        help='Name of the input file (fastq)')
    parser.add_argument('-o', '--out', type=str,
                        help='Name of merged output file')
    parser.add_argument('-l', '--length', type=int, help='fixed length of the reads')

    args = parser.parse_args()
    main(args.input, args.out, args.length)
