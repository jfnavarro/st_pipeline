#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" Script to generate a fastq file from ids files.
"""
import argparse
from main.common.utils import *
import sys

def main(input, out):
    
    if not os.path.isfile(input) or out == "":
        print "Wrong parameters"
        sys.exit(1)
        
    inputfile = safeOpenFile(input,"r")
    outfile = safeOpenFile(out,"w")
    
    for line in inputfile.readlines():
        split = line.split()
        head = str(split[0]) + ":" + str(split[1]) + ":" + str(split[2]) 
        seq = str(split[0])
        outfile.write(">" + head + "\n")
        outfile.write(seq + "\n")
        
    inputfile.close()
    outfile.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=str,
                        help='Name of the input file (ids)')
    parser.add_argument('-o', '--out', type=str,
                        help='Name of merged output file')

    args = parser.parse_args()
    main(args.input, args.out)
