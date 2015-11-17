#! /usr/bin/env python
#@author Jose Fernandez
""" 
Scripts that takes as input a file generated
in the demultiplexing step (taggd) of the ST pipeline
and computes a count-table of the UMIs presented
in the reads
"""

import sys
import os
import argparse
import pysam
from collections import defaultdict
from stpipeline.common.utils import fileOk
 
def main(sam_filename, outfile, mc_start_position, mc_end_position):

    start_pos = int(mc_start_position)
    end_pos = int(mc_end_position)
    
    if not fileOk(sam_filename) or not sam_filename.endswith(".sam") \
    or start_pos < 0 or start_pos >= end_pos:
        sys.stderr.write("Error, input file not present or invalid format for the parameters\n")
        sys.exit(-1)
    
    if not outfile:
        outfile = "output_counts.txt"
     
    # Count UMIs   
    sam_file = pysam.AlignmentFile(sam_filename, "r")
    umis = defaultdict(int)
    for rec in sam_file:
        seq = str(rec.query_sequence)
        umi = seq[start_pos:end_pos]
        umis[umi] += 1
        
    # Write counts to a file
    with open(outfile, "w") as filehandler:
        for umi,count in umis.iteritems():
            filehandler.write(str(umi) + "\t" + str(count) + "\n")
       
    print "Done"
     
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sam_file", 
                        help="SAM file generated from Taggd in the ST Pipeline")
    parser.add_argument("--outfile", default=None, help="Name of the output file")
    parser.add_argument('--mc-start-position', default=18,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=27,
                        help='Position (base wise) of the last base of the molecular barcodes')
    args = parser.parse_args()
    main(args.sam_file, args.outfile, args.mc_start_position, args.mc_end_position)
