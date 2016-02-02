#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script that takes the reads.bed file generated
from the ST pipeline and a bunch of selections
made in the ST Viewer and filters the reads.bed
with the the barcodes that are present in the selections.
"""

import argparse
import sys
import os
from collections import defaultdict
from stpipeline.common.utils import fileOk

def main(bed_file, barcodes_files, outfile):

    if not fileOk(bed_file) or len(barcodes_files) <= 0:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    if not outfile:
        outfile = "filtered_" + os.path.basename(bed_file)
           
    # loads all the barcodes
    barcodes = set()
    for barcode_file in barcodes_files:
        with open(barcode_file, "r") as filehandler:
            for line in filehandler.readlines()[1:]:
                tokens = line.split()
                barcodes.add(tokens[1])
        
    # writes entries that contain a barcode in the previous list
    with open(bed_file, "r") as filehandler_read:
        with open(outfile, "w") as filehandler_write:
            lines = filehandler_read.readlines()
            first_line = lines.pop(0).strip()
            filehandler_write.write(first_line + "\n")
            for line in lines:
                tokens = line.split()
                bc = tokens[7]
                if bc in barcodes:
                    filehandler_write.write(line)
                    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_file", 
                        help="Tab delimited file containing the clusters and the barcodes")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument("--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes from the viewer")
    args = parser.parse_args()
    main(args.bed_file, args.barcodes_files, args.outfile)

