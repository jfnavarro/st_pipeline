#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script that takes a file generated from formatToTable.py
and a file generated from stVi with barcodes to filter
out the entries that contains the same barcodes and also
adds an extra column with the class that the barcode belongs to
"""

import argparse
import sys
import os
from collections import defaultdict
from stpipeline.common.utils import fileOk

def main(bed_file, barcodes_files, outfile=None):

    if not fileOk(bed_file) and len(barcodes_files) > 0:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    if outfile is None:
        outfile = "reduced_" + bed_file
           
    # loads all the barcodes
    barcodes = defaultdict(int)
    index_to_class = defaultdict(str)
    index = 0
    for barcode_file in barcodes_files:
        index += 1
        with open(barcode_file, "r") as filehandler:
            class_name = os.path.splitext(barcode_file)[0]
            index_to_class[index] = class_name
            for line in filehandler.readlines():
                tokens = line.split()
                barcodes[tokens[1]] = index     
        
    # writes entries that contain a barcode in the previous list
    with open(bed_file, "r") as filehandler_read:
        with open(outfile, "w") as filehandler_write:
            for line in filehandler_read.readlines():
                tokens = line.split()
                index = barcodes[tokens[0]]
                if index != 0:
                    for ele in tokens:
                        filehandler_write.write(ele + "\t")
                    class_name = index_to_class[index]
                    filehandler_write.write(str(class_name) + "\n")
                    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_file", help="Tab delimited file containing the clusters and the barcodes")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument( "--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes from the viewer")
    args = parser.parse_args()
    main(args.bed_file, args.barcodes_files, args.outfile)
