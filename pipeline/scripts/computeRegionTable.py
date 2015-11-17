#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script that takes regions extracted
from the ST Viewer an generates a table
that contains the barcodes and the regions
they belong to (the regions are taken
from the file names)
"""

import argparse
import sys
import os

def main(barcodes_files, outfile):

    if len(barcodes_files) <= 0:
        sys.stderr.write("Error, input file not present\n")
        sys.exit(-1)
     
    if not outfile:
        outfile = "regions_table.txt"
           
    # reads barcodes and classes
    barcode_to_class = dict()
    for barcode_file in barcodes_files:
        class_name = os.path.splitext(os.path.basename(barcode_file))[0]
        with open(barcode_file, "r") as filehandler:
            for line in filehandler.readlines()[1:]:
                tokens = line.split()
                barcode = tokens[1]
                barcode_to_class[barcode] = class_name
                    
    # writes a BARCODE -> CLASS file
    with open(outfile, "w") as filehandler_write:
        for barcode,class_name in barcode_to_class.iteritems():
            filehandler_write.write(barcode + "\t" + class_name + "\n")
                      
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outfile", default=None, help="Name of the output file")
    parser.add_argument("--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes from the viewer")
    args = parser.parse_args()
    main(args.barcodes_files, args.outfile)

