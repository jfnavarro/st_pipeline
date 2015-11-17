#! /usr/bin/env python
#@Author Jose Fernandez
"""
Script that takes the file with a recomputed
ctts counts table from formatToTable.py and one
or more regions extracted from the viewer and
then filter out from the table the barcodes
that are not in the regions files.
"""

import argparse
import sys
import os
from stpipeline.common.utils import fileOk

def main(counts_file, barcodes_files, outfile):

    if not fileOk(counts_file) or len(barcodes_files) <= 0:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)

    if not outfile:
        outfile = "filtered_" + os.path.basename(counts_file)

    # loads all the barcodes
    barcodes = set()
    for barcode_file in barcodes_files:
        with open(barcode_file, "r") as filehandler:
            for line in filehandler.readlines()[1:]:
                tokens = line.split()
                barcodes.add(tokens[1])

    # writes entries that contain a barcode in the previous list
    with open(counts_file, "r") as filehandler_read:
        with open(outfile, "w") as filehandler_write:
            lines = filehandler_read.readlines()
            first_line = lines.pop(0).strip()
            filehandler_write.write(first_line + "\n")
            for line in lines:
                tokens = line.split()
                bc = tokens[0]
                if bc in barcodes:
                    filehandler_write.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("counts_file",
                        help="Tab delimited file containing the counts table")
    parser.add_argument("--outfile", default=None, help="Name of the output file")
    parser.add_argument("--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes extracted from the ST Viewer")
    args = parser.parse_args()
    main(args.counts_file, args.barcodes_files, args.outfile)