#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Scripts that parses the different files generated 
with countClusters to convert them to R data frame
and filter the entries with the barcodes given

"""

import argparse
import sys
import subprocess
from stpipeline.common.utils import fileOk

def createDataFrame(file, original_file, use_density=False, out_file=None, new_paraclu=False):
    args = ['formatToTable.py']
    if use_density:
        args += ["--use-density"]
    if out_file is not None:
        args += ["--outfile"]
        args += [out_file]
    if new_paraclu:
        args += ["--new-paraclu"]
    args += [file]
    args += [original_file]
    try:
        proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        sys.stderr.write("Error, executing formatToTable\n")
        sys.exit(-1)    
    
def filterBarcodes(file, barcodes_files):
    args = ['filterEntries.py']
    args += [file]
    args += ["--barcodes-files"]
    for bc_file in barcodes_files:
        args += [bc_file]
    try:
        proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        sys.stderr.write("Error, executing filterEntries\n")
        sys.exit(-1)    
           
def main(clusters_file, original_file, barcodes_files, new_paraclu=False):
    
    if not fileOk(clusters_file) or not fileOk(original_file) or len(barcodes_files) == 0:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    # Some constants
    file_table = "output_table.bed"
    file_density_table = "output_table_density.bed"
    
    # create combined data frame
    createDataFrame(clusters_file, original_file, False, file_table, new_paraclu)
    # create combined density data frame
    createDataFrame(clusters_file, original_file, True, file_density_table, new_paraclu)
    # Filter Combined table
    filterBarcodes(file_table, barcodes_files)
    # Filter Combined density table
    filterBarcodes(file_density_table, barcodes_files)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('clusters_file', 
                        help="The files with clusters and barcodes")
    parser.add_argument("--original-file", help="Original file with the ST data in BED format")
    parser.add_argument( "--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes from the viewer")
    parser.add_argument("--new-paraclu", action="store_true", 
                        default=False, help="Use new paraclu that outputs barcodes")
    args = parser.parse_args()
    main(args.clusters_file, args.original_file, args.barcodes_files, args.new_paraclu)

