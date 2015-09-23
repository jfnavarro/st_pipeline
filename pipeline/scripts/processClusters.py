#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Small pipelines that does the whole process of :

countClusters.py + formatTaTable.py + filterEntries.py

To obtain a file called ouput_table.bed in R table format
compatible with R, the file will contain an extra column
with the class of the barcodes that were filtered.

"""

import argparse
import sys
import os
import subprocess
import tempfile
from stpipeline.common.utils import fileOk

def countClusters(file_in, min_data, min_size, 
                  min_density, disable_filter, new_paraclu):
    args = ['countClusters.py']
    args += ['--min-data-value', min_data]
    args += ['--min-density-increase', min_density]
    args += ['--max-cluster-size', min_size]
    if disable_filter:
        args += ['--disable-filter']
    if new_paraclu:
        args += ['--new-paraclu']    
    args += [file_in]
    try:
        proc = subprocess.Popen([str(i) for i in args], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        sys.stderr.write("Error, executing countClusters\n")
        sys.exit(-1)   
        
def createDataFrame(file_in, original_file, 
                    use_density=False, out_file=None, new_paraclu=False):
    args = ['formatToTable.py']
    if use_density:
        args += ["--use-density"]
    if out_file is not None:
        args += ["--outfile"]
        args += [out_file]
    if new_paraclu:
        args += ["--new-paraclu"]
    args += [file_in]
    args += [original_file]
    try:
        proc = subprocess.Popen([str(i) for i in args],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        sys.stderr.write("Error, executing formatToTable\n")
        sys.exit(-1)    
    
def filterBarcodes(file_in, barcodes_files, output_file = None):
    args = ['filterEntries.py']
    args += [file_in]
    args += ["--barcodes-files"]
    for bc_file in barcodes_files:
        args += [bc_file]
    if output_file is not None:
        args += ["--outfile", output_file]
    try:
        proc = subprocess.Popen([str(i) for i in args], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        sys.stderr.write("Error, executing filterEntries\n")
        sys.exit(-1)    
           
def main(original_file, barcodes_files, 
         min_data, min_size, min_density, 
         disable_filter, new_paraclu=False):
    
    if not fileOk(original_file) or len(barcodes_files) == 0:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    # Some constants
    file_clusters = "output_clusters.bed"
    file_table = tempfile.mktemp(prefix="st_clusters_data_frame")
    file_table_reduced = "output_table.bed"
    
    # compute clusters
    print "Computing clusters..."
    countClusters(original_file, min_data, min_size, 
                  min_density, disable_filter, new_paraclu)
    # create combined data frame
    print "Converting to table format..."
    createDataFrame(file_clusters, original_file, False, file_table, new_paraclu)
    # Filter Combined table
    print "Filtering entries..."
    filterBarcodes(file_table, barcodes_files, file_table_reduced)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original-file", help="Original file with the ST data in BED format")
    parser.add_argument("--barcodes-files", nargs='+', type=str,
                        help="Tab delimited file containing barcodes from the viewer")
    parser.add_argument("--min-data-value", default=30, 
                        help="Omits grouped entries whose total count is lower than this")
    parser.add_argument("--disable-filter", action="store_true", 
                        default=False, help="Disable second filter(paraclu-cut)")
    parser.add_argument("--max-cluster-size", default=200, 
                        help="Discard clusters whose size in positions is bigger than this")
    parser.add_argument("--min-density-increase", default=2, 
                        help="Discard clusters whose density is lower than this")
    parser.add_argument("--new-paraclu", action="store_true", 
                        default=False, help="Use new paraclu that outputs barcodes")
    args = parser.parse_args()
    main(args.original_file, args.barcodes_files, args.min_data_value, 
         args.disable_filter, args.max_cluster_size, 
         args.min_density_increase, args.new_paraclu)

