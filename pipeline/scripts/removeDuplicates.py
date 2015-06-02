#! /usr/bin/env python
""" 
Scripts that parses a fastq file as input and remove
duplicates using UMI
It needs parameters for the UMI length and positions
@author Jose Fernandez
"""

import sys
import os
import argparse
from stpipeline.common.utils import *
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive
from stpipeline.common.fastq_utils import readfq, writefq

def main(filename, allowed_mismatches = 1, mc_start_position = 29,
         mc_end_position = 37, min_cluster_size = 2):
    
    if filename is None or not os.path.isfile(filename) \
    or not filename.endswith("fastq"):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
    
    if mc_start_position < 0  or mc_end_position < 0 or mc_end_position <= mc_start_position:
        sys.stderr.write("Error: UMI parameters are invalid\n")
        sys.exit(-1)
        
    out_file = replaceExtension(getCleanFileName(filename),'_clean.fastq')    
    out_handle = safeOpenFile(out_file, 'w')
    out_writer = writefq(out_handle)
    in_file = safeOpenFile(filename, "rU")
    
    reads = list()
    for line in readfq(in_file):
        reads.append(line)
        
    clusters = countMolecularBarcodesClustersNaive(reads, 
                                                   allowed_mismatches,
                                                   mc_start_position, 
                                                   mc_end_position, 
                                                   min_cluster_size)
    for read in clusters:
        out_writer.send(read)

    out_handle.close()
    in_file.close()

    print "Number of reads present " + str(len(reads))
    print "Number of reads present after removing duplicates " + str(len(clusters))
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in FASTQ')
    parser.add_argument('--mc-allowed-mismatches', default=1,
                        help='Number of allowed mismatches when applying the molecular barcodes PCR filter')
    parser.add_argument('--mc-start-position', default=29,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=37,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=2,
                        help='Min number of equal molecular barcodes to count as a cluster')

    args = parser.parse_args()
    main(args.input,
         int(args.mc_allowed_mismatches),
         int(args.mc_start_position), 
         int(args.mc_end_position), 
         int(args.min_cluster_size))
                                    
