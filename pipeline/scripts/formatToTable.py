#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script that parses the tab delimited file
generated from countClusters.py with clusters
and barcodes to a data frame format like this :

    cluster1.....clusterN
BC1 
BC2
...
BCN

It needs the original BED file with ST data to extract the reads count
If no output file is given the output will be : output_table_ctts.txt
"""

import argparse
import sys
from collections import defaultdict
from stpipeline.common.utils import fileOk

def main(input_files, outfile):
    
    if len(input_files) != 2:
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)

    bed_file = input_files[0]
    original_file = input_files[1]
    
    if not fileOk(bed_file) or not fileOk(original_file):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    if outfile is None:
        outfile = "output_table_ctts.txt"
           
    # load all the original barcode - gene coordinates
    map_original_clusters = defaultdict(list)
    with open(original_file, "r") as filehandler:
        for line in filehandler.readlines()[1:]:
            tokens = line.split()
            chromosome = str(tokens[0])
            start_site = str(tokens[1])
            end_site = str(tokens[2])
            strand = str(tokens[5])
            #gene = tokens[6]
            barcode = str(tokens[7])
            map_original_clusters[(chromosome,strand)].append((barcode,int(start_site),int(end_site)))
                            
    # loads all the clusters
    map_clusters = defaultdict(int)
    clusters = set()
    barcodes = set()
    with open(bed_file, "r") as filehandler:
        for line in filehandler.readlines()[1:]:
            tokens = line.split()
            chromosome = str(tokens[0])
            strand = str(tokens[1])
            start = int(tokens[2])
            end = int(tokens[3])
            # doing a full search of intersections over all barcodes
            # If we could rely on that no barcodes were missing doing the clustering we could
            # a faster approach not needing to iterate all the barcodes but only one   
            # this intersection method is prob overcounting
            for barcode_orig, start_orig, end_orig in map_original_clusters[chromosome, strand]:
                if strand == "-":
                    start_orig = end_orig
                if start_orig >= start and start_orig <= end:
                    map_clusters[(barcode_orig,chromosome,strand,start,end)] += 1
                    barcodes.add(barcode_orig) 
            clusters.add((chromosome,strand,start,end))    
    
    # write cluster count for each barcode 
    with open(outfile, "w") as filehandler:
        clusters_string = "\t".join("%s:%s-%s,%s" % cluster for cluster in clusters)
        filehandler.write(clusters_string + "\n")
        for bc in barcodes:
            filehandler.write(bc)
            for chro,strand,star,end in clusters:
                count = map_clusters[(bc,chro,strand,star,end)]
                filehandler.write("\t" + str(count))
            filehandler.write("\n")            
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_files', nargs=2, 
                        help="The tab delimited file containing the clusters and the ST original BED file")
    parser.add_argument("--outfile", default=None, help="Name of the output file")
    args = parser.parse_args()
    main(args.input_files, args.outfile)
