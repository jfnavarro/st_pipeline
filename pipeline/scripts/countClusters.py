#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script that takes as input a BED file generated
from the ST Pipeline and computes peaks clusters
based on the starting site and the strand. 
It uses paraclu to compute the clusters.
Output will be called output_clusters.bed
"""

import argparse
import sys
import os
from collections import defaultdict
from stpipeline.common.utils import fileOk, which
import subprocess
import tempfile

def paracluFilter(file_input, file_output, max_cluster_size, 
                  min_density_increase, single_cluster_density, new_paraclu):
    # sort before filter
    temp_filtered = tempfile.mktemp(prefix="st_countClusters_filtered")
    with open(temp_filtered, "w") as filehandler:
        args = ['sort']
        if new_paraclu:
            args += ["-k9,9"]
        args += ["-k1,1"]
        args += ["-k2,2"]
        args += ["-k3n,3"]
        args += ["-k4nr,4"]
        args += [file_input]
        try:
            proc = subprocess.Popen([str(i) for i in args], stdout=filehandler, stderr=subprocess.PIPE)
            (stdout, errmsg) = proc.communicate()
        except Exception as e:
            sys.stderr.write("Error, sorting\n")
            sys.exit(-1)
                
    # call paraclu-cut to filter + clusters
    with open(file_output, "w") as filehandler:
        args = ['paraclu-cut.sh']
        args += ["-l", max_cluster_size]
        args += ["-d", min_density_increase]
        args += [file_input]
        try:
            proc = subprocess.Popen([str(i) for i in args], stdout=filehandler, stderr=subprocess.PIPE)
            (stdout, errmsg) = proc.communicate()
        except Exception as e:
            sys.stderr.write("Error, executing paraclu-cut\n")
            sys.exit(-1)
            
def appendBarcodes(file_input, file_output, 
                   map_reads_barcodes, map_reads_extended, new_paraclu):
    # Append barcodes to paraclu-cut or put barcode from the end to the start
    with open(file_input, "r") as file_handler_read:
        with open(file_output, "w") as file_handler_write:
            file_handler_write.write("# sequence\tstrand\tstart\tend\tcount\tsum\tmin_density\tmax_density\tbarcode\n")
            for line in file_handler_read.readlines():
                if line.find("#") != -1:
                    continue
                tokens = line.split()
                chromosome = str(tokens[0])
                strand = str(tokens[1])
                star_site = int(tokens[2])
                end_site = str(tokens[3])
                value = str(tokens[4])
                density = str(tokens[5])
                min_density = str(tokens[6])
                max_density = str(tokens[7])
                barcodes = map_reads_barcodes[(chromosome,strand,star_site)]
                for barcode in barcodes:
                    new_value = map_reads_extended[barcode,chromosome,strand,star_site]
                    file_handler_write.write(chromosome + "\t" + strand
                                             + "\t" + str(star_site) + "\t" + end_site + "\t" 
                                             + str(new_value) + "\t" + density + "\t" + min_density 
                                             + "\t" + max_density + "\t" + barcode)
                    file_handler_write.write("\n")
                         
def main(bed_file, min_data_value=30, disable_filter=False, max_cluster_size=200, 
         min_density_increase=2, new_paraclu=False):

    if not fileOk(bed_file):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(-1)
     
    if which("paraclu") is None:
        sys.stderr.write("Error, paraclu was not found in your system\n")
        sys.exit(-1)  
          
    # First we parse the BED file to group by chr,strand and star_site
    # Chromosome      Start   End     Read    Score   Strand  Gene    Barcode
    print "Grouping entries by Chromosome, strand and start site..."
    map_reads = defaultdict(int)
    map_reads_extended = defaultdict(int)
    with open(bed_file, "r") as file_handler:
        for line in file_handler.readlines()[1:]:
            tokens = line.split()
            chromosome = str(tokens[0])
            star_site = int(tokens[1])
            end_site = int(tokens[2])
            strand = str(tokens[5])
            gene = str(tokens[6])
            barcode = str(tokens[7])
            if gene != "ENSMUSG00000092341" and gene.upper() != "MALAT1":
                if strand == "-":
                    star_site = end_site
                map_reads_extended[(barcode,chromosome,strand,star_site)] += 1
                map_reads[(chromosome,strand,star_site)] += 1
     
    print "Printing grouped entries to a temp file..."   
    temp_grouped_reads = tempfile.mktemp(prefix="st_countClusters_grouped_reads")
    with open(temp_grouped_reads, "w") as filehandler:
        # iterate the maps to write out the grouped entries with the barcodes
        if new_paraclu:
            for key,value in map_reads_extended.iteritems():
                barcode = key[0]
                chromosome = key[1]
                strand = key[2]
                star_site = key[3]
                filehandler.write(barcode + "\t" + chromosome + "\t" + strand 
                                  + "\t" + str(star_site) + "\t" + str(value))
                filehandler.write("\n")
        else:
            for key,value in map_reads.iteritems():
                chromosome = key[0]
                strand = key[1]
                star_site = key[2]
                filehandler.write(chromosome + "\t" + strand + "\t" 
                                  + str(star_site) + "\t" + str(value))
                filehandler.write("\n")
    
    # sort the files
    print "Sorting the grouped entries..."
    temp_grouped_reads_sorted = \
    tempfile.mktemp(prefix="st_countClusters_grouped_sorted_reads")
    with open(temp_grouped_reads_sorted, "w") as filehandler:
        args = ['sort']
        if new_paraclu:
            args += ["-k1,1"]
            args += ["-k2,2"]
            args += ["-k3,3"]
            args += ["-k4n,4"]
        else:
            args += ["-k1,1"]
            args += ["-k2,2"]
            args += ["-k3n,3"]
        args += [temp_grouped_reads]
        try:
            proc = subprocess.Popen([str(i) for i in args], 
                                    stdout=filehandler, stderr=subprocess.PIPE)
            (stdout, errmsg) = proc.communicate()
        except Exception as e:
            sys.stderr.write("Error, sorting\n")
            sys.exit(-1)    
    
    # call paraclu to get clusters
    print "Making the peaks calling with paraclu..."   
    clusters_file = "output_clusters.bed"
    with open(clusters_file, "w") as filehandler:
        if new_paraclu:
            args = ['paraclu2']
        else:
            args = ['paraclu']
        args += [min_data_value]
        args += [temp_grouped_reads_sorted]
        try:
            proc = subprocess.Popen([str(i) for i in args], 
                                    stdout=filehandler, stderr=subprocess.PIPE)
            (stdout, errmsg) = proc.communicate()
        except Exception as e:
            sys.stderr.write("Error, executing paraclu\n")
            sys.exit(-1)
            
    if not disable_filter:
        print "Filtering found clusters with paraclu-cut..."
        temp_filtered_clusters = \
        tempfile.mktemp(prefix="st_countClusters_clusters_filtered")
        paracluFilter(clusters_file, temp_filtered_clusters, 
                      max_cluster_size, min_density_increase, False, new_paraclu)          
       
        os.rename(temp_filtered_clusters, clusters_file)
        
    print "DONE!"
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bed_file", help="BED ST-data file")
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
    main(args.bed_file, args.min_data_value, args.disable_filter, args.max_cluster_size, args.min_density_increase, args.new_paraclu)

