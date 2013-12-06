#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" Script to converted demupltiplexed sam file into a tab-delimited file.
"""

import argparse
from main.core import mapping
from main.common.fastq_utils import *
from main.common.utils import *
import sys
import pysam

def createFastqMap(fastq):
    
    fw_file = safeOpenFile(fastq, "rU")
    mapped_reads = dict()
    
    #I create a map of clean read name to fastq record
    for line1 in readfq(fw_file):
        header = line1[0]
        clean_header = header.split()[0]
        mapped_reads[clean_header] = line1
    return mapped_reads
        
        
def main(input, out, fastq):
    
    if input is None or out is None or fastq is None or not os.path.isfile(input) or not os.path.isfile(fastq) or out == "":
        print "Wrong parameters"
        sys.exit(1)
        
    outfile = safeOpenFile(out,"w")    
    input = pysam.Samfile(input, "r")
    
    fastqmap = createFastqMap(fastq)
    
    for read in input:
        # filtering out unmapped reads
        # generating features records Name,chromosome,gene,barcode,x,y,Qul,Read 
        if not read.is_unmapped:
            #sequence = read.seq (original mapped seq)
            #quality = read.qual (original mapped quality)
            #the prunned header (sam format does not allow whitespaces)
            header = read.qname
            #get the original fastq reacord
            original_fastq_record = fastqmap[header]
            #get the name,chr and gene from the original fastq header
            complete_header = original_fastq_record[0]
            original_header_split = complete_header.split()
            name = original_header_split[0]
            chr = original_header_split[1]
            gene = original_header_split[2]
            #get original sequence and quality and remove the fake L
            sequence = original_fastq_record[1].replace("L","")
            quality = original_fastq_record[2].replace("L","")
            #get reference id to get the reference header to extract barcode,x,y
            ref_id = read.tid
            ref_header = input.getrname(ref_id)
            ref_header_split = ref_header.split(":")
            barcode = ref_header_split[0]
            x = ref_header_split[1]
            y = ref_header_split[2]
            
            #print name + "\t" + chr + "\t" + gene + "\t" + barcode + "\t" + x + "\t" + y + "\t" + quality + "\t" + sequence
            outfile.write(name + "\t" + chr + "\t" + gene + "\t" + barcode + "\t" + x + "\t" + y + "\t" + quality + "\t" + sequence + "\n")
        else:
            # not mapped stuff discard
            pass  
            
    input.close()
    outfile.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=str,
                        help='Name of the input file (sam)')
    parser.add_argument('-f', '--input-fastq', type=str,
                        help='Name oft the input file (fastq)')
    parser.add_argument('-o', '--out', type=str,
                        help='Name of merged output file')

    args = parser.parse_args()
    main(args.input, args.out, args.input_fastq)
