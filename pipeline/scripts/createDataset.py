#! /usr/bin/env python
#@author Jose Fernandez
""" 
Scripts that parses a SAM or BAM file containing
the aligned reads and the BARCODE,X,Y,UMI as Tags (BZ)
It generates JSON file with the ST data and a BED file with the reads. It does so
by aggregating the reads for unique gene-barcode tuples and removing
duplicates using the UMIs. 
"""
import sys
import json
import argparse
import pysam
import numpy as np
import os
import shelve
try:
    import gdbm
    found_gdbm = True
except ImportError:
    found_gdbm = False
from stpipeline.common.utils import getExtension
from collections import defaultdict
from blist import sorteddict
from stpipeline.common.clustering import countMolecularBarcodesClustersHierarchical, countMolecularBarcodesClustersNaive

class Transcript:
    """ 
    Simple container for the transcripts
    barcode and gene are str
    x, count and y are int
    reads is a list of tuples (read_name, sequence, quality, chromosome, start, end)
    """
    __slots__ = ('barcode','gene','x','y','count','reads')
    def __init__(self, barcode = None, gene = None, x = -1, y = -1, count = 0, reads = []):
        self.barcode = barcode
        self.gene = gene
        self.x = x
        self.y = y
        self.count = count
        self.reads = reads
        
    def __hash__(self):
        return hash((self.barcode, self.gene))
    
    def __eq__(self, other):
        return self.barcode == other.barcode and self.gene == other.gene
        
    def __add__(self, other):
        assert self == other
        self.count += other.count
        self.reads += other.reads
        return self
    
    def __cmp__(self, other):
        return self.barcode == other.barcode and self.gene == other.gene
    
    def __str__(self):
        return "Barcode: %s Gene: %s X: %s Y: %s Hits: %s" % \
            (self.barcode, self.gene, self.x, self.y, self.count)
            
    def toBarcodeDict(self):
        return {'barcode': self.barcode, 'gene': self.gene, 'x': self.x, 'y': self.y, 'hits': self.count}

def parseUniqueEvents(filename, low_memory=False, molecular_barcodes=False):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM file where the barcode and coordinates are present in the tags
    The output will be a list containing unique transcripts (gene,barcode) whose
    reads are aggregated
    :param filename the input file containing the SAM/BAM records
    :param low_memory if True a key-value DB will be used to save memory
    :param mc_start_position the start position of the molecular barcode in the read
    :param mc_end_position the end position of the molecular barcode in the read
    """
    if low_memory:
        #TODO use a global static name variable
        unique_events = shelve.Shelf(gdbm.open("st_temp_hash_create_dataset", "n")) 
    else:
        unique_events = sorteddict()
    sam_type = getExtension(filename).lower()
    flag = "r" if sam_type == "sam" else "rb"
    sam_file = pysam.AlignmentFile(filename, flag)
    for rec in sam_file.fetch(until_eof=True):
        clear_name = str(rec.query_name)
        mapping_quality = int(rec.mapping_quality)
        start = int(rec.reference_start)
        end = int(rec.reference_end)
        chrom = str(sam_file.getrname(rec.reference_id))
        strand = "-" if rec.is_reverse else "+"
        # Get taggd tags
        barcode,x,y,gene,seq = (None,None,None,None,None)
        for (k, v) in rec.tags:
            if k == "B0":
                barcode = str(v)
            elif k == "B1":
                x = int(v) ## The X coordinate
            elif k == "B2":
                y = int(v) ## The Y coordinate
            elif k == "XF":
                gene = str(v)
            elif molecular_barcodes and k == "B3":
                seq = str(v) ## The umi
            else:
                continue
        # Check that all tags are present
        if any(tag is None for tag in [barcode,x,y,gene]) or (molecular_barcodes and not seq):
            sys.stdout.write("Warning: Missing attributes for record %s\n" % clear_name)
            continue
        # Create a new transcript and assign it to dictionary
        transcript = Transcript(barcode=barcode, gene=gene, x=x, y=y, count=1, 
                                reads=[(chrom, start, end, clear_name, mapping_quality, strand, seq)])
        key = gene + barcode
        try:
            unique_events[key] += transcript
        except KeyError:
            unique_events[key] = transcript
            
    sam_file.close()
    unique_transcripts = unique_events.values()
    if low_memory: 
        unique_events.close()
        #TODO use a global static name variable
        if os.path.isfile("st_temp_hash_create_dataset"):
            os.remove("st_temp_hash_create_dataset")
    return unique_transcripts

def main(filename, 
         output_folder,
         output_file_template,
         molecular_barcodes,
         mc_cluster,
         allowed_mismatches,
         min_cluster_size,
         low_memory):
    
    if filename is None or not os.path.isfile(filename):
        sys.stderr.write("Error, input file not present or invalid: %s\n" % (filename))
        sys.exit(-1)

    sam_type = getExtension(filename).lower()
    if sam_type != "sam" and sam_type != "bam":
        sys.stderr.write("Error, invalid input format: %s\n" % (sam_type))
        sys.exit(-1)
        
    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = "."
    
    if mc_cluster not in ["naive","hierarchical"]:
        sys.stderr.write("Error: type of clustering algorithm is incorrect\n")
        sys.exit(-1)
    
    if low_memory and not found_gdbm:
        sys.stdout.write("Warning: low memory option is active but GDBM was not found in your system\n")
        low_memory = False
        
    if output_file_template:
        filenameJSON = str(output_file_template) + "_stdata.json"
        filenameReadsBED = str(output_file_template) + "_reads.bed"
    else:
        filenameJSON = "stdata.json"
        filenameReadsBED = "reads.bed"
    
    total_record = 0
    unique_genes = set()
    unique_barcodes = set()
    total_transcripts = 0
    discarded_reads = 0
    max_reads_unique_events = 0
    min_reads_unique_events = 10e6
    barcode_genes = defaultdict(int)
    barcode_reads = defaultdict(int)
    json_transcripts = list()
    # Obtain the clustering function
    if mc_cluster == "naive":
        clusters_func = countMolecularBarcodesClustersNaive
    else:
        clusters_func = countMolecularBarcodesClustersHierarchical
    # Parse unique events to generate JSON and BED files        
    unique_events = parseUniqueEvents(filename, low_memory,molecular_barcodes)
    
    with open(os.path.join(output_folder, filenameReadsBED), "w") as reads_handler:
        reads_handler.write("Chromosome\tStart\tEnd\tRead\tScore\tStrand\tGene\tBarcode\n")
        for transcript in unique_events:
            # Re-compute the read count accounting for PCR duplicates 
            # if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                # Extract the molecular barcodes
                mcs = [(read[6],read) for read in transcript.reads]
                del transcript.reads
                clusters = clusters_func(mcs,
                                         allowed_mismatches,
                                         min_cluster_size)
                transcript.reads = clusters
                num_clusters = len(clusters)
                discarded_reads += (transcript.count - num_clusters)
                transcript.count = num_clusters
             
            # Add a JSON entry for the transcript
            json_transcripts.append(transcript.toBarcodeDict())
             
            # Add a READS object to a list 
            for read in transcript.reads:
                reads_handler.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(read[0]),
                                                                          str(read[1]),
                                                                          str(read[2]),
                                                                          str(read[3]),
                                                                          str(read[4]),
                                                                          str(read[5]),
                                                                          transcript.gene,
                                                                          transcript.barcode))   
            # Some stats computation
            barcode_genes[transcript.barcode] += 1
            barcode_reads[transcript.barcode] += transcript.count
            unique_genes.add(transcript.gene)
            unique_barcodes.add(transcript.barcode)
            total_record += 1
            total_transcripts += transcript.count
            max_reads_unique_events = max(max_reads_unique_events, transcript.count)
            min_reads_unique_events = min(min_reads_unique_events, transcript.count)
            del transcript
            
    if total_record == 0:
        sys.stderr.write("Error, the number of transcripts present is 0\n")
        sys.exit(-1)
     
    del unique_events
    # To compute percentiles
    barcode_genes_array = np.sort(np.array(barcode_genes.values()))
    barcode_reads_array = np.sort(np.array(barcode_reads.values()))
    # Print some stats
    print "Number of unique transcripts present: " + str(total_transcripts) 
    print "Number of unique events (gene-barcode) present: " + str(total_record) 
    print "Number of unique barcodes present: " + str(len(unique_barcodes))
    print "Number of unique genes present: " + str(len(unique_genes))
    print "Barcode to genes percentiles: "
    print np.percentile(barcode_genes_array, [0,25,50,75,100])
    print "Barcode to reads percentiles: "
    print np.percentile(barcode_reads_array, [0,25,50,75,100])
    print "Max number of genes over all features: " + str(barcode_genes_array[-1])
    print "Min number of genes over all features: " + str(barcode_genes_array[1])
    print "Max number of reads over all features: " + str(barcode_reads_array[-1])
    print "Min number of reads over all features: " + str(barcode_reads_array[1])
    print "Max number of reads over all unique events: " + str(max_reads_unique_events)
    print "Min number of reads over all unique events: " + str(min_reads_unique_events)
    if molecular_barcodes:
        print "Number of discarded reads (possible PCR duplicates): " + str(discarded_reads)
        
    del barcode_genes
    del barcode_reads
    del unique_genes
    del unique_barcodes
     
    # Write JSON transcripts to file
    with open(os.path.join(output_folder, filenameJSON), "w") as json_handler:
        json.dump(json_transcripts, json_handler, indent=2, separators=(",",": "))  
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in SAM or BAM format containing barcodes and genes from taggd')
    parser.add_argument('--output-folder', type=str, default=None,
                        help='Path of the output folder (default is /.)')
    parser.add_argument('--output-file-template', type=str, default=None,
                        help="Describes how the output files will be called")
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, 
                        help="Activates the molecular barcodes duplicates filter")
    parser.add_argument('--mc-cluster', default="naive",
                        help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                        "Modes = {naive(default), hierarchical}")
    parser.add_argument('--mc-allowed-mismatches', default=1,
                        help='Number of allowed mismatches when applying the molecular barcodes filter')
    parser.add_argument('--min-cluster-size', default=2,
                        help='Min number of equal molecular barcodes to count as a cluster')
    parser.add_argument('--low-memory', default=False, action="store_true",
                        help="Writes temporary records into disk in order to save memory but gaining a speed penalty")

    args = parser.parse_args()
    main(args.input, 
         args.output_folder, 
         args.output_file_template, 
         args.molecular_barcodes, 
         args.mc_cluster,
         int(args.mc_allowed_mismatches), 
         int(args.min_cluster_size),
         args.low_memory)
                                    
