#! /usr/bin/env python
""" 
Scripts that parses a SAM or BAM file generated from Taggd and creates
a JSON file with the ST data and a BED file with the reads. It does so
by aggregating the reads for unique gene-barcode tuples. 
"""

import sys
import json
import argparse
import pysam
import numpy as np
import os
from stpipeline.common.utils import getExtension
from collections import defaultdict
from stpipeline.common.clustering import countMolecularBarcodesClustersHierarchical, countMolecularBarcodesClustersNaive, countMolecularBarcodesPrefixtrie

class Transcript:
    """ 
    Simple container for the transcripts
    barcode and gene are str
    x, count and y are int
    reads is a list of tuples (read_name, sequence, quality, chromosome, start, end)
    """
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

def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM file where the barcode and coordinates are present in the tags
    The output will be a list containing unique transcripts (gene,barcode) whose
    reads are aggregated
    """
    unique_events = dict()
    sam_type = getExtension(filename).lower()
    flag = "rb"
    if sam_type == "sam":
        flag = "r"
    sam_file = pysam.AlignmentFile(filename, flag)
    for rec in sam_file:
        clear_name = str(rec.query_name)
        seq = str(rec.query_sequence)
        qual = str(rec.query_qualities)
        mapping_quality = int(rec.mapping_quality)
        start = int(rec.reference_start)
        end = int(rec.reference_end)
        chrom = str(sam_file.getrname(rec.reference_id))
        strand = "+"
        if rec.is_reverse: strand = "-"
        
        for (k, v) in rec.tags:
            if k == "B0":
                barcode = str(v)
            elif k == "B1":
                x = int(v)
            elif k == "B2":
                y = int(v)
            elif k == "XF":
                gene = str(v)
            else:
                continue

        ##TODO should be a way to achieve this with defaultdict
        transcript = Transcript(barcode=barcode, gene=gene, x=x, y=y,
                                count=1, reads=[(seq, chrom, start, end, clear_name, 
                                                 mapping_quality, strand, gene, barcode)])
        try:
            unique_events[transcript] += transcript
        except KeyError:
            unique_events[transcript] = transcript
    
    sam_file.close()
    return unique_events.values()

def main(filename, output_folder,
         output_file_template = None,
         molecular_barcodes = False,
         mc_cluster = "naive",
         allowed_mismatches = 1, mc_start_position = 19, 
         mc_end_position = 27, min_cluster_size = 2):
    
    if filename is None or not os.path.isfile(filename):
        sys.stderr.write("Error, input file not present or invalid: %s\n" % (filename))
        sys.exit(-1)

    sam_type = getExtension(filename).lower()
    if sam_type != "sam" and sam_type != "bam":
        sys.stderr.write("Error, invalid input format: %s\n" % (sam_type))
        sys.exit(-1)
        
    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = "."
    
    if mc_cluster not in ["naive","hierarchical","counttrie"]:
        sys.stderr.write("Error: type of clutering algorithm is incorrect\n")
        sys.exit(-1)
         
    if molecular_barcodes and (mc_start_position < 0  
                               or mc_end_position < 0 or mc_end_position <= mc_start_position):
        sys.stderr.write("Error: Molecular Barcodes option is " \
                         "activated but the start/end positions parameters are incorrect\n")
        sys.exit(-1)
    
    if output_file_template is not None:
        filenameJSON = str(output_file_template) + "_stdata.json"
        filenameReadsBED = str(output_file_template) + "_reads.bed"
    else:
        filenameJSON = "stdata.json"
        filenameReadsBED = "reads.bed"
    
    total_record = 0
    unique_genes = set()
    unique_barcodes = set()
    total_barcodes = 0
    discarded_reads = 0
    max_reads_unique_events = 0
    min_reads_unique_events = 10e6
    barcode_genes = defaultdict(int)
    barcode_reads = defaultdict(int)
    json_transcripts = list()
    
    reads_handler = open(os.path.join(output_folder, filenameReadsBED), "w")
    reads_handler.write("Chromosome\tStart\tEnd\tRead\tScore\tStrand\tGene\tBarcode\n")
    
    unique_events = parseUniqueEvents(filename)
    for transcript in unique_events:                
            # Re-compute the read count accounting for PCR duplicates 
            # if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                ##TODO pass the type to the function so we avoid duplicated code
                if mc_cluster == "counttrie":
                    clusters = countMolecularBarcodesPrefixtrie(transcript.reads,
                                                                allowed_mismatches,
                                                                mc_start_position,
                                                                mc_end_position,
                                                                min_cluster_size)
                elif mc_cluster == "naive":
                    clusters = countMolecularBarcodesClustersNaive(transcript.reads,
                                                                   allowed_mismatches,
                                                                   mc_start_position,
                                                                   mc_end_position,
                                                                   min_cluster_size)
                else:
                    clusters = countMolecularBarcodesClustersHierarchical(transcript.reads,
                                                                          allowed_mismatches,
                                                                          mc_start_position,
                                                                          mc_end_position,
                                                                          min_cluster_size)
                num_clusters = len(clusters)
                transcript.reads = clusters
                discarded_reads += (transcript.count - num_clusters)
                transcript.count = num_clusters
            
            # Add a JSON entry for the transcript
            json_transcripts.append(transcript.toBarcodeDict())
            
            # Add a READS object to a list 
            for read in transcript.reads:
                reads_handler.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(read[1]),
                                                                          str(read[2]),
                                                                          str(read[3]),
                                                                          str(read[4]),
                                                                          str(read[5]),
                                                                          str(read[6]),
                                                                          str(read[7]),
                                                                          str(read[8])))
                
            # Some stats computation
            barcode_genes[transcript.barcode] += 1
            barcode_reads[transcript.barcode] += transcript.count
            unique_genes.add(transcript.gene)
            unique_barcodes.add(transcript.barcode)
            total_record += 1
            total_barcodes += transcript.count
            max_reads_unique_events = max(max_reads_unique_events, transcript.count)
            min_reads_unique_events = min(min_reads_unique_events, transcript.count)
    
    reads_handler.close()
            
    if total_record == 0:
        sys.stderr.write("Error, the number of transcripts present is 0\n")
        sys.exit(-1)
    
    del unique_events
    
    # To compute percentiles
    barcode_genes_array = np.sort(np.array(barcode_genes.values()))
    barcode_reads_array = np.sort(np.array(barcode_reads.values()))
    
    print "Number of Transcripts with Barcode present: " + str(total_barcodes) 
    print "Number of unique events present: " + str(total_record) 
    print "Number of unique Barcodes present: " + str(len(unique_barcodes))
    print "Number of unique Genes present: " + str(len(unique_genes))
    print "Barcode to Genes percentiles: "
    print np.percentile(barcode_genes_array, [0,25,50,75,100])
    print "Barcode to Reads percentiles: "
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
    del json_transcripts    
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in SAM or BAM format containing barcodes and genes from taggd')
    parser.add_argument('--output-folder', type=str,
                        help='Path of the output folder (default is /.)')
    parser.add_argument('--output-file-template', type=str, default=None,
                        help="Describes how the output files will be called")
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, 
                        help="Activates the molecular barcodes duplicates filter")
    parser.add_argument('--mc-cluster', default="naive",
                        help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                        "Modes = {naive(default), counttrie, hierarchical}")
    parser.add_argument('--mc-allowed-mismatches', default=1,
                        help='Number of allowed mismatches when applying the molecular barcodes filter')
    parser.add_argument('--mc-start-position', default=18,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=27,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=2,
                        help='Min number of equal molecular barcodes to count as a cluster')

    args = parser.parse_args()
    main(args.input, args.output_folder, args.output_file_template, 
         args.molecular_barcodes, args.mc_cluster,
         int(args.mc_allowed_mismatches), int(args.mc_start_position),
         int(args.mc_end_position), int(args.min_cluster_size))
                                    
