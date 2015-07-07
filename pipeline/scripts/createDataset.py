#! /usr/bin/env python
""" 
Scripts that parses a SAM or BAM file generated from Taggd and creates
JSON files containing all relevant information. It does so
by aggregating the reads for unique gene-barcode tuples. 
"""

import sys
import os
import json
import argparse
import pysam
import numpy as np
from stpipeline.common.utils import *
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive, countMolecularBarcodesPrefixtrie
from stpipeline.common.fastq_utils import *

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

        #@todo should be a way to achieve this with defaultdict
        transcript = Transcript(barcode=barcode, gene=gene, x=x, y=y,
                                count=1, reads=[(clear_name, seq, qual, chrom, 
                                                 start, end, strand, mapping_quality)])
        if unique_events.has_key(transcript):
            unique_events[transcript] += transcript
        else:
            unique_events[transcript] = transcript

    return unique_events.values()


def main(filename, output_folder, molecular_barcodes = False, use_prefix_tree = False,
         allowed_mismatches = 1, mc_start_position = 19, mc_end_position = 27, min_cluster_size = 2):
    
    if filename is None or not os.path.isfile(filename):
        sys.stderr.write("Error, input file not present or invalid : " + filename + "\n")
        sys.exit(-1)

    sam_type = getExtension(filename).lower()
    if sam_type != "sam" and sam_type != "bam":
        sys.stderr.write("Error, invalid input format : " + sam_type + "\n")
        sys.exit(-1)
        
    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = "."
    
    if molecular_barcodes and (mc_start_position < 0  
                               or mc_end_position < 0 or mc_end_position <= mc_start_position):
        sys.stderr.write("Error: Molecular Barcodes option is " \
                         "activated but the start/end positions parameters are incorrect\n")
        sys.exit(-1)
        
    total_record = 0
    json_barcodes = list()
    json_reads = list()
    unique_genes = set()
    unique_barcodes = set()
    total_barcodes = 0
    discarded_reads = 0
    bed_records = list()
    barcode_genes = dict()
    barcode_reads = dict()
    for transcript in parseUniqueEvents(filename):                
            # Re-compute the read count accounting for PCR duplicates 
            # if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                if use_prefix_tree:
                    clusters = countMolecularBarcodesPrefixtrie(transcript.reads, allowed_mismatches,
                                                                mc_start_position, mc_end_position, 
                                                                min_cluster_size)
                else:
                    clusters = countMolecularBarcodesClustersNaive(transcript.reads, allowed_mismatches,
                                                                   mc_start_position, mc_end_position,
                                                                   min_cluster_size)
                transcript.reads = clusters
                discarded_reads += (transcript.count - len(clusters))
                transcript.count = len(clusters)
            
            # add a JSON entry for the transcript  
            json_barcodes.append(transcript.toBarcodeDict())
            
            # get the reads that mapped to the transcript and generate a JSON file
            # and generates a list of BED records
            for read in transcript.reads:
                qula = read[2]
                seq = read[1]
                name = read[0]
                chrom = read[3]
                start = read[4]
                end = read[5]
                strand = read[6]
                quality_score = read[7]
                json_reads.append({'name': str(name), 
                                   'read': str(seq), 
                                   'quality': str(qula), 
                                   'barcode': transcript.barcode, 
                                   'gene': transcript.gene})
                bed_records.append((chrom, start, end, strand, transcript.gene, 
                                    transcript.barcode, str(name), quality_score))
                
            # Some stats computation
            if barcode_genes.has_key(transcript.barcode):
                barcode_genes[transcript.barcode] += 1
            else:
                barcode_genes[transcript.barcode] = 1
                
            if barcode_reads.has_key(transcript.barcode):
                barcode_reads[transcript.barcode] += transcript.count 
            else:
                barcode_reads[transcript.barcode] = transcript.count
                
            unique_genes.add(str(transcript.gene))
            unique_barcodes.add(str(transcript.barcode))
            total_record += 1
            total_barcodes += int(transcript.count)
    
    if total_record == 0:
        sys.stderr.write("Error, the number of transcripts present is 0\n")
        sys.exit(-1)
    
    # To compute percentiles
    barcode_genes_array = np.array(barcode_genes.values())
    barcode_reads_array = np.array(barcode_reads.values())
    
    print "Number of Transcripts with Barcode present : " + str(total_barcodes) 
    print "Number of unique events present : " + str(total_record) 
    print "Number of unique Barcodes present : " + str(len(unique_barcodes))
    print "Number of unique Genes present : " + str(len(unique_genes))
    print "Barcode to Genes percentiles :"
    print np.percentile(barcode_genes_array, [0,25,50,75,100])
    print "Barcode to Reads percentiles :"
    print np.percentile(barcode_reads_array, [0,25,50,75,100])
    
    if molecular_barcodes:
        print "Number of discarded reads (possible PCR duplicates) : " + str(discarded_reads)
        
    filename = "barcodes.json"
    filenameReads = "reads.json"
    filenameReadsBED = "reads.bed"
    
    # Dump the JSON files to the output files
    with open(os.path.join(output_folder, filename), "w") as filehandler:
        json.dump(json_barcodes, filehandler, indent=2, separators=(',', ': '))  
    with open(os.path.join(output_folder, filenameReads), "w") as filehandlerReads:
        json.dump(json_reads, filehandlerReads, indent=2, separators=(',', ': '))    
    # Dump the reads in BED format
    with open(os.path.join(output_folder, filenameReadsBED), "w") as filehandlerReadsBED:
        filehandlerReadsBED.write("Chromosome\tStart\tEnd\tRead\tScore\tStrand\tGene\tBarcode\n")
        for bed_record in bed_records:
            filehandlerReadsBED.write(str(bed_record[0]) + "\t" \
                                      + str(bed_record[1]) + "\t" + str(bed_record[2]) + "\t" \
                                      + str(bed_record[6]) + "\t" + str(bed_record[7]) + "\t" \
                                      + str(bed_record[3]) + "\t" + str(bed_record[4]) + "\t" \
                                      + str(bed_record[5]) + "\n")
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in SAM or BAM format containing barcodes and genes')
    parser.add_argument('--output-folder', type=str,
                        help='Path of the output folder (default is /.)')
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, 
                        help="Activates the molecular barcodes duplicates filter")
    parser.add_argument('--use-prefix-tree', 
                        action="store_true", default=False, 
                        help="Use a prefix tree instead of naive algorithm for molecular barcodes duplicates filter")
    parser.add_argument('--mc-allowed-mismatches', default=1,
                        help='Number of allowed mismatches when applying the molecular barcodes filter')
    parser.add_argument('--mc-start-position', default=19,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=30,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=2,
                        help='Min number of equal molecular barcodes to count as a cluster')

    args = parser.parse_args()
    main(args.input, args.output_folder, args.molecular_barcodes, args.use_prefix_tree,
         int(args.mc_allowed_mismatches), int(args.mc_start_position),
         int(args.mc_end_position), int(args.min_cluster_size))
                                    
