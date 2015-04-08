#! /usr/bin/env python
""" 
Scripts that parses a SAM or FASTQ file coming from Taggd to generate
JSON files containing all relevant information. It does so
by aggregating the reads for unique gene-barcode tuples. 
"""

import sys
import os
import json
import argparse
import pysam
from stpipeline.common.utils import *
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive
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
        return "Barcode: %s Gene: %s X: %s Y: %s Hits: %s NumReads: %s" % \
            (self.barcode, self.gene, self.x, self.y, self.count, len(self.reads))
            
    def toBarcodeDict(self):
        return {'barcode': self.barcode, 'gene': self.gene, 'x': self.x, 'y': self.y, 'hits': self.reads}

def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM file where the barcode and coordinates are present in the tags
    The output will be a list containing unique transcripts (gene,barcode) whose
    reads are aggregated
    """
    unique_events = dict()
    sam_file = pysam.AlignmentFile(filename, "r")
    for rec in sam_file:
        clear_name = str(rec.query_name)
        seq = str(rec.query_sequence)
        qual = str(rec.query_qualities)
        #mapping_quality = int(rec.mapping_quality)
        start = int(rec.reference_start)
        end = int(rec.reference_end)
        chrom = str(sam_file.getrname(rec.reference_id))
        strand = "+"
        if rec.mate_is_reverse: strand = "-"
        
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
                                count=1, reads=[(clear_name, seq, qual, chrom, start, end, strand)])
        if unique_events.has_key(transcript):
            unique_events[transcript] += transcript
        else:
            unique_events[transcript] = transcript

    return unique_events.values()


def main(filename, output_name, output_folder, trim_bases = 42, molecular_barcodes = False, 
         allowed_missmatches = 1, mc_start_position = 19, mc_end_position = 27, min_cluster_size = 2):
    
    if  filename is None or output_name is None or not os.path.isfile(filename):
        sys.stderr.write("Error, one of the input file/s not present\n")
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
    
    for transcript in parseUniqueEvents(filename):                
            #re-compute the read count accounting for PCR duplicates 
            #if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                clusters = countMolecularBarcodesClustersNaive(transcript.reads, allowed_missmatches, 
                                               mc_start_position, mc_end_position, min_cluster_size)
                transcript.reads = clusters
                transcript.count = len(clusters)
            
            # add a JSON entry for the transcript  
            json_barcodes.append(transcript.toBarcodeDict())
            
            #get the reads that mapped to the transcript and generate a JSON file
            for read in transcript.reads:
                qula = read[2]
                seq = read[1]
                name = read[0]
                chrom = read[3]
                start = read[4]
                end = read[5]
                strand = read[6]
                json_reads.append({'name': str(name), 
                                   'read': str(seq[trim_bases:]), 
                                   'quality': str(qula[trim_bases:]), 
                                   'barcode': transcript.barcode, 
                                   'gene': transcript.gene})
                bed_records.append((chrom, start, end, strand, transcript.gene, transcript.barcode))
                
            #some stats    
            unique_genes.add(str(transcript.gene))
            unique_barcodes.add(str(transcript.barcode))
            total_record += 1
            total_barcodes += int(transcript.count)
    
    if total_record == 0:
        sys.stderr.write("Error: the number of transcripts present is 0\n")
        sys.exit(-1)
    
    print "Number of Transcripts with Barcode present : " + str(total_barcodes) 
    print "Number of unique events present : " + str(total_record) 
    print "Number of unique Barcodes present : " + str(len(unique_barcodes))
    print "Number of unique Genes present : " + str(len(unique_genes))
    if molecular_barcodes:
        print "Number of discarded reads (possible PCR duplicates) : " + str(discarded_reads)
        
    filename = output_name + "_barcodes.json"
    filenameReads = output_name + "_reads.json"
    filenameReadsBED = output_name + "_reads.bed"
    
    #dump the JSON files to the output files
    with open(os.path.join(output_folder, filename), "w") as filehandler:
        json.dump(json_barcodes, filehandler, indent=2, separators=(',', ': '))  
    with open(os.path.join(output_folder, filenameReads), "w") as filehandlerReads:
        json.dump(json_reads, filehandlerReads, indent=2, separators=(',', ': '))    
    #dump the reads in BED format
    with open(os.path.join(output_folder, filenameReadsBED), "w") as filehandlerReadsBED:
        filehandlerReadsBED.write("Chromosome\tStart\tEnd\tStrand\tGene\tBarcode\n")
        for bed_record in bed_records:
            filehandlerReadsBED.write(str(bed_record[0]) + "\t" \
                                      + str(bed_record[1]) + "\t" + str(bed_record[2]) + "\t" \
                                      + str(bed_record[3]) + "\t" + str(bed_record[4]) + "\t" + str(bed_record[5]) + "\n")
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in FASTQ, SAM or BAM format, augmented with gene.')
    parser.add_argument('--output-folder', type=str,
                        help='Path of the output folder (default is /.)')
    parser.add_argument('--output-name', type=str,
                        help='Name of the output files')
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, help="Activates the molecular barcodes PCR duplicates filter")
    parser.add_argument('--mc-allowed-missmatches', default=1,
                        help='Number of allowed mismatches when applying the molecular barcodes PCR filter')
    parser.add_argument('--mc-start-position', default=19,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=30,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=2,
                        help='Min number of equal molecular barcodes to count as a cluster')
    parser.add_argument('--trim-bases', default=42,
                        help='Number of bases to trim from the output reads')

    args = parser.parse_args()
    main(args.input, args.output_name, args.output_folder,  int(args.trim_bases), args.molecular_barcodes, 
         int(args.mc_allowed_missmatches), int(args.mc_start_position), 
         int(args.mc_end_position), int(args.min_cluster_size))
                                    
