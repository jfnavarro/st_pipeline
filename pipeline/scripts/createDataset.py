#! /usr/bin/env python
""" 
    Scripts that parses a tab delimited file with the format (Name,chromosome,gene,barcode,x,y,Qul,Read)
    into two json files one containing the features and another one containing the reads
"""

import sys
import os
import json
import argparse
from stpipeline.common.utils import *
from stpipeline.common.clustering import countMolecularBarcodesClustersNaive
from stpipeline.common.fastq_utils import *

import pysam

class Transcript:
    """ 
    Simple container for the transcripts
    """
    def __init__(self, barcode = None, gene = None, x = -1, y = -1, reads = 0, 
                 sequence = [], quality = [], readName = []):
        self.barcode = barcode
        self.gene = gene
        self.x = x
        self.y = y
        self.reads = reads
        self.sequences = sequence
        self.qualities = quality
        self.readNames = readName
        
    def __hash__(self):
        return hash((self.barcode, self.gene))
    
    def __eq__(self, other):
        return self.barcode == other.barcode and self.gene == other.gene
        
    def __add__(self, other):
        assert self == other
        self.reads += other.reads
        self.sequences += other.sequences
        self.qualities += other.qualities
        self.readNames += other.readNames
        return self
    
    def __cmp__(self, other):
        return self.barcode == other.barcode and self.gene == other.gene
    
    def __str__(self):
        return "Barcode: %s Gene: %s X: %s Y: %s Hits: %s NumReads: %s" % \
            (self.barcode, self.gene, self.x, self.y, self.reads, len(self.sequences))
            
    def toBarcodeDict(self):
        return {'barcode': self.barcode, 'gene': self.gene, 'x': self.x, 'y': self.y, 'hits': self.reads}

def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    Expected fastq format as follows for the description line:
    @<annotation> Chr:<chromosome> Gene:<gene> B0:Z:<barcode> B1:Z:<x> B2:Z:<y>
    The output will be a list containing unique transcripts (gene,barcode) whose
    reads are aggregated
    """

    unique_events = dict()

    suffix = getExtension(filename).lower()

    if suffix == "fq" or suffix == "fastq":
        # FASTQ
        with open(filename, 'r') as filehandler:
            for (descr, seq, qual) in readfq(filehandler):
                cols = descr.split()
                assert len(cols) == 6
                clear_name = str(cols[0].replace("@",""))
                chromosome = str(cols[1].replace("Chr:",""))
                gene = str(cols[2].replace("Gene:",""))
                barcode = str(cols[3].replace("B0:Z:", ""))
                x = int(cols[4].replace("B1:Z:", ""))
                y = int(cols[5].replace("B2:Z:", ""))
                seq = str(seq)
                qual = str(qual)

                #@todo should be a way to achieve this with defaultdict
                transcript = Transcript(barcode=barcode, gene=gene, x=x, y=y,
                                        reads=1, sequence=[seq], quality=[qual], readName=[clear_name])
                if unique_events.has_key(transcript):
                    unique_events[transcript] += transcript
                else:
                    unique_events[transcript] = transcript

    elif suffix == "sam" or suffix == "bam":
        # SAM/BAM FILES.
        if suffix == "sam":
            flg = "r"
        else:
            flg = "rb"
        with pysam.AlignmentFile(filename, flg) as filehandler:
            for rec in filehandler:
                clear_name = str(rec.query_name)
                seq = str(rec.query_sequence)
                qual = str(rec.query_qualities)
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
                                        reads=1, sequence=[seq], quality=[qual], readName=[clear_name])
                if unique_events.has_key(transcript):
                    unique_events[transcript] += transcript
                else:
                    unique_events[transcript] = transcript
    else:
        raise ValueError("Invalid file format: expected FASTQ, SAM or BAM.")

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
    
    for transcript in parseUniqueEvents(filename):
            #re-compute the read count accounting for PCR duplicates 
            #if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                clusters = countMolecularBarcodesClustersNaive(transcript.sequences, allowed_missmatches, 
                                               mc_start_position, mc_end_position, min_cluster_size)
                reads_covered_by_clusters = sum([len(x) for x in clusters])
                #adjust the transcript's reads count by the difference of 
                #the original transcript's reads count minus the total number of reads clustered
                #and plus the number of clusters (each cluster counts as one read)
                adjusted_reads = len(clusters) + (transcript.reads - reads_covered_by_clusters)
                discarded_reads += (transcript.reads - adjusted_reads)
                transcript.reads = adjusted_reads
            
            # add a JSON entry for the transcript  
            json_barcodes.append(transcript.toBarcodeDict())
            
            #get the reads that mapped to the transcript and generate a JSON file
            for qula, read, name in zip(transcript.qualities, transcript.sequences, transcript.readNames):
                json_reads.append({'name': str(name), 'read': str(read[trim_bases:]), 
                                   'quality': str(qula[trim_bases:]), 'barcode': transcript.barcode, 'gene': transcript.gene})
                
            #some stats    
            unique_genes.add(str(transcript.gene))
            unique_barcodes.add(str(transcript.barcode))
            total_record += 1
            total_barcodes += int(transcript.reads)
    
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

    #dump the JSON files to the output files
    with open(os.path.join(output_folder, filename), "w") as filehandler:
        json.dump(json_barcodes, filehandler, indent=2, separators=(',', ': '))  
    with open(os.path.join(output_folder, filenameReads), "w") as filehandlerReads:
        json.dump(json_reads, filehandlerReads, indent=2, separators=(',', ': '))    
        
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
                                    
