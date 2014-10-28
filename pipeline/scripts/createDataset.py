#! /usr/bin/env python
""" 
    Scripts that parses a tab delimited file with the format (Name,chromosome,gene,barcode,x,y,Qul,Read)
    into two json files one containing the features and another one containing the reads
"""
import sys
import os
import json
import argparse

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
    Expected tab-delimited format as follows : 
    [read_name | chromosome | gene | barcode | x | y | read_quality | read_sequence]
    The output will be a list containing unique transcripts (gene,barcode) whose
    reads are aggregated
    @todo : add a proper CSV parser with column names
    """

    unique_events = dict()
    
    with open(filename, 'r') as filehandler:
        for line in filehandler:
            cols = line.split()
            assert len(cols) == 8
            gene = str(cols[2].replace("Gene:",""))
            clear_name = str(cols[0].replace("@",""))
            #chromosome = str(cols[1].replace("Chr:",""))
            barcode = str(cols[3])
            x = int(cols[4])
            y = int(cols[5])
            qual = str(cols[6])
            seq = str(cols[7])
            
            #@todo should be a way to achieve this with defaultdict
            transcript = Transcript(barcode=barcode, gene=gene, x=x, y=y, 
                                    reads=1, sequence=[seq], quality=[qual], readName=[clear_name])
            if unique_events.has_key(transcript):
                unique_events[transcript] += transcript
            else:
                unique_events[transcript] = transcript
            
    return unique_events.values()

def hamming_distance(s1, s2):
    """
    Returns the Hamming distance between equal-length sequences.
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def extractMolecularBarcodes(reads, mc_start_position, mc_end_position):
    """ 
    Extracts a list of molecular barcodes from the list of reads given their
    start and end positions
    """
    molecular_barcodes = list()
    for read in reads:
        if mc_end_position > len(read):
            sys.stderr.write("Error, molecular barcode could not be found in the read " + read + "\n")
            sys.exit(-1)
        molecular_barcodes.append(read[mc_start_position:mc_end_position])
    return molecular_barcodes

def numberOfClusters(molecular_barcodes, allowed_missmatches, min_cluster_size):
    """
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_missmatches)
    It will return a list with the all the clusters and their size
    """
    molecular_barcodes.sort()
    nclusters = list()
    local_cluster_size = 0
    
    for i in xrange(0, len(molecular_barcodes) - 1):
        distance = hamming_distance(molecular_barcodes[i], molecular_barcodes[i + 1]) 
        
        if distance < allowed_missmatches:
            local_cluster_size += 1
        else:
            if local_cluster_size >= min_cluster_size:
                nclusters.append(local_cluster_size + 1)
            local_cluster_size = 0
            
    if local_cluster_size > 0:
        nclusters.append(local_cluster_size + 1)
        
    return nclusters

def removePCRduplicates(reads, allowed_missmatches, mc_start_position, mc_end_position, min_cluster_size):
    """ 
    Returns a list wit the number of clusters and their sizes obtained from the molecular
    barcodes present in the reads
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    return numberOfClusters(molecular_barcodes, allowed_missmatches, min_cluster_size)

def main(filename, output_name, output_folder, molecular_barcodes = False, 
         allowed_missmatches = 3, mc_start_position = 19, mc_end_position = 27, min_cluster_size = 10):
    
    if  filename is None or output_name is None or not os.path.isfile(filename):
        sys.stderr.write("Error, one of the input file/s not present\n")
        sys.exit(-1)

    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = "."
    
    if molecular_barcodes and (mc_start_position < 0  or mc_end_position < 0 or mc_end_position <= mc_start_position):
        sys.stderr.write("Error: Molecular Barcodes option is activated but the start/end positions parameters are incorrect\n")
        sys.exit(-1)
        
    total_record = 0
    json_barcodes = list()
    json_reads = list()
    unique_genes = set()
    unique_barcodes = set()
    total_barcodes = 0
    discarded_reads = 0
    
    for transcript in parseUniqueEvents(filename):
            #re-compute the read count accounting for PCR duplicates if indicated (read sequence must contain molecular barcode)
            if molecular_barcodes:
                clusters = removePCRduplicates(transcript.sequences, allowed_missmatches, 
                                               mc_start_position, mc_end_position, min_cluster_size)
                reads_covered_by_clusters = sum(clusters)
                adjusted_reads = len(clusters) + (transcript.reads - reads_covered_by_clusters)
                discarded_reads += (transcript.reads - adjusted_reads)
                transcript.reads = adjusted_reads
            
            # add a JSON entry for the transcript  
            json_barcodes.append(transcript.toBarcodeDict())
            
            #get the reads that mapped to the transcript and generate a JSON file
            for qula, read, name in zip(transcript.qualities, transcript.sequences, transcript.readNames):
                json_reads.append({'name': str(name), 'read': str(read), 
                                   'quality': str(qula), 'barcode': transcript.barcode, 'gene': transcript.gene})
                
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
        json.dump(json_barcodes,filehandler, indent=2, separators=(',', ': '))  
    with open(os.path.join(output_folder, filenameReads), "w") as filehandlerReads:
        json.dump(json_reads,filehandlerReads, indent=2, separators=(',', ': '))    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=str,
                        help='Input file in CSV format [read_name | chromosome | gene | barcode | x | y | read_quality | read_sequence]')
    parser.add_argument('--output-folder', type=str,
                        help='Path of the output folder (default is /.)')
    parser.add_argument('--output-name', type=str,
                        help='Name of the output files')
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, help="Activates the molecular barcodes PCR duplicates filter")
    parser.add_argument('--mc-allowed-missmatches', default=2,
                        help='Number of allowed missmatches when applying the molecular barcodes PCR filter')
    parser.add_argument('--mc-start-position', default=19,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=30,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=10,
                        help='Min number of equal molecular barcodes to count as a cluster')

    args = parser.parse_args()
    main(args.input, args.output_name, args.output_folder, args.molecular_barcodes, 
         int(args.mc_allowed_missmatches), int(args.mc_start_position), int(args.mc_end_position), int(args.min_cluster_size))
                                    
