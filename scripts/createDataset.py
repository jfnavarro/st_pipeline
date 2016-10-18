#! /usr/bin/env python
""" 
Scripts that parses a SAM or BAM file containing
aligned/annotated reads and the spot coordinates,
gene and UMI (X,Y,gene,UMI) as Tags (BZ).

It then generates a data frame (genes as columns) 
with the ST data and a BED file with the unique transcripts. 
It does so by aggregating the reads to compute counts and it uses
for each gene-spot. 

It removes duplicates by using Unique Molecular Identifiers.
To enable this you must activate the flag and indicate
what type of algorithm to use for clustering and how many
similar UMIs must be to be counted as one unique transcript
(for this you can pass the number of mismatches).

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se> 
"""
import sys
import argparse
import pysam
import os
import numpy as np
from collections import defaultdict
import pandas as pd
from stpipeline.common.clustering import countMolecularBarcodesClustersHierarchical, countMolecularBarcodesClustersNaive

#TODO this function uses too much memory, optimize it. (Maybe Cython)
def parseUniqueEvents(filename, molecular_barcodes=False):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM/BAM file where the spot coordinates are present in the tags
    The output will be hash table [spot][gene] -> list of transcripts. 
    :param filename: the input file containing the SAM/BAM records
    :param molecular_barcodes: if True a the records contain UMIs
    :return: A map of spots(x,y) to a map of gene names to a list of transcript 
    (chrom, start, end, clear_name, mapping_quality, strand, sequence)
    As map[(x,y)][gene]->list((chrom, start, end, clear_name, mapping_quality, strand, sequence))
    """
    
    unique_events = defaultdict(lambda : defaultdict(list))
    sam_type = os.path.splitext(filename)[1].lower()
    flag = "r" if sam_type == ".sam" else "rb"
    sam_file = pysam.AlignmentFile(filename, flag)
    
    for rec in sam_file.fetch(until_eof=True):
        clear_name = rec.query_name
        mapping_quality = rec.mapping_quality
        start = rec.reference_start
        end = rec.reference_end
        chrom = sam_file.getrname(rec.reference_id)
        strand = "-" if rec.is_reverse else "+"
        # Get taggd tags
        x,y,gene,seq = (None,None,None,None)
        for (k, v) in rec.tags:
            if k == "B1":
                x = int(v) ## The X coordinate
            elif k == "B2":
                y = int(v) ## The Y coordinate
            elif k == "XF":
                gene = str(v) ## The gene name
            elif k == "B3":
                seq = str(v) ## The UMI (optional)
            else:
                continue
        # Check that all tags are present
        if any(tag is None for tag in [x,y,gene]) or (molecular_barcodes and not seq):
            sys.stdout.write("Warning: Missing attributes for record {}\n".format(clear_name))
            continue
        
        # Create a new transcript and add it to the dictionary
        transcript = (chrom, start, end, clear_name, mapping_quality, strand, seq)
        unique_events[(x,y)][gene].append(transcript)
      
    sam_file.close()
    return unique_events

def main(filename, 
         output_folder,
         output_file_template,
         molecular_barcodes,
         mc_cluster,
         allowed_mismatches,
         min_cluster_size):
    
    if filename is None or not os.path.isfile(filename):
        sys.stderr.write("Error, input file not present or invalid: {}\n".format(filename))
        sys.exit(1)

    sam_type = os.path.splitext(filename)[1].lower()
    if sam_type not in [".sam",".bam"]:
        sys.stderr.write("Error, invalid input format: {}\n".format(sam_type))
        sys.exit(1)
        
    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = os.getcwd()
    
    if mc_cluster not in ["naive","hierarchical"]:
        sys.stderr.write("Error, type of clustering algorithm is incorrect\n")
        sys.exit(1)
        
    if output_file_template:
        filenameDataFrame = "{}_stdata.tsv".format(output_file_template)
        filenameReadsBED = "{}_reads.bed".format(output_file_template)
    else:
        filenameDataFrame = "stdata.tsv"
        filenameReadsBED = "reads.bed"
        
    # Some counters
    total_record = 0
    discarded_reads = 0
    
    # Obtain the clustering function
    if mc_cluster == "naive":
        clusters_func = countMolecularBarcodesClustersNaive
    else:
        clusters_func = countMolecularBarcodesClustersHierarchical
        
    # Containers needed to create the data frame
    list_row_values = list()
    list_indexes = list()   
    
    # Parse unique events to generate the unique counts and the BED file    
    unique_events = parseUniqueEvents(filename, molecular_barcodes)
    with open(os.path.join(output_folder, filenameReadsBED), "w") as reads_handler:
        # Unique events is a dict() [spot][gene] -> list(transcripts)
        for spot, value in unique_events.iteritems():
            x = spot[0]
            y = spot[1]
            for gene, transcripts in value.iteritems():
                # Re-compute the read count accounting for duplicates 
                # (read sequences must contain a molecular barcode)
                if molecular_barcodes:
                    # New list of unique transcripts
                    unique_transcripts = list()
                    # Get the original number of transcripts
                    num_transcripts = len(transcripts)
                    # Group transcripts by start pos and strand
                    grouped_transcripts = defaultdict(list)
                    for transcript in transcripts:
                        grouped_transcripts[(transcript[1], transcript[5])].append((transcript[6],transcript)) 
                    # For each group of transcripts
                    # cluster the transcripts based on the molecular barcode
                    # and returns the unique transcripts
                    for mcs in grouped_transcripts.values():
                        unique_transcripts += clusters_func(mcs,
                                                            allowed_mismatches,
                                                            min_cluster_size)
                    # Update the discarded transcripts count
                    transcripts = unique_transcripts
                    discarded_reads += (num_transcripts - len(transcripts))
                    
                # Update read counts in the container (replace the list
                # of transcripts for a number so it can be exported as a data frame)
                value[gene] = len(transcripts)
             
                # Write every unique transcript to the BED output (adding spot coordinate and gene name)
                for read in transcripts:
                    reads_handler.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(read[0],
                                                                                               read[1],
                                                                                               read[2],
                                                                                               read[3],
                                                                                               read[4],
                                                                                               read[5],
                                                                                               gene,
                                                                                               x,y)) 
                # keep a counter of the number of unique events (spot - gene) processed
                total_record += 1
                
            # Add spot and dict [gene] -> count to containers
            list_indexes.append("{0}x{1}".format(x, y))
            list_row_values.append(value)
            
    if total_record == 0:
        sys.stderr.write("Error, the number of transcripts present is 0\n")
        sys.exit(1)
    
    # Create the data frame
    counts_table = pd.DataFrame(list_row_values, index=list_indexes)
    counts_table.fillna(0, inplace=True)
    
    # Compute some statistics
    total_transcripts = np.sum(counts_table.values, dtype=np.int32)
    number_genes = len(counts_table.columns)
    max_count = counts_table.values.max()
    min_count = counts_table.values.min()
    aggregated_spot_counts = counts_table.sum(axis=0).values
    aggregated_gene_counts = counts_table.sum(axis=1).values
    
    # Print some statistics
    print "Number of unique transcripts present: {}".format(total_transcripts)
    print "Number of unique events (gene-barcode) present: {}".format(total_record)
    print "Number of unique genes present: {}".format(number_genes)
    # TODO I really have no idea why this does not work {:f} would not work either
    #print "Max number of genes over all features: {}".format(aggregated_gene_counts.max())
    #print "Min number of genes over all features: {}".format(aggregated_gene_counts.min())
    #print "Max number of reads over all features: {}".format(aggregated_spot_counts.max())
    #print "Min number of reads over all features: {}".format(aggregated_spot_counts.min())
    #print "Max number of reads over all unique events: {}".format(max_count)
    #print "Min number of reads over all unique events: {}".format(min_count)
    if molecular_barcodes:
        print "Number of discarded reads (possible PCR duplicates): {}".format(discarded_reads)
     
    # Write data frame to file
    counts_table.to_csv(os.path.join(output_folder,filenameDataFrame), sep="\t", na_rep=0)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', type=str, required=True,
                        help="Input file in SAM or BAM format containing annotated " \
                        "records with the spot coordinates and the gene as tags")
    parser.add_argument('--output-folder', type=str, default=None,
                        help="Path of the output folder (default is /.)")
    parser.add_argument('--output-file-template', type=str, default=None,
                        help="Describes how the output files will be called")
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", default=False, 
                        help="Activates the UMIs duplicates filter")
    parser.add_argument('--mc-cluster', default="naive", metavar="[STR]", 
                        type=str, choices=["naive", "hierarchical"],
                        help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                        "Modes = {naive(default), hierarchical}")
    parser.add_argument('--mc-allowed-mismatches', default=1,  metavar="[INT]", type=int, choices=range(1,10),
                        help='Number of allowed mismatches when applying the UMIs filter (default: %(default)s)')
    parser.add_argument('--min-cluster-size', default=2,  metavar="[INT]", type=int, choices=range(1, 100),
                        help='Min number of equal UMIs to count as an unique transcript (default: %(default)s)')

    args = parser.parse_args()
    main(args.input, 
         args.output_folder, 
         args.output_file_template, 
         args.molecular_barcodes, 
         args.mc_cluster,
         int(args.mc_allowed_mismatches), 
         int(args.min_cluster_size))
                                    
