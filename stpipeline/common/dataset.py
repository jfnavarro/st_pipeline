""" 
This module contains routines to create
a ST dataset and some statistics. The dataset
will contain several files with the ST data in different
formats
"""
import sys
import pysam
import os
import numpy as np
from collections import defaultdict
import pandas as pd
from stpipeline.common.clustering import *
from stpipeline.common.sam_utils import parseUniqueEvents
import logging

class RangeKey:
    def __init__(self, strand, location, offset):
        self.strand = strand
        self.location = location
        self.offset = offset
    def __hash__(self):
        return hash((self.strand, self.offset))
    def __eq__(self, other):
        return self.strand == other.strand \
            and abs(self.location - other.location) <= self.offset
    
def createDataset(input_file,
                  qa_stats,
                  umi_cluster_algorithm="naive",
                  umi_allowed_mismatches=1,
                  umi_counting_offset=150,
                  output_folder=None,
                  output_template=None,
                  verbose=True):
    """
    The functions parses the reads in SAM/BAM format
    that had been annotated and demultiplexed (containing spatial
    barcode).
    It then groups them by gene-barcode to count reads. 
    It outputs the records in a matrix of counts in TSV format and BED format and it also 
    writes out some statistics.
    It will only count unique molecules using the UMI present in each read.
    :param input_file: the file with the annotated-demultiplexed records
    :param qa_stats: the Stats object to add some stats (THIS IS PASSED BY REFERENCE)
    :param umi_cluster_algorithm: the clustering algorithm to cluster UMIs
    :param umi_allowed_mismatches: the number of miss matches allowed to remove
                                  duplicates by UMIs
    :param umi_counting_offset: the number of bases allowed as offset when counting UMIs
    :param output_folder: path to place the output files
    :param output_template: the name of the dataset
    :param verbose: True if we can to collect the stats in the logger
    :type input_file: str
    :type umi_cluster_algorithm: str
    :type umi_allowed_mismatches: boolean
    :type umi_counting_offset: integer
    :type output_folder: str
    :type output_template: str
    :type verbose: bool
    :raises: RuntimeError,ValueError,OSError,CalledProcessError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(input_file):
        error = "Error creating dataset, input file not present {}\n".format(input_file)
        logger.error(error)
        raise RuntimeError(error)
      
    if output_template:
        filenameDataFrame = "{}_stdata.tsv".format(output_template)
        filenameReadsBED = "{}_reads.bed".format(output_template)
    else:
        filenameDataFrame = "stdata.tsv"
        filenameReadsBED = "reads.bed"
         
    # Some counters
    total_record = 0
    discarded_reads = 0
    
    # Obtain the clustering function
    if umi_cluster_algorithm == "naive":
        group_umi_func = countUMINaive
    elif umi_cluster_algorithm == "hierarchical":
        group_umi_func = countUMIHierarchical
    elif umi_cluster_algorithm == "Adjacent":
        group_umi_func = dedup_adj
    elif umi_cluster_algorithm == "AdjacentBi":
        group_umi_func = dedup_dir_adj
    else:
        error = "Error creating dataset.\n" \
        "Incorrect clustering algorithm {}".format(umi_cluster_algorithm)
        logger.error(error)
        raise RuntimeError(error)
 
    # Containers needed to create the data frame
    list_row_values = list()
    list_indexes = list()   
    
    # Parse unique events to generate the unique counts and the BED file    
    unique_events = parseUniqueEvents(input_file)
    with open(os.path.join(output_folder, filenameReadsBED), "w") as reads_handler:
        # Unique events is a dict() [spot][gene] -> list(transcripts)
        for (x,y), value in unique_events.iteritems():
            for gene, transcripts in value.iteritems():
                # Re-compute the read count accounting for duplicates 
                # (read sequences must contain a UMI)
                # Get the original number of transcripts (reads)
                gene_count = len(transcripts)
                # Group transcripts by start position (allowing a certain offset), UMI and strand 
                grouped_transcripts = defaultdict(lambda : defaultdict(list))
                for transcript in transcripts:
                    strand = str(transcript[5])
                    start = int(transcript[1]) if strand == "+" else int(transcript[2])
                    umi = transcript[6]
                    grouped_transcripts[RangeKey(strand, start, 
                                                 umi_counting_offset)][umi].append(transcript)
                assert len(grouped_transcripts) > 0
                # For each group of transcripts
                # cluster the transcripts based on the UMI
                # and returns the unique transcripts
                # The new gene count will be the number of unique transcripts (UMIs)
                new_gene_count = 0
                unique_transcripts = list()
                for umis in grouped_transcripts.values():
                    unique_umis = group_umi_func(umis.keys(), umi_allowed_mismatches)
                    new_gene_count += len(unique_umis)
                    # Choose 1 random transcript for the clustered transcripts (by UMI)
                    unique_transcripts += [random.choice(umis[u_umi]) for u_umi in unique_umis]
                assert new_gene_count > 0 and new_gene_count <= gene_count   
                # Update the discarded transcripts count
                discarded_reads += (gene_count - new_gene_count)
                # Update read counts in the container (replace the list
                # of transcripts for a number so it can be exported as a data frame)
                value[gene] = new_gene_count
                # Write every unique transcript to the BED output (adding spot coordinate and gene name)
                for read in unique_transcripts:
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
        error = "Error creating dataset, input file did not contain any transcript\n"
        logger.error(error)
        raise RuntimeError(error)
    
    # Create the data frame
    counts_table = pd.DataFrame(list_row_values, index=list_indexes)
    counts_table.fillna(0, inplace=True)
    
    # Compute some statistics
    total_barcodes = len(counts_table.index)
    total_transcripts = np.sum(counts_table.values, dtype=np.int32)
    number_genes = len(counts_table.columns)
    max_count = counts_table.values.max()
    min_count = counts_table.values.min()
    aggregated_spot_counts = counts_table.sum(axis=1).values
    aggregated_gene_counts = (counts_table != 0).sum(axis=1).values
    max_genes_feature = aggregated_gene_counts.max()
    min_genes_feature = aggregated_gene_counts.min()
    max_reads_feature = aggregated_spot_counts.max()
    min_reads_feature = aggregated_spot_counts.min()
    average_reads_feature = np.mean(aggregated_spot_counts)
    average_genes_feature = np.mean(aggregated_gene_counts)
    std_reads_feature = np.std(aggregated_spot_counts)
    std_genes_feature = np.std(aggregated_gene_counts)
        
    # Print some statistics
    if verbose:
        logger.info("Number of unique molecules present: {}".format(total_transcripts))
        logger.info("Number of unique events (gene-feature) present: {}".format(total_record))
        logger.info("Number of unique genes present: {}".format(number_genes))
        logger.info("Max number of genes over all features: {}".format(max_genes_feature))
        logger.info("Min number of genes over all features: {}".format(min_genes_feature))
        logger.info("Max number of unique molecules over all features: {}".format(max_reads_feature))
        logger.info("Min number of unique molecules over all features: {}".format(min_reads_feature))
        logger.info("Average number genes per feature: {}".format(average_genes_feature))
        logger.info("Average number unique molecules per feature: {}".format(average_reads_feature))
        logger.info("Std number genes per feature: {}".format(std_genes_feature))
        logger.info("Std number unique molecules per feature: {}".format(std_reads_feature))
        logger.info("Max number of unique molecules over all unique events: {}".format(max_count))
        logger.info("Min number of unique molecules over all unique events: {}".format(min_count))
        logger.info("Number of discarded reads (possible duplicates): {}".format(discarded_reads))
        
    # Update the QA object
    qa_stats.reads_after_duplicates_removal = (total_transcripts - discarded_reads)
    qa_stats.unique_events = total_record
    qa_stats.barcodes_found = total_barcodes
    qa_stats.genes_found = number_genes
    qa_stats.duplicates_found = discarded_reads
    qa_stats.max_genes_feature = max_genes_feature
    qa_stats.min_genes_feature = min_genes_feature
    qa_stats.max_reads_feature = max_reads_feature
    qa_stats.min_reads_feature = min_reads_feature
    qa_stats.max_reads_unique_event = max_count
    qa_stats.min_reads_unique_event = min_count
    qa_stats.avergage_gene_feature = average_genes_feature
    qa_stats.average_reads_feature = average_reads_feature
     
    # Write data frame to file
    counts_table.to_csv(os.path.join(output_folder, filenameDataFrame), sep="\t", na_rep=0)       
