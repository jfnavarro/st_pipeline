""" 
This module contains routines to create
a ST dataset and some statistics. The dataset
will contain several files with the ST data in different
formats
"""
import sys
import os
import numpy as np
from collections import defaultdict
import pandas as pd
from stpipeline.common.clustering import *
from stpipeline.common.unique_events_parser import parse_unique_events
import logging
import sys

def computeUniqueUMIs(transcripts, umi_counting_offset, umi_allowed_mismatches, group_umi_func):
    """ Helper function to compute unique transcripts UMIs from
    a given list of transcripts
    """
    # Sort transcripts by strand and start position
    sorted_transcripts = sorted(transcripts, key = lambda x: (x[5], x[1]))
    # Group transcripts by strand and start-position allowing an offset
    # And then performs the UMI clustering in each group to finally
    # compute the gene count as the sum of the unique UMIs for each group (strand,start,offset)
    grouped_transcripts = defaultdict(list)
    # TODO A probably better approach is to get the mean of all the start positions
    # and then make mean +- 300bp (user defined) a group to account for the library
    # size variability and then group the rest of transcripts normally by (strand, start, position).
    unique_transcripts = list()
    num_transcripts = len(transcripts)
    for i in range(num_transcripts - 1):
        current = sorted_transcripts[i]
        nextone = sorted_transcripts[i + 1]
        grouped_transcripts[current[6]].append(current)
        if abs(current[1] - nextone[1]) > umi_counting_offset or current[5] != nextone[5]:
            # A new group has been reached (strand, start-pos, offset)
            # Compute unique UMIs by hamming distance
            unique_umis = group_umi_func(list(grouped_transcripts.keys()), umi_allowed_mismatches)
            # Choose 1 random transcript for the clustered transcripts (by UMI)
            unique_transcripts += [random.choice(grouped_transcripts[u_umi]) for u_umi in unique_umis]
            # Reset the container
            grouped_transcripts = defaultdict(list)
    # We process the last one and more transcripts if they were not processed
    lastone = sorted_transcripts[num_transcripts - 1]
    grouped_transcripts[lastone[6]].append(lastone)
    unique_umis = group_umi_func(list(grouped_transcripts.keys()), umi_allowed_mismatches)
    unique_transcripts += [random.choice(grouped_transcripts[u_umi]) for u_umi in unique_umis]
    return unique_transcripts

def createDataset(input_file,
                  qa_stats,
                  gff_filename=None,
                  umi_cluster_algorithm="hierarchical",
                  umi_allowed_mismatches=1,
                  umi_counting_offset=250,
                  diable_umi=False,
                  output_folder=None,
                  output_template=None,
                  verbose=True):
    """
    The functions parses the reads in BAM format
    that have been annotated and demultiplexed (containing spatial barcode).
    It then groups them by gene-barcode to count reads accounting for duplicates
    using the UMIs (clustering them suing the strand and start position). 
    It outputs the records in a matrix of counts in TSV format and BED format and it also 
    writes out some statistics.
    :param input_file: the file with the annotated-demultiplexed records in BAM format
    :param qa_stats: the Stats object to add some stats (THIS IS PASSED BY REFERENCE)
    :param umi_cluster_algorithm: the clustering algorithm to cluster UMIs
    :param umi_allowed_mismatches: the number of miss matches allowed to remove
                                  duplicates by UMIs
    :param umi_counting_offset: the number of bases allowed as offset (start position) when counting UMIs
    :param diable_umi: when True the UMI filtering step will not be performed
    :param output_folder: path to place the output files
    :param output_template: the name of the dataset
    :param verbose: True if we can to collect the stats in the logger
    :type input_file: str
    :type umi_cluster_algorithm: str
    :type umi_allowed_mismatches: boolean
    :type umi_counting_offset: integer
    :type diable_umi: bool
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
    elif umi_cluster_algorithm == "Affinity":
        group_umi_func = affinity_umi_removal
    else:
        error = "Error creating dataset.\n" \
        "Incorrect clustering algorithm {}".format(umi_cluster_algorithm)
        logger.error(error)
        raise RuntimeError(error)
 
    # Containers needed to create the data frame
    list_row_values = list()
    list_indexes = list()   

    # Parse unique events to generate the unique counts and the BED file
    unique_events = parse_unique_events(input_file, gff_filename)
    with open(os.path.join(output_folder, filenameReadsBED), "w") as reads_handler:
        # this is the generator returning a dictionary with spots for each gene
        for gene, spots in unique_events:
            transcript_counts_by_spot = {}
            for spot_coordinates, reads in list(spots.items()):
                (x,y) = spot_coordinates
                # Re-compute the read count accounting for duplicates using the UMIs
                # Transcripts is the list of transcripts (chrom, start, end, clear_name, mapping_quality, strand, UMI)
                # First:
                # Get the original number of transcripts (reads)
                reads = list(reads)
                read_count = len(reads)
                if not diable_umi:
                    # Compute unique transcripts (based on UMI, strand and start position +- threshold)
                    unique_transcripts = computeUniqueUMIs(reads, umi_counting_offset, 
                                                           umi_allowed_mismatches, group_umi_func)
                else:
                    unique_transcripts = reads
                # The new transcript count
                transcript_count = len(unique_transcripts)
                assert transcript_count > 0 and transcript_count <= read_count
                # Update the discarded reads count
                discarded_reads += (read_count - transcript_count)
                # Update read counts in the container (replace the list
                # of transcripts for a number so it can be exported as a data frame)
                transcript_counts_by_spot["{0}x{1}".format(x, y)] = transcript_count
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
            list_indexes.append(gene)
            list_row_values.append(transcript_counts_by_spot)
            
    if total_record == 0:
        error = "Error creating dataset, input file did not contain any transcript\n"
        logger.error(error)
        raise RuntimeError(error)
    
    # Create the data frame
    counts_table = pd.DataFrame(list_row_values, index=list_indexes)
    counts_table.fillna(0, inplace=True)
    counts_table=counts_table.T # Transpose the dictionary to still get the spots as rows and genes as columns in the final tsv
    
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
    qa_stats.reads_after_duplicates_removal = int(total_transcripts)
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
    qa_stats.average_gene_feature = average_genes_feature
    qa_stats.average_reads_feature = average_reads_feature
     
    # Write data frame to file
    counts_table.to_csv(os.path.join(output_folder, filenameDataFrame), sep="\t", na_rep=0)       
