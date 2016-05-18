#!/usr/bin/env python
""" 
This module contains routines
to compute saturation points on a
set of annotated reads in BAM/SAM format
"""
import pysam
import random
import math
import os
from collections import defaultdict
import logging
from stpipeline.common.dataset import createDataset
from stpipeline.common.stats import Stats
from stpipeline.common.utils import safeRemove, getExtension

def computeSaturation(nreads, 
                      annotated_reads,
                      molecular_barcodes,
                      mc_cluster,
                      mc_allowed_mismatches,
                      min_cluster_size,
                      expName,
                      low_memory,
                      temp_folder=None):
    """
    It splits the input file up into sub-files containing
    random reads from the input file up to the saturation point. 
    It then calls createDataset.py and retrieve some stats to
    compute the saturation points information that is then added
    to the log file.
    :param nreads: the number of reads present in the annotated_reads file
    :param annotated_reads: path to a SAM/BAM file with the annotated reads
    :param molecular_barcodes: True is the reads contain UMIs
    :param mc_allowed_mismatches: the number of miss matches allowed to remove
                                  duplicates by UMIs
    :param min_cluster_size: the min size of the clusters to remove duplicates by UMIs
    :param expName: the name of the dataset
    :param low_memory: True if we want to use a slower but more memory efficient approach
    :param temp_folder: the path where to put the output files
    :type nreads: integer
    :type annotated_reads: str
    :type molecular_barcodes: boolean
    :type mc_allowed_mismatches: boolean
    :type min_cluster_size: integer
    :type expName: str
    :type low_memory: boolean
    :type temp_folder: str
    :raises: RuntimeError
    """
      
    logger = logging.getLogger("STPipeline")
              
    saturation_points = list()
    for x in xrange(0,15):
        spoint = int(math.floor(1e5 + (math.exp(x) * 1e5)))
        if spoint >= nreads:
            break
        saturation_points.append(spoint)
        
    files = dict()
    file_names = dict()
    subsampling = dict()
    file_ext = getExtension(annotated_reads).lower()
    flag_read = "rb"
    flag_write = "wb"
    if file_ext == "sam":
        flag_read = "r"
        flag_write = "wh"
                 
    annotated_sam = pysam.AlignmentFile(annotated_reads, flag_read)
             
    # Generate subsamples and files
    for spoint in saturation_points:
        # Create a file for the sub sample point
        file_name = "subsample_" + str(spoint) + "." + file_ext
        if temp_folder is not None and os.path.isdir(temp_folder):
            file_name = os.path.join(temp_folder, file_name)
        output_sam = pysam.AlignmentFile(file_name, flag_write, template=annotated_sam)
        file_names[spoint] = file_name
        files[spoint] = output_sam
        # Generate a list of indexes in the sam file to extract sub samples 
        indices = list(xrange(nreads))
        random.shuffle(indices)
        subbed = indices[0:spoint]
        subbed.sort()
        subsampling[spoint] = subbed
                 
    # Write subsamples
    index = 0
    sub_indexes = defaultdict(int)
    for read in annotated_sam.fetch(until_eof=True):
        for spoint in saturation_points:
            sub_index = sub_indexes[spoint]
            if sub_index < len(subsampling[spoint]) and subsampling[spoint][sub_index] == index:
                files[spoint].write(read)
                sub_indexes[spoint] += 1
        index += 1  
                 
    # Close the files
    annotated_sam.close()
    for file_sam in files.itervalues():
        file_sam.close()
                 
    # Compute saturation points
    saturation_points_values_unique_events = list()
    saturation_points_values_reads = list()
    saturation_points_values_genes = list()
    for spoint in saturation_points:
        stats = Stats()
        input_file = file_names[spoint]
        try:
            createDataset(input_file,
                          stats,
                          molecular_barcodes,
                          mc_cluster,
                          mc_allowed_mismatches,
                          min_cluster_size,
                          expName,
                          False, # Verbose
                          low_memory)
        except Exception:
            error = "Error computing saturation curve: createDataset execution failed\n"
            logger.error(error)
            print "Error", error
            raise
            
        saturation_points_values_unique_events.append(stats.unique_events)
        saturation_points_values_reads.append(stats.reads_after_duplicates_removal)
        saturation_points_values_genes.append(stats.genes_found)
                     
    # Remove the files
    for file_sam in file_names.itervalues():
        safeRemove(file_sam)
                         
    # Generate points for the plots
    logger.info("Saturation points:")
    logger.info(', '.join(str(a) for a in saturation_points))
    logger.info("Unique events per saturation point")
    logger.info(', '.join(str(a) for a in saturation_points_values_unique_events))
    logger.info("Unique transcripts per saturation point")
    logger.info(', '.join(str(a) for a in saturation_points_values_reads))
    logger.info("Unique genes per saturation point")
    logger.info(', '.join(str(a) for a in saturation_points_values_genes))
        