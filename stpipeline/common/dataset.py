#!/usr/bin/env python
""" 
This module contains routines to create
a ST dataset and some statistics. The dataset
will contain several files with the ST data in different
formats
"""
import subprocess
import gc
import logging

def createDataset(input_file,
                  qa_stats,
                  molecular_barcodes=False,
                  mc_cluster="naive",
                  allowed_mismatches=1,
                  min_cluster_size=2,
                  output_folder=None,
                  output_template=None,
                  verbose=True,
                  low_memory=False):
    """
    The script createDataset.py parses reads in SAM/BAM format
    that had been annotated and demultiplexed (containing spatial
    barcode).
    It then groups them by gene-barcode to count reads. 
    It outputs the records in JSON format and BED format and it also 
    writes out some statistics.
    It also allows to remove PCR Duplicates using molecular barcodes
    This function is a wrapper around createdDataset.py
    :param input_file: the file with the annotated-demultiplexed records
    :param qa_stats: the Stats() object to store different stats
    :param molecular_barcodes: True if the reads contain UMIs
    :param mc_cluster: the type of clustering (naive or hierarchical)
    :param allowed_mismatches: how many mismatches allowed when clustering UMIS
    :param min_cluster_size: the min size of a cluster when clustering UMIs
    :param output_folder: path to place the output files
    :param output_template: the name of the dataset
    :param verbose: True if we can to collect the stats in the logger
    :param low_memory: True if we want to run a slower but more memory efficient algorithm
    :type input_file: str
    :type molecular_barcodes: bool
    :type allowed_mismatches: integer
    :type min_cluster_size: integer
    :type output_folder: str
    :type output_template: str
    :type verbose: bool
    :type low_memory: bool
    :raises: RuntimeError
    """
    logger = logging.getLogger("STPipeline")
    
    args = ['createDataset.py', '--input', str(input_file)]
        
    if molecular_barcodes:
        args += ['--molecular-barcodes',
                '--mc-allowed-mismatches', allowed_mismatches,
                '--min-cluster-size', min_cluster_size,
                '--mc-cluster', mc_cluster]
            
    if output_folder: args += ['--output-folder', output_folder]
    if output_template: args += ['--output-file-template', output_template]
    if low_memory: args += ['--low-memory']
        
    gc.collect()      
    try:
        proc = subprocess.Popen([str(i) for i in args], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                shell=False, close_fds=True)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        error = "Error creating dataset: createDataset execution failed\n"
        logger.error(error)
        logger.error(e)
        raise RuntimeError(error) 
        
    if len(errmsg) > 0:
        error = "Error, There was an error creating the dataset (%s\n%s)" % (stdout,errmsg)
        logger.error(error)
        raise RuntimeError(error)    
              
    procOut = stdout.split("\n")
    for line in procOut:
        # Collect and write QA stats
        # TODO find a cleaner way to do this
        if line.find("Number of unique transcripts present:") != -1:
            qa_stats.reads_after_duplicates_removal = int(line.split()[-1])
        if line.find("Number of unique events (gene-barcode) present:") != -1:
            qa_stats.unique_events = int(line.split()[-1])
        if line.find("Number of unique barcodes present:") != -1:
            qa_stats.barcodes_found = int(line.split()[-1])
        if line.find("Number of unique genes present:") != -1:
            qa_stats.genes_found = int(line.split()[-1])
        if line.find("Number of discarded reads (possible PCR duplicates):") != -1:
            qa_stats.duplicates_found = int(line.split()[-1])
        if line.find("Max number of genes over all features:") != -1:
            qa_stats.max_genes_feature = int(line.split()[-1])
        if line.find("Min number of genes over all features:") != -1:
            qa_stats.min_genes_feature = int(line.split()[-1])
        if line.find("Max number of reads over all features:") != -1:
            qa_stats.max_reads_feature = int(line.split()[-1])
        if line.find("Min number of reads over all features:") != -1:
            qa_stats.min_reads_feature = int(line.split()[-1])
        if line.find("Max number of reads over all unique events:") != -1:
            qa_stats.max_reads_unique_event = int(line.split()[-1])
        if line.find("Min number of reads over all unique events:") != -1:
            qa_stats.min_reads_unique_event = int(line.split()[-1])
        if verbose:   
            logger.info(str(line))