#!/usr/bin/env python
""" 
This module contains some functions and utilities for SAM/BAM files
"""

import pysam
from stpipeline.common.utils import *
import logging 
import pysam

def sortSamFile(input_sam, outputFolder=None):
    """
    :input is a SAM/BAM file with mapped reads
    :outputFolder the location where to place the output file
    It simply sorts a sam/bam file containing mapped reads by position
    """
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start SAM sorting")
    
    output_sam = replaceExtension(getCleanFileName(input_sam),'_sorted.sam')
    if outputFolder is not None and os.path.isdir(outputFolder):
        output_sam = os.path.join(outputFolder, output_sam)
        
    pysam.sort("-o", output_sam, "-O", "sam", "-T", output_sam, input_sam)
    
    if not fileOk(output_sam):
        error = "Error: output file is not present " + output_sam
        logger.error(error)
        raise RuntimeError(error + "\n")
        
    logger.info("Finish SAM sorting")
    return output_sam


def filterAnnotatedReads(annot_reads, htseq_no_ambiguous=False, 
                         outputFolder=None, keep_discarded_files=False):
    """ 
    :param annot_reads SAM file obtained from HTSEQ-Count
    :param htseq_no_ambiguous true if we want to discard ambiguous annotations
    :param outputFolder if we want to specify where to put the output file
    :param keep_discarded_files true if we want to write the un-annotated reads to a file
    This function will iterate a SAM file coming from HTSEQ-Count to discard
    un-annotated reads
    """
    
    #TODO this function should make sure that for a pair only the reverse
    #or forward are returned
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start filtering annotated reads")
    
    filter_htseq = ["__no_feature",
              "__too_low_aQual",
              "__not_aligned",
              "__alignment_not_unique"]
    
    file_output = replaceExtension(getCleanFileName(annot_reads),'_filtered.sam')
    file_output_discarded = replaceExtension(getCleanFileName(annot_reads),'_discarded.sam')
    if outputFolder is not None and os.path.isdir(outputFolder):
        file_output = os.path.join(outputFolder, file_output)
        file_output_discarded = os.path.join(outputFolder, file_output_discarded)
    
    infile = pysam.AlignmentFile(annot_reads, "r")
    outfile = pysam.AlignmentFile(file_output, "wh", template=infile)
    if keep_discarded_files:
        outfile_discarded = pysam.AlignmentFile(file_output_discarded, "wh", template=infile)
       
    dropped = 0
    for sam_record in infile:
        gene_name = str(sam_record.get_tag("XF"))
        if gene_name in filter_htseq or \
            (htseq_no_ambiguous and gene_name.find("__ambiguous") != -1):
            dropped += 1
            if keep_discarded_files:
                outfile_discarded.write(sam_record)
        else:
            outfile.write(sam_record)
                  
    infile.close()
    outfile.close()
    if keep_discarded_files:
        outfile_discarded.close()
        
    logger.info("Finish filtering annotated reads, dropped : " + str(dropped) + " reads")  
    return file_output