#!/usr/bin/env python
""" 
This module contains some functions and utilities for SAM/BAM files
"""

from stpipeline.common.utils import *
from stpipeline.common.fastq_utils import reverse_complement
import logging 
import pysam

def sortSamFile(input_sam, outputFolder=None):
    """
    :param input is a SAM/BAM file with mapped reads
    :param outputFolder the location where to place the output file
    It simply sorts a sam/bam file containing mapped reads by position
    """
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start SAM sorting")
    
    output_sam = 'mapped_filtered_sorted.sam'
    if outputFolder is not None and os.path.isdir(outputFolder):
        output_sam = os.path.join(outputFolder, output_sam)
        
    pysam.sort("-n", "-o", output_sam, "-O", "sam", "-T", output_sam, input_sam)
    
    if not fileOk(output_sam):
        error = "Error annotating: output file is not present " + output_sam
        logger.error(error)
        raise RuntimeError(error + "\n")
        
    logger.info("Finish SAM sorting")
    return output_sam

def filterMappedReads(mapped_reads, min_length=28, pair_mode_keep="reverse",
                      outputFolder=None, keep_discarded_files=False):
    """ 
    :param annot_reads SAM file obtained from STAR
    :param action to perform if both strands are mapped (keep reverse, forward or both)
    :param outputFolder if we want to specify where to put the output file
    :param keep_discarded_files true if we want to write the un-annotated reads to a file
    Iterate the alignments and discards un-mapped and if both strands
    map, keep only the reverse
    """
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start filtering mapped reads")

    assert(pair_mode_keep in ["reverse", "forward", "both"])
    
    file_output = 'mapped_filtered.sam'
    file_output_discarded = 'mapped_discarded.sam'
    if outputFolder is not None and os.path.isdir(outputFolder):
        file_output = os.path.join(outputFolder, file_output)
        file_output_discarded = os.path.join(outputFolder, file_output_discarded)
    
    infile = pysam.AlignmentFile(mapped_reads, "r")
    outfile = pysam.AlignmentFile(file_output, "wh", template=infile)
    if keep_discarded_files:
        outfile_discarded = pysam.AlignmentFile(file_output_discarded, "wh", template=infile)
       
    keep_only_forward = (pair_mode_keep == "forward")
    keep_only_reverse = (pair_mode_keep == "reverse")
    dropped = 0
    #to remove secondary alignments sam_record.is_secondary
    for sam_record in infile:
        if sam_record.is_reverse:
            sam_record.query_sequence = reverse_complement(sam_record.query_sequence)
           
        mapped_bases = 0
        for cigar_tuple in sam_record.cigartuples:
            if cigar_tuple[0] == 0:
                mapped_bases += cigar_tuple[1]
           
        if not sam_record.is_secondary:
            sam_record.set_tag("NH", None)
                 
        if sam_record.is_unmapped or \
        sam_record.is_secondary or \
        mapped_bases < min_length or \
        (sam_record.is_paired and sam_record.is_proper_pair and \
         ((sam_record.is_read1 and keep_only_reverse) or (sam_record.is_read2 and keep_only_forward))):
            dropped += 1
            if keep_discarded_files:
                outfile_discarded.write(sam_record)
        else:
            outfile.write(sam_record)
                  
    infile.close()
    outfile.close()
    if keep_discarded_files:
        outfile_discarded.close()

    if not fileOk(file_output):
        error = "Error filtering mapped reads: output file is not present " + file_output
        logger.error(error)
        raise RuntimeError(error + "\n")
            
    logger.info("Finish filtering mapped reads, dropped : " + str(dropped) + " reads")  
    return file_output

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
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start filtering annotated reads")
    
    filter_htseq = ["__no_feature",
              "__too_low_aQual",
              "__not_aligned",
              "__alignment_not_unique"]
    
    file_output = 'annotated_filtered.sam'
    file_output_discarded = 'annotated_discarded.sam'
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
        if gene_name in filter_htseq or sam_record.is_unmapped or \
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

    if not fileOk(file_output):
        error = "Error filtering annotated reads: output file is not present " + file_output
        logger.error(error)
        raise RuntimeError(error + "\n")
            
    logger.info("Finish filtering annotated reads, dropped : " + str(dropped) + " reads")  
    return file_output