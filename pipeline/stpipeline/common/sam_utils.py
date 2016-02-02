#!/usr/bin/env python
""" 
This module contains some functions and utilities for SAM/BAM files
"""

from stpipeline.common.utils import *
from stpipeline.common.fastq_utils import reverse_complement
import logging 
import pysam
from stpipeline.common.stats import Stats
from HTSeq import BAM_Reader, SAM_Reader, BAM_Writer

def sortSamFile(input_sam, outputFolder=None):
    """
    :param input is a SAM/BAM file with mapped reads
    :param outputFolder the location where to place the output file
    It simply sorts a sam/bam file containing mapped reads by position
    """
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start SAM sorting")
    
    sam_type = getExtension(input_sam).lower()
    output_sam = 'mapped_filtered_sorted.' + sam_type
        
    if outputFolder is not None and os.path.isdir(outputFolder):
        output_sam = os.path.join(outputFolder, output_sam)
        
    pysam.sort("-n", "-o", output_sam, "-O", sam_type, "-T", output_sam, input_sam)
    
    if not fileOk(output_sam):
        error = "Error annotating: output file is not present " + output_sam
        logger.error(error)
        raise RuntimeError(error + "\n")
        
    logger.info("Finish SAM sorting")
    return output_sam

def filterMappedReads(mapped_reads, 
                      qa_stats, 
                      hash_reads,
                      min_length=28, 
                      pair_mode_keep="reverse",
                      outputFolder=None, 
                      keep_discarded_files=False):
    """ 
    :param mapped_reads SAM file obtained from STAR containing the aligmnents
    :param hash_reads a hash table with read names that contain barcode
    :param min_length the min number of mapped bases we enforce in an alignment
    :param pair_mode_keep if both pairs are aligned what to do (reverse, forward, both)
    :param outputFolder if we want to specify where to put the output file
    :param keep_discarded_files true if we want to write the un-annotated reads to a file
    Iterate the alignments and discards un-mapped and if both strands
    map, keep only the reverse
    """
    
    logger = logging.getLogger("STPipeline")
    assert(pair_mode_keep in ["reverse", "forward", "both"])
    
    sam_type = getExtension(mapped_reads).lower()
    file_output = 'mapped_filtered.' + sam_type
    file_output_discarded = 'mapped_discarded.' + sam_type
        
    if outputFolder is not None and os.path.isdir(outputFolder):
        file_output = os.path.join(outputFolder, file_output)
        file_output_discarded = os.path.join(outputFolder, file_output_discarded)
    
    flag_read = "rb"
    flag_write = "wb"
    if sam_type == "sam":
        flag_read = "r"
        flag_write = "wh"
        
    infile = pysam.AlignmentFile(mapped_reads, flag_read)
    outfile = pysam.AlignmentFile(file_output, flag_write, template=infile)
    if keep_discarded_files:
        outfile_discarded = pysam.AlignmentFile(file_output_discarded, flag_write, template=infile)
       
    keep_only_forward = (pair_mode_keep == "forward")
    keep_only_reverse = (pair_mode_keep == "reverse")
    dropped_unmapped = 0
    dropped_secondary = 0
    dropped_short = 0
    dropped_both_pairs = 0
    dropped_barcode = 0
    present = 0
    
    # to remove secondary alignments use sam_record.is_secondary
    for sam_record in infile.fetch(until_eof=True):
        present += 1
        discard_read = False
        
        # Add the barcode and coordinates info if present otherwise discard
        if hash_reads.has_key(sam_record.query_name):
            for tag in hash_reads[sam_record.query_name]:
                tag_tokens = tag.split(":")
                sam_record.set_tag(tag_tokens[0], tag_tokens[2], tag_tokens[1])
        else:
            dropped_barcode += 1
            continue
        
        # Reverse the read if It mapped to the reverse complement
        if sam_record.is_reverse:
            qual = sam_record.query_qualities
            sam_record.query_sequence = reverse_complement(sam_record.query_sequence)
            # Apparently updating the sequence resets the quality scores
            sam_record.query_qualities = qual
            
        # Get how many bases were mapped
        mapped_bases = 0
        for cigar_tuple in sam_record.cigartuples:
            if cigar_tuple[0] == 0:
                mapped_bases += cigar_tuple[1]
        
        # We need this so we don't duplicate reads
        if not sam_record.is_secondary:
            sam_record.set_tag("NH", None)
            
        # Collect stats for different types of reads
        # and discard when necessary  
        if sam_record.is_unmapped:
            dropped_unmapped += 1
            discard_read = True
        elif sam_record.is_secondary:
            dropped_secondary += 1
            discard_read = True
        elif mapped_bases != 0 and mapped_bases < min_length:
            dropped_short += 1
            discard_read = True
        elif (sam_record.is_paired and sam_record.is_proper_pair and \
         ((sam_record.is_read1 and keep_only_reverse) or \
          (sam_record.is_read2 and keep_only_forward))):
            dropped_both_pairs += 1
            discard_read = True
            
        if discard_read:
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
            
    logger.info("Finish filtering mapped reads, stats:" \
                "\nPresent: " + str(present) + \
                "\nDropped - unmapped : " + str(dropped_unmapped) + \
                "\nDropped - secondary alignment : " + str(dropped_secondary) + \
                "\nDropped - too short : " + str(dropped_short) + \
                "\nDropped - barcode : " + str(dropped_barcode) + \
                "\nDropped - two pair aligned : " + str(dropped_both_pairs))  
    
    # Update QA object 
    qa_stats.reads_after_mapping = present - (dropped_unmapped + dropped_secondary + 
                                              dropped_short + dropped_both_pairs)
    return file_output