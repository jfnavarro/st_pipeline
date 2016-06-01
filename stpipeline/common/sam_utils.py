""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import os
import logging 
import pysam

def sortSamFile(input_sam, outputFolder=None):
    """
    It simply sorts by position a sam/bam file containing mapped reads 
    :param input: is a SAM/BAM file with mapped reads
    :param outputFolder: the location where to place the output file (optional)
    :type input: str
    :type outputFolder: str
    :returns: the path to the sorted file
    :raises: RuntimeError
    """
    
    logger = logging.getLogger("STPipeline")
    
    sam_type = os.path.splitext(input_sam)[1].lower()
    output_sam = 'mapped_filtered_sorted{}'.format(sam_type)
        
    if outputFolder is not None and os.path.isdir(outputFolder):
        output_sam = os.path.join(outputFolder, output_sam)
        
    pysam.sort("-n", "-o", output_sam, "-O", sam_type, 
               "-T", output_sam, input_sam)
    
    if not fileOk(output_sam):
        error = "Error sorting the SAM/BAM file.\n" \
        "Output file is not present\n {}".format(output_sam)
        logger.error(error)
        raise RuntimeError(error)
        
    return output_sam

def filterMappedReads(mapped_reads,
                      hash_reads,
                      file_output,
                      file_output_discarded=None,
                      min_length=28):
    """ 
    Iterate a SAM/BAM file containing mapped reads 
    and discards reads that are secondary or too short.
    It also discards reads that do not contain a valid barcode.
    It will add the barcode, coordinates and umi as extra tags
    to the output SAM/BAM file. The UMI will be added only if it is present.
    It assumes all the reads are mapped (do not contain un-aligned reads).
    :param mapped_reads: path to a SAM/BAM file containing the alignments
    :param hash_reads: a hash table of read_names to (x,y,umi) tags
    :param min_length: the min number of mapped bases we enforce in an alignment
    :param file_output: the path where to put the records
    :param file_output_discarded: the path where to put discarded files
    :type mapped_reads: str
    :type hash_reads: dict
    :type min_length: integer
    :type file_output: str
    :type file_output_discarded: str
    :raises: RuntimeError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(mapped_reads):
        error = "Error, input file not present {}\n".format(mapped_reads)
        logger.error(error)
        raise RuntimeError(error)
    
    # Create output files handlers
    flag_read = "rb"
    flag_write = "wb"
    infile = pysam.AlignmentFile(mapped_reads, flag_read)
    outfile = pysam.AlignmentFile(file_output, flag_write, template=infile)
    if file_output_discarded is not None:
        outfile_discarded = pysam.AlignmentFile(file_output_discarded, 
                                                flag_write, template=infile)
    # Create some counters and loop the records
    dropped_secondary = 0
    dropped_short = 0
    dropped_barcode = 0
    present = 0
    for sam_record in infile.fetch(until_eof=True):
        present += 1
        discard_read = False
        
        # Add the barcode and coordinates info if present otherwise discard
        try:
            # The probability of a collision is very very low
            key = hash(sam_record.query_name)
            for tag in hash_reads[key]:
                tag_tokens = tag.split(":")
                sam_record.set_tag(tag_tokens[0], tag_tokens[2], tag_tokens[1])
        except KeyError:
            dropped_barcode += 1
            continue
            
        # Get how many bases were mapped
        mapped_bases = 0
        for cigar_tuple in sam_record.cigartuples:
            if cigar_tuple[0] == 0:
                mapped_bases += cigar_tuple[1]
        
        # We need this so we don't duplicate reads
        if not sam_record.is_secondary:
            sam_record.set_tag("NH", None)
            
        # Discard if secondary alignment or only few bases mapped  
        if sam_record.is_secondary:
            dropped_secondary += 1
            discard_read = True
        elif mapped_bases != 0 and mapped_bases < min_length:
            dropped_short += 1
            discard_read = True

        if discard_read:
            if file_output_discarded is not None:
                outfile_discarded.write(sam_record)
        else:
            outfile.write(sam_record)
    
    # Close handlers           
    infile.close()
    outfile.close()
    if file_output_discarded is not None:
        outfile_discarded.close()

    if not fileOk(file_output):
        error = "Error filtering mapped reads.\n" \
        "Output file is not present\n {}".format(file_output)
        logger.error(error)
        raise RuntimeError(error)
            
    logger.info("Finish filtering mapped reads, stats:" \
                "\nPresent: {0}" \
                "\nDropped - secondary alignment: {1}" \
                "\nDropped - too short: {2}" \
                "\nDropped - barcode: {3}".format(present,
                                                  dropped_secondary,
                                                  dropped_short,
                                                  dropped_barcode))
    
    # Update QA object 
    qa_stats.reads_after_mapping = present - (dropped_secondary + dropped_short)