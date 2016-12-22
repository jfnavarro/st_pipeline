""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import os
import logging 
import pysam
from collections import defaultdict

#TODO this function uses too much memory, optimize it. (Maybe Cython or C++)
def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM/BAM file where the spot coordinates are present in the tags
    The output will be hash table [spot][gene] -> list of transcripts. 
    :param filename: the input file containing the annotated SAM/BAM records
    :return: A map of spots(x,y) to a map of gene names to a list of transcript 
    (chrom, start, end, clear_name, mapping_quality, strand, sequence)
    As map[(x,y)][gene]->list((chrom, start, end, clear_name, mapping_quality, strand, UMI))
    """
    
    logger = logging.getLogger("STPipeline")
    
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
        # Get TAGGD tags
        x,y,gene,seq = (None,None,None,None)
        for (k, v) in rec.tags:
            if k == "B1":
                x = int(v) ## The X coordinate
            elif k == "B2":
                y = int(v) ## The Y coordinate
            elif k == "XF":
                gene = str(v) ## The gene name
            elif k == "B3":
                umi = str(v) ## The UMI
            else:
                continue
        # Check that all tags are present
        if None in [x,y,gene,umi]:
            logger.warning("Warning parsing annotated reads.\n" \
                           "Missing attributes for record {}\n".format(clear_name))
            continue
        
        # Create a new transcript and add it to the dictionary
        transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
        unique_events[(x,y)][gene].append(transcript)
      
    sam_file.close()
    return unique_events

def sortSamFile(input_sam, outputFolder=None):
    """
    It simply sorts by position a sam/bam file containing mapped reads 
    :param input: is a SAM/BAM file with mapped reads
    :param outputFolder: the location where to place the output file (optional)
    :type input: str
    :type outputFolder: str
    :return: the path to the sorted file
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
    and discards the reads that are secondary or too short.
    It also discards reads that do not contain a valid barcode.
    It will add the barcode, coordinates and UMI as extra tags
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
    sam_type = os.path.splitext(mapped_reads)[1].lower()
    flag_read = "r" if sam_type == ".sam" else "rb"
    flag_write = "w" if sam_type == ".sam" else "wb"
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
            # Using as key the read name as it was used to generate the dictionary
            for tag in hash_reads[sam_record.query_name]:
                # TODO add error check here
                tag_tokens = tag.split(":")
                sam_record.set_tag(tag_tokens[0], tag_tokens[2], tag_tokens[1])
        except KeyError:
            dropped_barcode += 1
            continue
            
        # Get how many bases were mapped
        mapped_bases = sum([cigar[1] for cigar in sam_record.cigartuples if cigar[0] == 0])
            
        # Discard if secondary alignment or only few bases mapped  
        if sam_record.is_secondary:
            dropped_secondary += 1
            discard_read = True
        elif mapped_bases != 0 and mapped_bases < min_length:
            dropped_short += 1
            discard_read = True
        else:
            # We need this so htseq-count
            # does not discard the read thinking that it is secondary
            sam_record.set_tag("NH", 1)
                               
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
    qa_stats.reads_after_mapping = present - \
    (dropped_secondary + dropped_short + dropped_barcode)