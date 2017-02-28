""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import os
import logging 
import pysam
from collections import defaultdict

# TODO this function uses too much memory, optimize it. (Maybe Cython or C++)
def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a SAM/BAM file where the spot coordinates, gene and UMI are present as extra tags
    The output will be a dictionary [spot][gene] -> (chrom, start, end, clear_name, mapping_quality, strand, umi). 
    :param filename: the input file containing the annotated SAM/BAM records
    :return: A dictionary of spots(x,y) to a map of gene names to a list of transcripts 
    (chrom, start, end, clear_name, mapping_quality, strand, umi)
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
        # Account for soft-clipped bases when retrieving the start/end coordinates
        start = rec.reference_start - rec.query_alignment_start
        end = rec.reference_end + (rec.query_length - rec.query_alignment_end)
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

def filterMappedReads(mapped_reads,
                      hash_reads,
                      file_output,
                      file_output_discarded=None,
                      min_length=30):
    """ 
    Iterates a SAM/BAM file containing mapped reads 
    and discard the reads that are secondary or too short.
    It also discards reads that are not demultiplexed with TaggD
    (for that a dictionary with the read name as key and the X,Y and UMI)
    as values must be given.
    This function will add the X,Y coordinates and UMI as extra tags
    to the output SAM/BAM file. 
    It assumes all the reads are aligned (do not contain un-aligned reads).
    :param mapped_reads: path to a SAM/BAM file containing the alignments
    :param hash_reads: a dictionary of read_names to (x,y,umi) SAM tags
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
        # Add the UMI and X,Y coordinates as extra SAM tags
        try:
            # Using as key the read name as it was used to generate the dictionary
            # In order to save memory we truncate the read
            # name to only keep the unique part (lane, tile, x_pos, y_pos)
            # TODO this procedure is specific to only Illumina technology
            key = "".join(sam_record.query_name.split(":")[-4:])
            for tag in hash_reads[key]:
                # TODO add error check here
                tag_tokens = tag.split(":")
                sam_record.set_tag(tag_tokens[0], tag_tokens[2], tag_tokens[1])
        except KeyError:
            dropped_barcode += 1
            if file_output_discarded is not None:
                outfile_discarded.write(sam_record)
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