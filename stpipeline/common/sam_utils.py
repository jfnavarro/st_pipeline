""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

import pysam
        
def convert_to_AlignedSegment(header, sequence, quality, 
                              barcode_sequence, umi_sequence):
    """
    This function converts the input variables 
    (header,sequence,quality,barcode_sequence,umi_sequence)
    to a unaligned pysam.AlignedSegment with the umi and barcode 
    informations as the following tags:
        Tag  Value
        "B0" barcode_sequence
        "B3" umi_sequence
    :param header: string with the header information
    :param sequence: string with the DNA/RNA sequence
    :param quality: string with the base calling quality values
    :param barcode_sequence: string with the barcode sequence
    :param umi_sequence: string with the unique molecular identifier sequence
    """

    # create
    aligned_segment = pysam.AlignedSegment()

    # Set the standard values
    # Header must not contain empty spaces
    aligned_segment.query_name = header.split()[0]
    aligned_segment.query_sequence = sequence
    aligned_segment.query_qualities = pysam.qualitystring_to_array(quality)

    # setting the flag to un_mapped
    aligned_segment.flag |= pysam.FUNMAP

    # Set the tags
    aligned_segment.set_tag('B0', barcode_sequence)
    aligned_segment.set_tag('B3', umi_sequence)
    aligned_segment.set_tag('RG', '0')

    return aligned_segment

def merge_bam(merged_file_name, files_to_merge, ubam=False):
    """
    Function for merging partial bam files after annotation.
    also counts the number of reads for different types of annotations (the XF tags of the reads)
    and returns these counts as a dict: anotation=>count
    :param merged_file_name: name of the merged output bamfile
    :param files_to_merge: list with names of the partial bamfiles to merge
    :param ubam: indicates unaligned bam file (True or False, default False)
    :returns: the number of annotated records
    """
    annotations = {}
    with pysam.AlignmentFile(files_to_merge[0], mode='rb', check_sq=(not ubam)) as input_bamfile:
        merged_file = pysam.AlignmentFile(merged_file_name,
                                          mode="wb", template=input_bamfile)
    # Simply merges the BAM files and creates a counter of annotated records
    for file_name in files_to_merge:
        input_bamfile = pysam.AlignmentFile( file_name, mode='rb', check_sq=(not ubam) )
        for record in input_bamfile.fetch(until_eof=True):
            merged_file.write(record)
            if ubam: annotation = None
            else: annotation = record.get_tag("XF")
            try:
                annotations[annotation] += 1
            except KeyError:
                annotations[annotation] = 1
    return sum(annotations.values())