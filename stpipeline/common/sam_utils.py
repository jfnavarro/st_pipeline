""" 
This module contains some functions and utilities for SAM/BAM files
"""

import pysam
import os
import math


def split_bam(input_bamfile_name, temp_dir, threads):
    """
    Splits a BAM file in to parts with equal read counts.
    The number of parts to split the BAM File into equals the 
    number of cores given as input.
    :param input_bamfile_name: path to the BAM file to be splitted
    :param temp_dir: path to the folder where to put the created files
    :param threads: the number of CPU cores to use
    :retuns: the list of splitted BAM files
    """
    # Index and open the input BAM
    pysam.index(input_bamfile_name,
                os.path.join(temp_dir, '{0}.bai'.format(input_bamfile_name)))
    input_bamfile = pysam.AlignmentFile(input_bamfile_name, mode='rb')
    assert input_bamfile.check_index()

    output_file_names = {part: os.path.join(temp_dir,
                                            "{0}.part_{1}.bam".format(input_bamfile_name, part))
                         for part in range(threads)}
    # Open the output bam files
    output_bamfiles = {
        part: pysam.AlignmentFile(file_name, mode="wbu", template=input_bamfile) \
        for part, file_name in output_file_names.iteritems()
    }

    # Split the BAM file
    total_read_count = input_bamfile.mapped + input_bamfile.unmapped
    reads_per_part = math.ceil(float(total_read_count) / threads)
    _tmp_read_counter = 0
    part = 0
    for record in input_bamfile.fetch(until_eof=True):
        output_bamfiles[part].write(record)
        _tmp_read_counter += 1
        if _tmp_read_counter == reads_per_part:
            part += 1
            _tmp_read_counter = 0
    input_bamfile.close()
    return output_file_names.values()


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
    Function for merging partial BAM files into one.
    :param merged_file_name: name of the merged output bam file
    :param files_to_merge: list with names of the partial bam files to merge
    :param ubam: indicates unaligned bam file (True or False, default False)
    :returns: the total number of records
    """
    assert files_to_merge is not None and len(files_to_merge) > 0
    num_ele = 0
    with pysam.AlignmentFile(files_to_merge[0], mode='rb',
                             check_sq=(not ubam)) as input_bamfile:
        merged_file = pysam.AlignmentFile(merged_file_name,
                                          mode="wb", template=input_bamfile)
        # Simply merges the BAM files and creates a counter of annotated records
        for file_name in files_to_merge:
            input_bamfile = pysam.AlignmentFile(file_name, mode='rb', check_sq=(not ubam))
            for record in input_bamfile.fetch(until_eof=True):
                merged_file.write(record)
                num_ele += 1
            input_bamfile.close()
        merged_file.close()
    return num_ele
