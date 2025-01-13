"""
This module contains some functions and utilities for SAM/BAM files
"""

import math
import os
from typing import List

import pysam


def split_bam(input_bam: str, temp_dir: str, threads: int) -> List[str]:
    """
    Splits a BAM file into chunks with equal read counts. The number of chunks
    equals the number of CPU cores specified.

    Args:
        input_bam: Path to the BAM file to be split.
        temp_dir: Directory where the created files will be stored.
        threads: Number of CPU cores to use for splitting.

    Returns:
        List of paths to the split BAM files.
    """
    pysam.index(input_bam, os.path.join(temp_dir, f"{input_bam}.bai"))  # type: ignore
    input_bamfile = pysam.AlignmentFile(input_bam, mode="rb")
    assert input_bamfile.check_index()

    output_file_names = {part: os.path.join(temp_dir, f"{input_bam}.part_{part}.bam") for part in range(threads)}

    output_bamfiles = {
        part: pysam.AlignmentFile(file_name, mode="wb", template=input_bamfile)
        for part, file_name in output_file_names.items()
    }

    total_read_count = input_bamfile.mapped + input_bamfile.unmapped
    reads_per_part = math.ceil(total_read_count / threads)
    read_counter = 0
    part = 0
    for record in input_bamfile.fetch(until_eof=True):
        output_bamfiles[part].write(record)
        read_counter += 1
        if read_counter == reads_per_part:
            part += 1
            read_counter = 0

    input_bamfile.close()
    return list(output_file_names.values())


def convert_to_AlignedSegment(
    header: str, sequence: str, quality: str, barcode_sequence: str, umi_sequence: str
) -> pysam.AlignedSegment:
    """
    Converts input variables to an unaligned `pysam.AlignedSegment` with UMI and
    barcode information as tags.

    Args:
        header: Header information for the segment.
        sequence: DNA/RNA sequence.
        quality: Base calling quality values.
        barcode_sequence: Barcode sequence.
        umi_sequence: Unique molecular identifier sequence.

    Returns:
        A new AlignedSegment object with the provided data.
    """
    aligned_segment = pysam.AlignedSegment()
    aligned_segment.query_name = header.split()[0]
    aligned_segment.query_sequence = sequence
    aligned_segment.query_qualities = pysam.qualitystring_to_array(quality)
    aligned_segment.flag |= pysam.FUNMAP
    aligned_segment.set_tag("B0", barcode_sequence)
    aligned_segment.set_tag("B3", umi_sequence)
    aligned_segment.set_tag("RG", "0")
    return aligned_segment


def merge_bam(merged_file_name: str, files_to_merge: List[str], ubam: bool = False) -> int:
    """
    Merges multiple partial BAM files into a single file.

    Args:
        merged_file_name: Path to the output merged BAM file.
        files_to_merge: List of paths to the partial BAM files.
        ubam: Indicates if the files are unaligned BAM (uBAM). Default is False.

    Returns:
        Total number of records in the merged BAM file.
    """
    assert files_to_merge, "The list of files to merge cannot be empty."
    num_records = 0

    with pysam.AlignmentFile(files_to_merge[0], mode="rb", check_sq=not ubam) as input_bamfile:
        merged_file = pysam.AlignmentFile(merged_file_name, mode="wb", template=input_bamfile)
        for file_name in files_to_merge:
            with pysam.AlignmentFile(file_name, mode="rb", check_sq=not ubam) as input_file:
                for record in input_file.fetch(until_eof=True):
                    merged_file.write(record)
                    num_records += 1
        merged_file.close()

    return num_records
