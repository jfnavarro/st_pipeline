"""
A module that contains functions to parse and filter input reads for ST data processing.
"""
import os
import pysam
import logging
import dnaio
from stpipeline.common.fastq_utils import trim_quality, remove_adaptor, check_umi_template, has_sufficient_content
from stpipeline.common.sam_utils import convert_to_AlignedSegment
from typing import Optional, Tuple

logger = logging.getLogger("STPipeline")

bam_header = {"HD": {"VN": "1.5", "SO": "unsorted"}, "RG": [{"ID": "0", "SM": "unknown_sample", "PL": "ILLUMINA"}]}


def filter_input_data(
    fw_file: str,
    rv_file: str,
    out_file: str,
    out_file_discarded: Optional[str],
    barcode_length: int,
    start_position: int,
    filter_AT_content: float,
    filter_GC_content: float,
    umi_start: int,
    umi_end: int,
    min_qual: int,
    min_length: int,
    polyA_min_distance: int,
    polyT_min_distance: int,
    polyG_min_distance: int,
    polyC_min_distance: int,
    polyN_min_distance: int,
    qual64: bool,
    umi_filter: bool,
    umi_filter_template: str,
    umi_quality_bases: int,
    adaptor_missmatches: int,
    overhang: int,
    disable_umi: bool,
    disable_barcode: bool,
) -> Tuple[int, int]:
    """
    Handles input read filtering and quality trimming for sequencing data (paired FASTQ files).
    - It performs a sanity check (forward and reverse reads same length and order)
    - It performs a BWA-based quality trimming discarding very short reads
    - It removes adaptors from the reads (optional)
    - It checks for AT and GC content (optional)
    - It performs a sanity check on the UMI (optional)
    Reads that do not pass the filters are discarded (both R1 and R2)
    Reads that pass the filter are written as BAM (R2)

    Args:
        fw_file: Path to the FASTQ file containing forward reads (R1).
        rv_file: Path to the FASTQ file containing reverse reads (R2).
        out_file: Path to the output BAM file.
        out_file_discarded: Path to the output FASTQ file for discarded reads.
        barcode_length: Length of the barcode sequence.
        start_position: Starting position of the barcode in the sequence.
        filter_AT_content: Maximum allowed percentage of A and T bases in a read for filtering.
        filter_GC_content: Maximum allowed percentage of G and C bases in a read for filtering.
        umi_start: Starting position of the UMI in the sequence.
        umi_end: Ending position of the UMI in the sequence.
        min_qual: Minimum quality threshold for quality trimming.
        min_length: Minimum valid length for a read after trimming.
        polyA_min_distance: Minimum distance for PolyA adaptor trimming.
        polyT_min_distance: Minimum distance for PolyT adaptor trimming.
        polyG_min_distance: Minimum distance for PolyG adaptor trimming.
        polyC_min_distance: Minimum distance for PolyC adaptor trimming.
        polyN_min_distance: Minimum distance for PolyN adaptor trimming.
        qual64: True if quality scores are in Phred64 format, False for Phred33.
        umi_filter: If True, applies UMI quality template filtering.
        umi_filter_template: Template for UMI quality filtering.
        umi_quality_bases: Maximum number of low-quality bases allowed in a UMI.
        adaptor_missmatches: Number of mismatches allowed when removing adaptors.
        overhang: Overhang for barcode extraction.
        disable_umi: If True, skips UMI filtering.
        disable_barcode: If True, skips barcode extraction.

    Returns:
        Total reads processed and remaining reads after filtering.

    Raises:
        RuntimeError: If input files are missing or errors occur during processing.
    """
    if not os.path.isfile(fw_file) or not os.path.isfile(rv_file):
        error = f"Error doing quality trimming, input file/s not present {fw_file} {rv_file}"
        logger.error(error)
        raise RuntimeError(error)

    keep_discarded_files = out_file_discarded is not None

    # Build fake sequence adaptors with the parameters given
    adaptorA = "".join("A" for _ in range(polyA_min_distance))
    adaptorT = "".join("T" for _ in range(polyT_min_distance))
    adaptorG = "".join("G" for _ in range(polyG_min_distance))
    adaptorC = "".join("C" for _ in range(polyC_min_distance))
    adaptorN = "".join("N" for _ in range(polyN_min_distance))

    # Quality format
    phred = 64 if qual64 else 33

    # Some counters
    total_reads = 0
    dropped_umi = 0
    dropped_umi_template = 0
    dropped_AT = 0
    dropped_GC = 0
    dropped_adaptor = 0
    too_short_after_trimming = 0

    bam_file = pysam.AlignmentFile(out_file, "wb", header=bam_header)
    if keep_discarded_files:
        out_writer_discarded = dnaio.open(out_file_discarded, mode="w")  # type: ignore

    with dnaio.open(fw_file, rv_file) as reader:
        for r1, r2 in reader:
            header_fw, sequence_fw, quality_fw = r1.name, r1.sequence, r1.qualities
            header_rv, sequence_rv, quality_rv = r2.name, r2.sequence, r2.qualities
            orig_sequence_rv, orig_quality_rv = sequence_rv, quality_rv
            discard_read = False
            total_reads += 1

            if header_fw.split()[0] != header_rv.split()[0]:
                logger.warning(f"Pair reads found with different names {header_fw} and {header_rv}.")

            if not disable_barcode:
                barcode = sequence_fw[max(0, start_position - overhang) : (start_position + barcode_length + overhang)]
            else:
                barcode = None

            if not disable_umi:
                umi_seq = sequence_fw[umi_start:umi_end]
                if umi_filter and not check_umi_template(umi_seq, umi_filter_template):
                    dropped_umi_template += 1
                    discard_read = True

                umi_qual = quality_fw[umi_start:umi_end]
                if not discard_read and len([b for b in umi_qual if (ord(b) - phred) < min_qual]) > umi_quality_bases:
                    dropped_umi += 1
                    discard_read = True
            else:
                umi_seq = None

            if (
                not discard_read
                and filter_AT_content > 0
                and has_sufficient_content(sequence_rv, "AT", filter_AT_content)
            ):
                dropped_AT += 1
                discard_read = True

            if (
                not discard_read
                and filter_GC_content > 0
                and has_sufficient_content(sequence_rv, "GC", filter_GC_content)
            ):
                dropped_GC += 1
                discard_read = True

            if not discard_read:
                if polyA_min_distance >= 5:
                    sequence_rv, quality_rv = remove_adaptor(sequence_rv, quality_rv, adaptorA, adaptor_missmatches)
                if polyT_min_distance >= 5:
                    sequence_rv, quality_rv = remove_adaptor(sequence_rv, quality_rv, adaptorT, adaptor_missmatches)
                if polyG_min_distance >= 5:
                    sequence_rv, quality_rv = remove_adaptor(sequence_rv, quality_rv, adaptorG, adaptor_missmatches)
                if polyC_min_distance >= 5:
                    sequence_rv, quality_rv = remove_adaptor(sequence_rv, quality_rv, adaptorC, adaptor_missmatches)
                if polyN_min_distance >= 5:
                    sequence_rv, quality_rv = remove_adaptor(sequence_rv, quality_rv, adaptorN, adaptor_missmatches)

                if len(sequence_rv) < min_length:
                    dropped_adaptor += 1
                    discard_read = True

            if not discard_read:
                sequence_rv, quality_rv = trim_quality(sequence_rv, quality_rv, min_qual, min_length, phred)
                if not sequence_rv or not quality_rv:
                    too_short_after_trimming += 1
                    discard_read = True

            if not discard_read:
                bam_file.write(convert_to_AlignedSegment(header_rv, sequence_rv, quality_rv, barcode, umi_seq))
            elif keep_discarded_files:
                out_writer_discarded.write(dnaio.SequenceRecord(header_rv, orig_sequence_rv, orig_quality_rv))

    bam_file.close()
    if keep_discarded_files:
        out_writer_discarded.close()

    dropped_rv = (
        dropped_umi + dropped_umi_template + dropped_AT + dropped_GC + dropped_adaptor + too_short_after_trimming
    )
    logger.info(f"Trimming stats total reads (pair): {total_reads}")
    logger.info(f"Trimming stats {dropped_rv} reads have been dropped!")
    logger.info(f"Trimming stats you just lost about {(dropped_rv / total_reads):.2%} of your data")
    logger.info(f"Trimming stats reads remaining: {total_reads - dropped_rv}")

    return total_reads, total_reads - dropped_rv
