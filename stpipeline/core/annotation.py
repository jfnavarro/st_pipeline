"""
This module contains a modified version of htseq-count
with modifications and improvements to perform annotation
of ST mapped reads (BAM file).
"""
import logging
import os
from typing import Optional, List
import pysam
from stpipeline.common.utils import file_ok
import HTSeq  # type: ignore


class UnknownChrom(Exception):
    """
    Exception raised for unknown chromosome errors.
    """

    pass


def invert_strand(iv: HTSeq.GenomicInterval) -> HTSeq.GenomicInterval:
    """
    Inverts the strand of a genomic interval.

    Args:
        iv: A genomic interval.

    Returns:
        A copy of the genomic interval with the strand inverted.

    Raises:
        ValueError: If the strand is not '+' or '-'.
    """
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def count_reads_in_features(
    sam_filename: str,
    gff_filename: str,
    samtype: str,
    stranded: str,
    overlap_mode: str,
    feature_type: List[str],
    id_attribute: str,
    minaqual: int,
    samout: str,
    include_non_annotated: bool,
    htseq_no_ambiguous: bool,
    outputDiscarded: Optional[str],
) -> int:
    """
    Counts reads in features using a modified version of HTSeq.

    Args:
        sam_filename: Path to the input SAM/BAM file.
        gff_filename: Path to the GFF file with features.
        samtype: Input file type ('sam' or 'bam').
        stranded: Strand specificity ('yes', 'no', 'reverse').
        overlap_mode: Overlap mode for counting ('union', 'intersection-strict', 'intersection-nonempty').
        feature_type: List of feature types to include (e.g., ['exon']).
        id_attribute: Attribute to use as feature ID (e.g., 'gene_id').
        minaqual: Minimum mapping quality for reads.
        samout: Path to the output SAM/BAM file.
        include_non_annotated: Whether to include non-annotated reads.
        htseq_no_ambiguous: Whether to exclude ambiguous reads.
        outputDiscarded: Path to the output file for discarded reads.

    Returns:
        The number of annotated reads.

    Raises:
        RuntimeError: If an error occurs during processing.
        ValueError: If an invalid overlap mode or strand information is encountered.
    """
    # Set up filters
    filter_htseq = ["__too_low_aQual", "__not_aligned"]
    if not include_non_annotated:
        filter_htseq.append("__no_feature")
    filter_htseq_no_ambiguous = htseq_no_ambiguous

    # Open SAM/BAM output file
    flag_write = "wb" if samtype == "bam" else "wh"
    flag_read = "rb" if samtype == "bam" else "r"
    saminfile = pysam.AlignmentFile(sam_filename, flag_read)  # type: ignore
    samoutfile = pysam.AlignmentFile(samout, flag_write, template=saminfile)  # type: ignore
    if outputDiscarded is not None:
        samdiscarded = pysam.AlignmentFile(outputDiscarded, flag_write, template=saminfile)  # type: ignore
    saminfile.close()

    # Counter for annotated records
    annotated = 0

    def write_to_samout(read: HTSeq.SAM_Alignment, assignment: str) -> None:
        """
        Writes a read and its assignment to the SAM output file.
        """
        sam_record = read.to_pysam_AlignedSegment(samoutfile)
        sam_record.set_tag("XF", assignment, "Z")
        if (
            read is not None
            and assignment not in filter_htseq
            and not (filter_htseq_no_ambiguous and "__ambiguous" in assignment)
        ):
            samoutfile.write(sam_record)
        elif outputDiscarded is not None:
            samdiscarded.write(sam_record)

    # Load features
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    counts = {}
    gff = HTSeq.GFF_Reader(gff_filename)

    if feature_type is None:
        feature_type = ["exon"]

    for f in gff:
        if f.type in feature_type:
            feature_id = f.attr.get(id_attribute)
            if feature_id is None:
                raise ValueError(f"Feature {f.name} does not contain a {id_attribute} attribute.")
            if stranded != "no" and f.iv.strand == ".":
                raise ValueError(f"Feature {f.name} at {f.iv} does not have strand information.")
            features[f.iv] += feature_id
            counts[feature_id] = 0

    if not counts:
        raise RuntimeError(f"No features of type {','.join(feature_type)} found.")

    SAM_or_BAM_Reader = HTSeq.SAM_Reader if samtype == "sam" else HTSeq.BAM_Reader

    try:
        read_seq = SAM_or_BAM_Reader(sam_filename)
    except Exception as e:
        raise RuntimeError("Error reading SAM/BAM file.") from e

    for r in read_seq:
        if not r.aligned:
            write_to_samout(r, "__not_aligned")
            continue
        if r.aQual < minaqual:
            write_to_samout(r, "__too_low_aQual")
            continue
        if stranded != "reverse":
            iv_seq = (co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0)
        else:
            iv_seq = (invert_strand(co.ref_iv) for co in r.cigar if co.type in ["M", "=", "X"] and co.size > 0)

        fs = set()  # type: ignore
        for iv in iv_seq:
            if iv.chrom not in features.chrom_vectors:
                write_to_samout(r, "__no_feature")
                break
            for _, fs2 in features[iv].steps():
                fs = fs.union(fs2) if overlap_mode == "union" else fs.intersection(fs2) if fs else fs2.copy()

        if not fs:
            write_to_samout(r, "__no_feature")
        elif len(fs) > 1:
            write_to_samout(r, f"__ambiguous[{'+'.join(fs)}]")
            annotated += 1
        else:
            write_to_samout(r, list(fs)[0])
            annotated += 1

    samoutfile.close()
    if outputDiscarded is not None:
        samdiscarded.close()

    return annotated


def annotateReads(
    mappedReads: str,
    gtfFile: str,
    outputFile: str,
    outputDiscarded: Optional[str],
    mode: str,
    strandness: str,
    htseq_no_ambiguous: bool,
    include_non_annotated: bool,
    feature_types: List[str],
) -> int:
    """
    Annotates a BAM file with mapped reads using HTSeq and writes annotated records to a file.

    Args:
        mappedReads: Path to a BAM file with mapped reads sorted by coordinate.
        gtfFile: Path to an annotation file in GTF format.
        outputFile: Path to write the annotated records (BAM).
        outputDiscarded: Path to write the non-annotated records (BAM).
        mode: Overlapping mode ('union', 'intersection-strict', 'intersection-nonempty').
        strandness: Strand specificity ('yes', 'no', 'reverse').
        htseq_no_ambiguous: Whether to discard ambiguous annotations.
        include_non_annotated: Whether to include non-annotated reads as '__no_feature' in the output.
        feature_types: List of feature types to use for annotation (default is ['exon']).

    Returns:
        The number of annotated reads.

    Raises:
        RuntimeError: If input files are missing or errors occur during annotation.
    """
    logger = logging.getLogger("STPipeline")

    if not os.path.isfile(mappedReads):
        raise RuntimeError(f"Input file not found: {mappedReads}")

    if feature_types is None:
        feature_types = ["exon"]
    try:
        annotated = count_reads_in_features(
            mappedReads,
            gtfFile,
            "bam",  # Type BAM for files
            strandness,  # Strand yes/no/reverse
            mode,  # intersection_nonempty, union, intersection_strict
            feature_types,  # feature types in GFF
            "gene_id",  # gene_id or gene_name
            0,  # Min quality score
            outputFile,
            include_non_annotated,
            htseq_no_ambiguous,
            outputDiscarded,
        )
    except Exception as e:
        error = "Error during annotation. HTSEQ execution failed"
        logger.error(error)
        raise e

    if not file_ok(outputFile) or annotated == 0:
        error = f"Error during annotation. Output file not present or empty {outputFile}"
        logger.error(error)
        raise RuntimeError(error)

    logger.info(f"Annotated reads: {annotated}")
    return annotated
