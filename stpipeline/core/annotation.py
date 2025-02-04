"""
This module contains a modified version of htseq-count
with modifications and improvements to perform annotation
of ST mapped reads (BAM file).
"""

import logging
import os
from typing import List, Optional, Tuple, Dict

import HTSeq  # type: ignore
import pysam

from stpipeline.common.utils import file_ok

logger = logging.getLogger("STPipeline")


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
        raise ValueError(f"Illegal strand {iv}")
    return iv2


class ReadCounter:
    def __init__(
        self,
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
        output_discarded: Optional[str] = None,
    ) -> None:
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
        """
        self.sam_filename = sam_filename
        self.gff_filename = gff_filename
        self.samtype = samtype
        self.stranded = stranded
        self.overlap_mode = overlap_mode
        self.feature_type = feature_type
        self.id_attribute = id_attribute
        self.minaqual = minaqual
        self.samout = samout
        self.include_non_annotated = include_non_annotated
        self.htseq_no_ambiguous = htseq_no_ambiguous
        self.output_discarded = output_discarded
        self.annotated = 0
        self._validate_inputs()
        self.features, self.counts = self._load_features()
        self.samoutfile, self.samdiscarded = self._open_sam_files()
        self.filter_htseq = ["__too_low_aQual", "__not_aligned"]
        if not self.include_non_annotated:
            self.filter_htseq.append("__no_feature")

    def __del__(self) -> None:
        """Ensures file handles are closed when the object is destroyed."""
        if self.samoutfile:
            self.samoutfile.close()
        if self.samdiscarded:
            self.samdiscarded.close()

    def _validate_inputs(self) -> None:
        """Validates input parameters."""
        if self.samtype not in ["bam", "sam"]:
            raise ValueError(f"Incorrect value for samtype {self.samtype}")
        if self.stranded not in ["yes", "no", "reverse"]:
            raise ValueError(f"Incorrect value for stranded option {self.stranded}")
        if self.overlap_mode not in ["union", "intersection-strict", "intersection-nonempty"]:
            raise ValueError(f"Incorrect value for overlap_mode option {self.overlap_mode}")
        if not self.feature_type:
            raise ValueError("Value types cannot be empty")

    def _load_features(self) -> Tuple[HTSeq.GenomicArrayOfSets, Dict[str, int]]:
        """Loads genomic features from the GFF file."""
        features = HTSeq.GenomicArrayOfSets("auto", self.stranded != "no")
        counts = {}
        gff = HTSeq.GFF_Reader(self.gff_filename)
        for f in gff:
            if f.type in self.feature_type:
                feature_id = f.attr.get(self.id_attribute)
                if feature_id is None:
                    raise ValueError(f"Feature {f.name} does not contain a {self.id_attribute} attribute.")
                if self.stranded != "no" and f.iv.strand == ".":
                    raise ValueError(f"Feature {f.name} at {f.iv} does not have strand information.")
                features[f.iv] += feature_id
                counts[feature_id] = 0

        if not counts:
            raise RuntimeError(f"No features of type {','.join(self.feature_type)} found.")

        return features, counts

    def _open_sam_files(self) -> Tuple[pysam.AlignmentFile, Optional[pysam.AlignmentFile]]:
        """Opens SAM/BAM files for output and discarded reads."""
        flag_write = "wb" if self.samtype == "bam" else "wh"
        flag_read = "rb" if self.samtype == "bam" else "r"
        saminfile = pysam.AlignmentFile(self.sam_filename, flag_read)  # type: ignore
        samoutfile = pysam.AlignmentFile(self.samout, flag_write, template=saminfile)  # type: ignore
        samdiscarded = None
        if self.output_discarded is not None:
            samdiscarded = pysam.AlignmentFile(self.output_discarded, flag_write, template=saminfile)  # type: ignore
        saminfile.close()
        return samoutfile, samdiscarded

    def _write_to_samout(self, read: HTSeq.SAM_Alignment, assignment: str) -> None:
        """Writes a read and its assignment to the SAM output file."""
        sam_record = read.to_pysam_AlignedSegment(self.samoutfile)
        sam_record.set_tag("XF", assignment, "Z")
        if (
            read is not None
            and assignment not in self.filter_htseq
            and not (self.htseq_no_ambiguous and "__ambiguous" in assignment)
        ):
            self.samoutfile.write(sam_record)
            self.annotated += 1
        elif self.output_discarded is not None:
            self.samdiscarded.write(sam_record)  # type: ignore

    def count_reads(self) -> int:
        """Counts reads in features using a modified version of HTSeq."""
        self.annotated = 0
        try:
            read_seq = (
                HTSeq.SAM_Reader(self.sam_filename) if self.samtype == "sam" else HTSeq.BAM_Reader(self.sam_filename)
            )
            com = ("M", "=", "X")
            for r in read_seq:
                if not r.aligned:
                    self._write_to_samout(r, "__not_aligned")
                    continue
                if r.aQual < self.minaqual:
                    self._write_to_samout(r, "__too_low_aQual")
                    continue
                iv_seq = (
                    (co.ref_iv for co in r.cigar if co.type in com and co.size > 0)
                    if self.stranded != "reverse"
                    else (invert_strand(co.ref_iv) for co in r.cigar if co.type in com and co.size > 0)
                )
                fs = set()  # type: ignore
                for iv in iv_seq:
                    if iv.chrom not in self.features.chrom_vectors:
                        fs = set()
                        break
                    for _, fs2 in self.features[iv].steps():
                        if self.overlap_mode == "union":
                            fs = fs.union(fs2)
                        elif len(fs2) > 0 or self.overlap_mode == "intersection-strict":
                            fs = fs.intersection(fs2) if fs else fs2.copy()
                if not fs:
                    self._write_to_samout(r, "__no_feature")
                elif len(fs) > 1:
                    self._write_to_samout(r, f"__ambiguous[{'+'.join(fs)}]")
                else:
                    self._write_to_samout(r, list(fs)[0])
        except Exception as e:
            raise RuntimeError("Error encountered during read counting") from e
        return self.annotated


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

    if not os.path.isfile(mappedReads):
        raise RuntimeError(f"Input file not found: {mappedReads}")

    try:
        annotated = ReadCounter(
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
        ).count_reads()
    except Exception as e:
        error = "Error during annotation: HTSEQ execution failed"
        logger.error(error)
        raise e

    if not file_ok(outputFile) or annotated == 0:
        error = f"Error during annotation: Output file not present or empty {outputFile}"
        logger.error(error)
        raise RuntimeError(error)

    logger.info(f"Annotated reads: {annotated}")
    return annotated
