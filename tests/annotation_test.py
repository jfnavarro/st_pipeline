#! /usr/bin/env python
"""
Unit-test the package annotation
"""

import os

import HTSeq
import pysam
import pytest

from stpipeline.core.annotation import ReadCounter, annotateReads, invert_strand


@pytest.fixture
def mock_gff_file(tmp_path):
    gff_content = (
        "chr1\tsource\texon\t1\t1000\t.\t+\t.\tgene_id=gene1;\n"
        "chr1\tsource\texon\t1001\t2000\t.\t+\t.\tgene_id=gene2;\n"
        "chr2\tsource\texon\t1\t1500\t.\t-\t.\tgene_id=gene3;\n"
    )
    gff_file = tmp_path / "mock.gff"
    with open(gff_file, "w") as f:
        f.write(gff_content)
    return str(gff_file)


@pytest.fixture
def mock_bam_file(tmp_path):
    bam_file = tmp_path / "mock.bam"
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": 2000, "SN": "chr1"},
            {"LN": 1500, "SN": "chr2"},
        ],
    }
    with pysam.AlignmentFile(bam_file, "wb", header=header) as f:
        for i in range(5):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{i}"
            segment.query_sequence = "ACTG" * 25
            segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = i * 100
            segment.set_tag("B1", i)
            segment.set_tag("B2", i * 2)
            segment.set_tag("XF", "gene1")
            segment.cigar = [(0, len(segment.query_sequence))]  # 0: MATCH
            f.write(segment)
    return str(bam_file)


def test_invert_strand():
    iv = HTSeq.GenomicInterval("chr1", 0, 1000, "+")
    inverted = invert_strand(iv)
    assert inverted.strand == "-"

    iv = HTSeq.GenomicInterval("chr1", 0, 1000, "-")
    inverted = invert_strand(iv)
    assert inverted.strand == "+"

    with pytest.raises(ValueError):
        iv = HTSeq.GenomicInterval("chr1", 0, 1000, ".")
        invert_strand(iv)


def test_count_reads_in_features(mock_bam_file, mock_gff_file, tmp_path):
    output_file = tmp_path / "output.bam"
    discarded_file = tmp_path / "discarded.bam"

    annotated_count = ReadCounter(
        sam_filename=mock_bam_file,
        gff_filename=mock_gff_file,
        samtype="bam",
        stranded="yes",
        overlap_mode="union",
        feature_type=["exon"],
        id_attribute="gene_id",
        minaqual=0,
        samout=str(output_file),
        include_non_annotated=True,
        htseq_no_ambiguous=False,
        output_discarded=str(discarded_file),
    ).count_reads()

    assert annotated_count > 0
    assert os.path.exists(output_file)
    assert os.path.exists(discarded_file)


def test_annotate_reads(mock_bam_file, mock_gff_file, tmp_path):
    output_file = tmp_path / "output.bam"
    discarded_file = tmp_path / "discarded.bam"

    annotateReads(
        mappedReads=mock_bam_file,
        gtfFile=mock_gff_file,
        outputFile=str(output_file),
        outputDiscarded=str(discarded_file),
        mode="union",
        strandness="yes",
        htseq_no_ambiguous=False,
        include_non_annotated=True,
        feature_types=["exon"],
    )

    assert os.path.exists(output_file)
    assert os.path.exists(discarded_file)
