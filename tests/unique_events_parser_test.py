#! /usr/bin/env python
"""
Unit-test the package unique_events_parser
"""

import pysam
import pytest

from stpipeline.common.unique_events_parser import GeneBuffer, Transcript, parse_unique_events


@pytest.fixture
def mock_gff_file(tmp_path):
    gff_content = (
        "chr1\tsource\tgene\t1\t1000\t.\t+\t.\tgene_id=gene1;\n"
        "chr1\tsource\tgene\t1001\t2000\t.\t+\t.\tgene_id=gene2;\n"
        "chr2\tsource\tgene\t1\t1500\t.\t-\t.\tgene_id=gene3;\n"
    )
    gff_file = tmp_path / "mock.gff"
    with open(gff_file, "w") as f:
        f.write(gff_content)
    return str(gff_file)


@pytest.fixture
def mock_bam_file(tmp_path):
    bam_file = tmp_path / "mock.bam"
    with pysam.AlignmentFile(bam_file, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 2000, "SN": "chr1"}]}) as f:
        for i in range(5):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{i}"
            segment.query_sequence = "ACTG" * 25
            segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = i * 100
            segment.cigar = [(0, len(segment.query_sequence))]  # 0: MATCH
            segment.set_tag("B1", i)
            segment.set_tag("B2", i * 2)
            segment.set_tag("XF", "gene1")
            segment.set_tag("B3", "UMI1")
            f.write(segment)
    return str(bam_file)


def test_compute_gene_end_coordinates(mock_gff_file):
    buffer = GeneBuffer(mock_gff_file)
    assert buffer.gene_end_coordinates["gene1"] == ("chr1", 1000)
    assert buffer.gene_end_coordinates["gene2"] == ("chr1", 2000)
    assert buffer.gene_end_coordinates["gene3"] == ("chr2", 1500)
    assert buffer.gene_end_coordinates["__no_feature"] == (None, -1)


def test_add_transcript(mock_gff_file):
    buffer = GeneBuffer(mock_gff_file)
    transcript = Transcript("chr1", 100, 200, "read1", 60, "+", "UMI1")
    buffer.add_transcript("gene1", (1, 2), transcript, 100)

    assert "gene1" in buffer.buffer
    assert (1, 2) in buffer.buffer["gene1"]
    assert buffer.buffer["gene1"][(1, 2)][0] == transcript


def test_check_and_clear_buffer(mock_gff_file):
    buffer = GeneBuffer(mock_gff_file)
    transcript = Transcript("chr1", 100, 200, "read1", 60, "+", "UMI1")
    buffer.add_transcript("gene1", (1, 2), transcript, 100)
    buffer.last_position = 300

    cleared_genes = list(buffer.check_and_clear_buffer(empty=True))
    assert len(cleared_genes) == 1
    assert cleared_genes[0][0] == "gene1"
    assert buffer.buffer == {}


def test_check_and_clear_buffer_no_feature(mock_gff_file):
    buffer = GeneBuffer(mock_gff_file)
    transcript = Transcript("chr1", 100, 200, "read1", 60, "+", "UMI1")
    buffer.add_transcript("__no_feature", (1, 2), transcript, 100)
    buffer.last_position = 300

    cleared_genes = list(buffer.check_and_clear_buffer(empty=True))
    assert len(cleared_genes) == 1
    assert cleared_genes[0][0] == "__no_feature"
    assert buffer.buffer == {}


def test_get_gene_end_position_ambiguous(mock_gff_file):
    buffer = GeneBuffer(mock_gff_file)
    ambiguous_gene = "__ambiguous[gene1+gene2]"
    chrom, end_position = buffer.get_gene_end_position(ambiguous_gene)
    assert chrom == "chr1"
    assert end_position == 2000


def test_parse_unique_events(mock_bam_file, mock_gff_file):
    unique_events = list(parse_unique_events(mock_bam_file, mock_gff_file))

    assert len(unique_events) == 1
    gene, spots = unique_events[0]
    assert gene == "gene1"
    assert len(spots) == 5
    assert spots[(0, 0)][0].clear_name == "read0"


def test_parse_unique_events_no_annotation(mock_bam_file):
    unique_events = list(parse_unique_events(mock_bam_file))

    assert len(unique_events) == 1
    gene, spots = unique_events[0]
    assert gene == "gene1"
    assert len(spots) == 5
    assert spots[(0, 0)][0].clear_name == "read0"
