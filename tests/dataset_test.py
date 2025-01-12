#! /usr/bin/env python
""" 
Unit-test the package dataset
"""
import pytest
from typing import List
from stpipeline.common.dataset import Transcript, compute_unique_umis, createDataset
from collections import Counter
from unittest.mock import Mock
import pysam

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

# Test for Transcript Dataclass
def test_transcript_dataclass():
    transcript = Transcript(
        chrom="chr1",
        start=100,
        end=200,
        clear_name="test_transcript",
        mapping_quality=60,
        strand="+",
        umi="ATGC"
    )

    assert transcript.chrom == "chr1"
    assert transcript.start == 100
    assert transcript.end == 200
    assert transcript.clear_name == "test_transcript"
    assert transcript.mapping_quality == 60
    assert transcript.strand == "+"
    assert transcript.umi == "ATGC"

# Test for compute_unique_umis
def mock_group_umi_func(umis: List[str], mismatches: int) -> List[str]:
    return umis[:1]  # Simplified mock implementation for testing

def test_compute_unique_umis():
    transcripts = [
        Transcript("chr1", 100, 200, "t1", 60, "+", "UMI1"),
        Transcript("chr1", 105, 205, "t2", 60, "+", "UMI2"),
        Transcript("chr1", 110, 210, "t3", 60, "+", "UMI3"),
    ]

    unique_transcripts = compute_unique_umis(
        transcripts, umi_counting_offset=10, umi_allowed_mismatches=1, group_umi_func=mock_group_umi_func
    )

    assert len(unique_transcripts) == 1
    assert unique_transcripts[0].umi == "UMI1"

# Test for createDataset with mocked dependencies
def test_create_dataset(tmp_path, monkeypatch, mock_bam_file, mock_gff_file):
    # Mock inputs
    output_folder = tmp_path
    umi_cluster_algorithm = "hierarchical"

    t1 = Transcript("chr1", 100, 200, "t1", 60, "+", "UMI1")
    t2 = Transcript("chr2", 300, 400, "t2", 60, "-", "UMI2")
    # Mock parse_unique_events
    mock_parse_unique_events = Mock(return_value=[
        ("gene1", {(10, 10): [t1, t2]}),
        ("gene2", {(20, 20): [t1, t2]})
    ])
    monkeypatch.setattr("stpipeline.common.dataset.parse_unique_events", mock_parse_unique_events)

    # Mock dedup_hierarchical
    mock_dedup_compute_unique_umis = Mock(return_value=[t1])
    monkeypatch.setattr("stpipeline.common.dataset.compute_unique_umis", mock_dedup_compute_unique_umis)

    stats = createDataset(
        input_file=mock_bam_file,
        output_folder=str(output_folder),
        gff_filename=mock_gff_file,
        umi_cluster_algorithm=umi_cluster_algorithm,
        umi_allowed_mismatches=1,
        umi_counting_offset=10,
        disable_umi=False,
        output_template="output",
        verbose=False
    )

    assert stats["genes_found"] == 2
    assert stats["reads_after_duplicates_removal"] == 2