#! /usr/bin/env python
""" 
Unit-test the package dataset
"""
import pytest
from typing import List
from stpipeline.common.dataset import Transcript, compute_unique_umis, createDataset
from collections import Counter
from unittest.mock import Mock

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
def test_create_dataset(tmp_path, monkeypatch):
    # Mock inputs
    input_file = tmp_path / "test_input.bam"
    gff_filename = tmp_path / "test.gff"
    output_folder = tmp_path
    umi_cluster_algorithm = "hierarchical"

    input_file.write_text("")  # Create an empty mock file

    # Mock parse_unique_events
    mock_parse_unique_events = Mock(return_value=[
        ("gene1", {(10, 10): [Transcript("chr1", 100, 200, "t1", 60, "+", "UMI1")]}),
        ("gene2", {(20, 20): [Transcript("chr2", 300, 400, "t2", 60, "-", "UMI2")]})
    ])
    monkeypatch.setattr("your_module_name.parse_unique_events", mock_parse_unique_events)

    # Mock dedup_hierarchical
    mock_dedup_hierarchical = Mock(return_value=["UMI1"])
    monkeypatch.setattr("stpipeline.common.clustering.dedup_hierarchical", mock_dedup_hierarchical)

    stats = createDataset(
        input_file=str(input_file),
        gff_filename=str(gff_filename),
        umi_cluster_algorithm=umi_cluster_algorithm,
        umi_allowed_mismatches=1,
        umi_counting_offset=10,
        disable_umi=False,
        output_folder=str(output_folder),
        output_template="output",
        verbose=False
    )

    assert stats["genes_found"] == 2
    assert stats["reads_after_duplicates_removal"] == 2
    assert "UMI1" in stats