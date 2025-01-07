#! /usr/bin/env python
""" 
Unit-test the package sam_utils
"""
import pytest
import os
import pysam
from unittest.mock import Mock, patch
from stpipeline.common.sam_utils import split_bam, convert_to_AlignedSegment, merge_bam

@pytest.fixture
def mock_bam_file(tmp_path):
    """
    Creates a mock BAM file for testing.

    Args:
        tmp_path: pytest fixture for temporary directory.

    Returns:
        Path to the mock BAM file.
    """
    bam_path = tmp_path / "mock.bam"
    with pysam.AlignmentFile(bam_path, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}) as bam_file:
        for i in range(100):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{i}"
            segment.query_sequence = "ACTG" * 25
            segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = i * 10
            bam_file.write(segment)
    return str(bam_path)

@patch("your_module_name.pysam.AlignmentFile")
def test_split_bam(mock_alignment_file, mock_bam_file, tmp_path):
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    threads = 4

    # Mock pysam.AlignmentFile
    mock_alignment_file.return_value.__enter__.return_value = Mock()

    split_files = split_bam(mock_bam_file, str(temp_dir), threads)

    assert len(split_files) == threads
    for split_file in split_files:
        assert os.path.exists(split_file)

@patch("your_module_name.pysam.AlignmentFile")
def test_convert_to_aligned_segment(mock_alignment_file):
    header = "read1"
    sequence = "ACTGACTGACTG"
    quality = "IIIIIIIIIIII"
    barcode_sequence = "ACGT"
    umi_sequence = "TGCA"

    aligned_segment = convert_to_AlignedSegment(header, sequence, quality, barcode_sequence, umi_sequence)

    assert aligned_segment.query_name == "read1"
    assert aligned_segment.query_sequence == sequence
    assert aligned_segment.query_qualities.tolist() == [40] * len(sequence)
    assert aligned_segment.get_tag("B0") == barcode_sequence
    assert aligned_segment.get_tag("B3") == umi_sequence

@patch("your_module_name.pysam.AlignmentFile")
def test_merge_bam(mock_alignment_file, tmp_path):
    bam1 = tmp_path / "part1.bam"
    bam2 = tmp_path / "part2.bam"
    merged_bam = tmp_path / "merged.bam"

    # Mock BAM files
    with pysam.AlignmentFile(bam1, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}) as f:
        pass
    with pysam.AlignmentFile(bam2, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}) as f:
        pass

    mock_alignment_file.return_value.__enter__.return_value = Mock()

    total_records = merge_bam(str(merged_bam), [str(bam1), str(bam2)])

    assert os.path.exists(merged_bam)
    assert total_records > 0
