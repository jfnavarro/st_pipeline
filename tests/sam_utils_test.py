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
def mock_bam_files(tmp_path):
    """Creates mock BAM files for testing."""
    bam_files = []
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 2000, "SN": "chr1"}]}

    for i in range(3):
        bam_file = tmp_path / f"mock{i}.bam"
        bam_files.append(str(bam_file))
        with pysam.AlignmentFile(bam_file, "wb", header=header) as f:
            for j in range(5):
                segment = pysam.AlignedSegment()
                segment.query_name = f"read{i}_{j}"
                segment.query_sequence = "ACTG" * 25
                segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
                segment.flag = 0
                segment.reference_id = 0
                segment.reference_start = j * 100
                segment.cigar = [(0, len(segment.query_sequence))]  # Match
                f.write(segment)

    return bam_files


@pytest.fixture
def mock_bam_file_for_split(tmp_path):
    """Creates a mock BAM file for testing split_bam."""
    bam_file = tmp_path / "mock.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 2000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_file, "wb", header=header) as f:
        for i in range(20):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{i}"
            segment.query_sequence = "ACTG" * 25
            segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = i * 100
            segment.cigar = [(0, len(segment.query_sequence))]  # Match
            f.write(segment)

    return str(bam_file)


def test_split_bam(mock_bam_file_for_split, tmp_path):
    """Test the split_bam function with a mocked BAM file."""
    threads = 4
    temp_dir = tmp_path / "split_bam_output"
    temp_dir.mkdir()

    # Call the split_bam function
    split_files = split_bam(mock_bam_file_for_split, str(temp_dir), threads)

    # Verify the number of split files
    assert len(split_files) == threads

    # Verify the contents of each split file
    total_records = 0
    for part, split_file in enumerate(split_files):
        with pysam.AlignmentFile(split_file, "rb") as f:
            records = list(f.fetch(until_eof=True))
            total_records += len(records)

            # Assert each file has roughly equal reads
            if part < threads - 1:
                assert len(records) == 5  # 20 reads / 4 threads = 5 reads per part

    # Verify the total records match the original file
    assert total_records == 20


@patch("stpipeline.common.sam_utils.pysam.AlignmentFile")
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


def test_merge_bam(mock_bam_files, tmp_path):
    merged_file = tmp_path / "merged.bam"

    # Call the merge_bam function
    total_records = merge_bam(str(merged_file), mock_bam_files)

    # Assert the total record count
    assert total_records == 15  # 3 files x 5 records each

    # Verify the contents of the merged BAM file
    with pysam.AlignmentFile(merged_file, "rb") as f:
        records = list(f.fetch(until_eof=True))
        assert len(records) == 15
        for i, record in enumerate(records):
            assert record.query_name == f"read{i // 5}_{i % 5}"
            assert record.query_sequence == "ACTG" * 25
            assert record.flag == 0
            assert record.reference_start == (i % 5) * 100
