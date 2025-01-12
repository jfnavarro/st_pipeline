#! /usr/bin/env python
""" 
Unit-test the package filter
"""
import pytest
import dnaio
from unittest.mock import Mock, patch
from stpipeline.common.filter import filter_input_data

def generate_test_fastq(filepath, records):
    """
    Generates a mock FASTQ file for testing.

    Args:
        filepath (str): Path to the file to create.
        records (list): List of tuples (header, sequence, quality) for the FASTQ records.
    """
    with open(filepath, "w") as f:
        for header, sequence, quality in records:
            f.write(f"@{header}\n{sequence}\n+\n{quality}\n")

@pytest.fixture
def setup_fastq_files(tmp_path):
    fw_records = [
        ("read1", "ACTGACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read2", "TTTTTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read3", "GGGGGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read4", "CCCCCCCCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read5", "ACTGACTGACTGACTGACTGACTG", "!!!!IIIIIIIIIIIIIIIIIIII"),  # Low-quality UMI
        ("read6", "ACTGACTGACTGACTGACTGACTG", "!!!!!!!!!!!!!!IIIIIIIIII")  # Too short after trimming
    ]
    rv_records = [
        ("read1", "ACTGACTGACTGACTGACTGACTG", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read2", "TTTTTTTTTTTTTTTTTTTTTTTT", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read3", "GGGGGGGGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read4", "CCCCCCCCCCCCCCCCCCCCCCCC", "IIIIIIIIIIIIIIIIIIIIIIII"),
        ("read5", "ACTGACTGACTGACTGACTGACTG", "!!!!IIIIIIIIIIIIIIIIIIII"),  # Low-quality UMI
        ("read6", "ACTGACTGACTGACTGACTGACTG", "!!!!!!!!!!!!!!IIIIIIIIII")  # Too short after trimming
    ]
    fw_file = tmp_path / "fw.fastq"
    rv_file = tmp_path / "rv.fastq"

    generate_test_fastq(fw_file, fw_records)
    generate_test_fastq(rv_file, rv_records)

    return str(fw_file), str(rv_file)

@patch("stpipeline.common.filter.pysam.AlignmentFile")
def test_filter_input_data(mock_alignment_file, setup_fastq_files, tmp_path):
    fw_file, rv_file = setup_fastq_files
    out_file = tmp_path / "output.bam"
    out_file_discarded = tmp_path / "discarded.fastq"

    mock_alignment_file.return_value.__enter__.return_value = Mock()

    total_reads, remaining_reads = filter_input_data(
        fw_file=fw_file,
        rv_file=rv_file,
        out_file=str(out_file),
        out_file_discarded=str(out_file_discarded),
        barcode_length=10,
        start_position=0,
        filter_AT_content=50.0,
        filter_GC_content=50.0,
        umi_start=0,
        umi_end=4,
        min_qual=20,
        min_length=20,
        polyA_min_distance=5,
        polyT_min_distance=5,
        polyG_min_distance=5,
        polyC_min_distance=5,
        polyN_min_distance=5,
        qual64=False,
        umi_filter=True,
        umi_filter_template=r"[ACGT]{4}",
        umi_quality_bases=1,
        adaptor_missmatches=2,
        overhang=2,
        disable_umi=False,
        disable_barcode=False
    )

    assert total_reads == 6
    assert remaining_reads < total_reads
    mock_alignment_file.assert_called_once_with(str(out_file), "wb")