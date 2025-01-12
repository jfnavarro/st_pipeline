#! /usr/bin/env python
"""
Unit-test the package filter
"""
import pytest
import gzip
from stpipeline.common.gff_reader import gff_lines, gff_parse


# Mock GTF/GFF data for testing
def generate_mock_gff_file(file_path: str, content: str):
    """
    Generates a mock GFF file for testing.

    Args:
        file_path (str): Path to the file to create.
        content (str): Content to write into the file.
    """
    with open(file_path, "w") as f:
        f.write(content)


@pytest.fixture
def mock_gff_file(tmp_path):
    content = (
        "##gff-version 3\n"
        "chr1\tsource\tfeature\t100\t200\t.\t+\t.\tID=gene1;Name=Gene1\n"
        "chr1\tsource\tfeature\t300\t400\t.\t-\t.\tID=gene2;Name=Gene2\n"
        "chr2\tsource\tfeature\t500\t600\t.\t+\t.\tID=gene3;Name=Gene3\n"
    )
    file_path = tmp_path / "test.gff"
    generate_mock_gff_file(file_path, content)
    return str(file_path)


@pytest.fixture
def mock_gzipped_gff_file(tmp_path):
    content = (
        "##gff-version 3\n"
        "chr1\tsource\tfeature\t100\t200\t.\t+\t.\tID=gene1;Name=Gene1\n"
        "chr1\tsource\tfeature\t300\t400\t.\t-\t.\tID=gene2;Name=Gene2\n"
        "chr2\tsource\tfeature\t500\t600\t.\t+\t.\tID=gene3;Name=Gene3\n"
    )
    file_path = tmp_path / "test.gff.gz"
    with gzip.open(file_path, "wt") as f:
        f.write(content)
    return str(file_path)


def test_gff_lines_plain(mock_gff_file):
    parsed_lines = list(gff_lines(mock_gff_file))
    assert len(parsed_lines) == 3
    assert parsed_lines[0]["seqname"] == "chr1"
    assert parsed_lines[0]["start"] == "100"
    assert parsed_lines[0]["end"] == "200"
    assert parsed_lines[0]["ID"] == "gene1"


def test_gff_lines_gzipped(mock_gzipped_gff_file):
    parsed_lines = list(gff_lines(mock_gzipped_gff_file))
    assert len(parsed_lines) == 3
    assert parsed_lines[0]["seqname"] == "chr1"
    assert parsed_lines[0]["start"] == "100"
    assert parsed_lines[0]["end"] == "200"
    assert parsed_lines[0]["ID"] == "gene1"


def test_gff_parse():
    line = "chr1\tsource\tfeature\t100\t200\t.\t+\t.\tID=gene1;Name=Gene1"
    parsed = gff_parse(line)
    assert parsed["seqname"] == "chr1"
    assert parsed["start"] == "100"
    assert parsed["end"] == "200"
    assert parsed["ID"] == "gene1"
    assert parsed["Name"] == "Gene1"
