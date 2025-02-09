#! /usr/bin/env python
"""
Unit-test the package saturation
"""

import os
from unittest.mock import patch

import pysam
import pytest

from stpipeline.common.saturation import (
    _cleanup_files,
    _compute_saturation_metrics,
    _determine_saturation_points,
    _generate_subsamples,
    _write_subsamples_to_files,
)


@pytest.fixture
def mock_bam_file(tmp_path):
    """
    Creates a mock BAM file for testing.

    Args:
        tmp_path: pytest fixture for temporary directory.

    Returns:
        Path to the mock BAM file and number of reads.
    """
    bam_path = tmp_path / "mock.bam"
    with pysam.AlignmentFile(
        bam_path, "wb", header={"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    ) as bam_file:
        for i in range(100):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{i}"
            segment.query_sequence = "ACTG" * 25
            segment.query_qualities = pysam.qualitystring_to_array("IIII" * 25)
            segment.flag = 0
            segment.reference_id = 0
            segment.reference_start = i * 10
            segment.cigar = [(0, len(segment.query_sequence))]  # 0: MATCH
            segment.set_tag("B1", i)
            segment.set_tag("B2", i * 2)
            segment.set_tag("XF", "gene1")
            segment.set_tag("B3", "UMI1")
            bam_file.write(segment)
    return str(bam_path), 100


def test_determine_saturation_points():
    nreads = 10000
    saturation_points = [100, 500, 1000, 20000]
    points = _determine_saturation_points(nreads, saturation_points)
    assert points == [100, 500, 1000]

    # Test with None
    points = _determine_saturation_points(nreads, None)
    assert len(points) > 0
    assert all(p < nreads for p in points)


def test_generate_subsamples(mock_bam_file, tmp_path):
    bam_file, nreads = mock_bam_file
    saturation_points = [10, 50, 100]
    temp_folder = tmp_path

    files, file_names, subsampling = _generate_subsamples(nreads, bam_file, saturation_points, temp_folder)

    assert len(files) == len(saturation_points)
    assert len(file_names) == len(saturation_points)
    for spoint in saturation_points:
        assert spoint in subsampling
        assert len(subsampling[spoint]) == spoint

    # Cleanup
    for file in files.values():
        file.close()


def test_write_subsamples_to_files(mock_bam_file, tmp_path):
    bam_file, nreads = mock_bam_file
    saturation_points = [10, 50, 100]
    temp_folder = tmp_path

    files, file_names, subsampling = _generate_subsamples(nreads, bam_file, saturation_points, temp_folder)
    _write_subsamples_to_files(files, subsampling, bam_file, saturation_points)

    for spoint, file_name in file_names.items():
        with pysam.AlignmentFile(file_name, "rb") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
            assert count == spoint

    # Cleanup
    for file_name in file_names.values():
        os.remove(file_name)


def test_compute_saturation_metrics(mock_bam_file, tmp_path):
    bam_file, nreads = mock_bam_file
    saturation_points = [10, 50, 100]
    temp_folder = tmp_path
    gff_filename = tmp_path / "mock.gff"
    gff_filename.write_text("chr1\tsource\tfeature\t1\t1000\t.\t+\t.\tgene_id=gene1\n")

    files, file_names, subsampling = _generate_subsamples(nreads, bam_file, saturation_points, temp_folder)
    _write_subsamples_to_files(files, subsampling, bam_file, saturation_points)

    with patch("stpipeline.common.dataset.createDataset") as mock_createDataset:
        mock_createDataset.return_value = {
            "reads_after_duplicates_removal": 10,
            "genes_found": 5,
            "average_gene_feature": 2.5,
            "average_reads_feature": 1.0,
        }

        results = _compute_saturation_metrics(
            file_names,
            saturation_points,
            str(gff_filename),
            "AdjacentBi",
            True,
            10,
            False,
            str(temp_folder),
            "test_exp",
        )

        assert len(results["reads"]) == len(saturation_points)
        assert len(results["genes"]) == len(saturation_points)

    # Cleanup
    for file_name in file_names.values():
        os.remove(file_name)


def test_cleanup_files(tmp_path):
    temp_files = [tmp_path / f"temp_file_{i}.bam" for i in range(5)]
    for file in temp_files:
        file.touch()
    _cleanup_files({i: str(file) for i, file in enumerate(temp_files)})
    for file in temp_files:
        assert not os.path.exists(file)
