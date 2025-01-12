#! /usr/bin/env python
"""
Unit-test the package stats
"""
import pytest
import os
import json
from stpipeline.common.stats import Stats


@pytest.fixture
def sample_stats():
    return Stats(
        input_reads_forward=1000,
        input_reads_reverse=900,
        reads_after_trimming_forward=800,
        reads_after_trimming_reverse=850,
        reads_after_rRNA_trimming=700,
        reads_after_mapping=600,
        reads_after_annotation=550,
        reads_after_demultiplexing=500,
        reads_after_duplicates_removal=450,
        genes_found=100,
        duplicates_found=50,
        pipeline_version="1.0",
        mapper_tool="bwa",
        annotation_tool="gffread",
        demultiplex_tool="umi_tools",
        input_parameters=["--trim", "--map"],
        max_genes_feature=50,
        min_genes_feature=5,
        max_reads_feature=200,
        min_reads_feature=10,
        average_gene_feature=25.0,
        average_reads_feature=100.0,
    )


def test_stats_str(sample_stats):
    stats_str = str(sample_stats)
    assert "input_reads_forward: 1000" in stats_str
    assert "pipeline_version: 1.0" in stats_str
    assert "average_reads_feature: 100.0" in stats_str


def test_write_json(sample_stats, tmp_path):
    json_file = tmp_path / "stats.json"
    sample_stats.write_json(str(json_file))

    assert os.path.exists(json_file)

    with open(json_file, "r") as file:
        data = json.load(file)
        assert data["input_reads_forward"] == 1000
        assert data["pipeline_version"] == "1.0"
        assert data["average_reads_feature"] == 100.0


def test_from_json(tmp_path):
    json_file = tmp_path / "stats.json"
    sample_data = {
        "input_reads_forward": 1000,
        "input_reads_reverse": 900,
        "reads_after_trimming_forward": 800,
        "reads_after_trimming_reverse": 850,
        "reads_after_rRNA_trimming": 700,
        "reads_after_mapping": 600,
        "reads_after_annotation": 550,
        "reads_after_demultiplexing": 500,
        "reads_after_duplicates_removal": 450,
        "genes_found": 100,
        "duplicates_found": 50,
        "pipeline_version": "1.0",
        "mapper_tool": "bwa",
        "annotation_tool": "gffread",
        "demultiplex_tool": "umi_tools",
        "input_parameters": ["--trim", "--map"],
        "max_genes_feature": 50,
        "min_genes_feature": 5,
        "max_reads_feature": 200,
        "min_reads_feature": 10,
        "average_gene_feature": 25.0,
        "average_reads_feature": 100.0,
    }

    with open(json_file, "w") as file:
        json.dump(sample_data, file)

    stats = Stats.from_json(str(json_file))
    assert stats.input_reads_forward == 1000
    assert stats.pipeline_version == "1.0"
    assert stats.average_reads_feature == 100.0
