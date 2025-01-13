#! /usr/bin/env python
"""
Integration tests for each script.
"""

import os
import tempfile
from shutil import copyfile, rmtree
from subprocess import check_call

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def setup_pipeline():
    """Set up temporary directories and required files for the pipeline tests."""
    # Obtain paths and files
    testdir = os.path.abspath("tests")
    infile_fw = os.path.join(testdir, "input/arrayjet_1002/testdata_R1.fastq.gz")
    infile_rv = os.path.join(testdir, "input/arrayjet_1002/testdata_R2.fastq.gz")
    annotfile = os.path.join(testdir, "config/annotations/Homo_sapiens.GRCh38.79_chr19.gtf")
    chipfile = os.path.join(testdir, "config/idfiles/150204_arrayjet_1000L2_probes.txt")

    # Temporary directories for output and logs
    tmpdir = tempfile.mkdtemp(prefix="st_pipeline_test_temp")
    outdir = tempfile.mkdtemp(prefix="st_pipeline_test_output")
    logfile = tempfile.mktemp(prefix="st_pipeline_test_log")

    # Create genome index dirs
    genomedir = os.path.join(tmpdir, "config/genomes/mouse_grcm38")
    contamdir = os.path.join(tmpdir, "config/contaminant_genomes/R45S5_R5S1")
    os.makedirs(genomedir)
    os.makedirs(contamdir)

    genomefasta = os.path.join(genomedir, "human_grcm38_chromosome19.fasta")
    genomefastagz = os.path.join(genomedir, "human_grcm38_chromosome19.fasta.gz")

    # Download and unpack genome files
    copyfile(os.path.join(testdir, "config/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz"), genomefastagz)
    check_call(["gunzip", genomefastagz])

    # Generate genome indexes
    check_call(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeSAindexNbases",
            "11",
            "--runThreadN",
            "4",
            "--genomeDir",
            genomedir,
            "--genomeFastaFiles",
            genomefasta,
        ]
    )

    check_call(
        [
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeSAindexNbases",
            "8",
            "--runThreadN",
            "4",
            "--genomeDir",
            contamdir,
            "--genomeFastaFiles",
            os.path.join(testdir, "config/contaminant_genomes/R45S5_R5S1/Rn45s_Rn5s.fasta"),
        ]
    )

    # Verify existence of input files
    assert os.path.exists(infile_fw)
    assert os.path.exists(infile_rv)
    assert os.path.isdir(genomedir)
    assert os.path.isdir(contamdir)
    assert os.path.exists(annotfile)
    assert os.path.exists(chipfile)
    assert os.path.isdir(outdir)
    assert os.path.isdir(tmpdir)

    # Yield setup data
    yield {
        "infile_fw": infile_fw,
        "infile_rv": infile_rv,
        "annotfile": annotfile,
        "chipfile": chipfile,
        "tmpdir": tmpdir,
        "outdir": outdir,
        "genomedir": genomedir,
        "contamdir": contamdir,
        "logfile": logfile,
    }

    # Cleanup
    rmtree(tmpdir, ignore_errors=True)
    rmtree(outdir, ignore_errors=True)


@pytest.fixture
def test_pipeline_run(setup_pipeline):
    """Test the full pipeline run."""
    data = setup_pipeline

    try:
        check_call(
            [
                "st_pipeline_run",
                "--verbose",
                "--no-clean-up",
                "--star-two-pass-mode",
                "--htseq-no-ambiguous",
                "--keep-discarded-files",
                "--threads",
                "4",
                "--log-file",
                data["logfile"],
                "--expName",
                "test",
                "--ids",
                data["chipfile"],
                "--ref-map",
                data["genomedir"],
                "--contaminant-index",
                data["contamdir"],
                "--output-folder",
                data["outdir"],
                "--temp-folder",
                data["tmpdir"],
                "--ref-annotation",
                data["annotfile"],
                data["infile_fw"],
                data["infile_rv"],
            ]
        )
    except Exception as e:
        pytest.fail(f"st_pipeline_run execution failed: {e}")

    # Verify output files
    datafile = os.path.join(data["outdir"], "test_stdata.tsv")
    readsfile = os.path.join(data["outdir"], "test_reads.bed")
    # statsfile = os.path.join(data["outdir"], "test_qa_stats.json")
    assert os.path.exists(datafile)
    assert os.path.getsize(datafile) > 1024
    assert os.path.exists(readsfile)
    assert os.path.getsize(readsfile) > 1024
    # assert os.path.exists(statsfile)
    return datafile


def test_st_qa(test_pipeline_run, tmpdir):
    """Test the st_qa."""
    datafile = test_pipeline_run

    try:
        check_call(["st_qa", "--outdir", str(tmpdir), datafile])
    except Exception as e:
        pytest.fail(f"st_qa execution failed: {e}")

    clean_name = os.path.basename(datafile).split(".")[0]
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_qa_stats.txt"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_density_genes_by_spot.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_density_spots_by_gene.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_scatter_reads_vs_genes.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_heatmap_counts.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_heatmap_genes.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_reads_spot.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_genes_spot.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_genes_spots_1.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_genes_spots_2.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_spots_gene.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_spots_gene_1.pdf"))
    assert os.path.exists(os.path.join(tmpdir, f"{clean_name}_hist_spots_gene_2.pdf"))


def test_multi_qa(test_pipeline_run, tmpdir):
    """Test the multi_qa."""
    datafile = test_pipeline_run

    try:
        check_call(["multi_qa", "--outdir", str(tmpdir), datafile, datafile])
    except Exception as e:
        pytest.fail(f"multi_qa execution failed: {e}")

    assert os.path.exists(os.path.join(tmpdir, "violin_plot_reads.pdf"))
    assert os.path.exists(os.path.join(tmpdir, "violin_plot_genes.pdf"))
    assert os.path.exists(os.path.join(tmpdir, "gene_correlations.png"))
    assert os.path.exists(os.path.join(tmpdir, "gene_correlations.tsv"))
    assert os.path.exists(os.path.join(tmpdir, "gene_similarities.tsv"))
    assert os.path.exists(os.path.join(tmpdir, "pca.pdf"))


@pytest.fixture
def setup_test_files(tmp_path):
    """Set up fake data for the test."""
    # Create a fake counts matrix
    counts_matrix = tmp_path / "counts_matrix.tsv"
    with counts_matrix.open("w") as f:
        f.write("gene1\tgene2\tgene3\n")
        f.write("1x1\t5\t0\t10\n")
        f.write("2x2\t0\t2\t3\n")
        f.write("3x3\t1\t0\t0\n")
        f.write("4x4\t0\t0\t0\n")

    # Create a fake coordinates file
    coordinates_file = tmp_path / "coordinates.tsv"
    with coordinates_file.open("w") as f:
        f.write("x\ty\tnew_x\tnew_y\n")
        f.write("1\t1\t10\t10\n")
        f.write("2\t2\t20\t20\n")
        f.write("3\t3\t30\t30\n")
        # Spot 4x4 will not be included because it's not in the coordinates file

    # Define the output file
    outfile = tmp_path / "adjusted_counts_matrix.tsv"

    return counts_matrix, coordinates_file, outfile


def test_adjust_matrix_coordinates(setup_test_files, tmp_path):
    counts_matrix, coordinates_file, outfile = setup_test_files

    # Call the adjust_matrix_coordinates script
    try:
        check_call(
            [
                "adjust_matrix_coordinates",
                str(counts_matrix),
                "--coordinates-file",
                str(coordinates_file),
                "--update-coordinates",
                "--outfile",
                str(outfile),
            ]
        )
    except Exception as e:
        pytest.fail(f"adjust_matrix_coordinates execution failed: {e}")

    # Verify the output file
    assert outfile.exists(), "Output file should be created"

    # Load the output matrix
    output_matrix = pd.read_csv(outfile, sep="\t", index_col=0)

    # Check that the output contains the expected data
    expected_data = {"gene1": [5, 0, 1], "gene2": [0, 2, 0], "gene3": [10, 3, 0]}
    expected_index = ["10.0x10.0", "20.0x20.0", "30.0x30.0"]  # Adjusted coordinates
    expected_df = pd.DataFrame(expected_data, index=expected_index)

    pd.testing.assert_frame_equal(output_matrix, expected_df, check_dtype=False)

    # Ensure genes with total count zero are removed
    assert "4x4" not in output_matrix.index, "Spots with no coordinates should be removed"


@pytest.fixture
def setup_fastq_files(tmp_path):
    """Set up fake FASTQ files and directories for the test."""
    # Create a run directory with mock FASTQ files
    run_path = tmp_path / "run"
    run_path.mkdir()

    # Create output directory
    out_path = tmp_path / "output"
    out_path.mkdir()

    # Create mock FASTQ files for two identifiers
    identifiers = ["S1", "S2"]
    for idx in identifiers:
        for r in ["R1", "R2"]:
            with open(run_path / f"{idx}_{r}_001.fastq", "w") as f:
                f.write(f"@SEQ_ID\n{idx}{r}SEQUENCE\n+\nIII\n")

    return run_path, out_path, identifiers


def test_merge_fastq(setup_fastq_files, tmp_path):
    run_path, out_path, identifiers = setup_fastq_files

    # Run the script
    try:
        check_call(
            ["merge_fastq", "--run-path", str(run_path), "--out-path", str(out_path), "--identifiers", *identifiers]
        )
    except Exception as e:
        pytest.fail(f"merge_fastq execution failed: {e}")

    # Verify that the merged files are created in the output directory
    for idx in identifiers:
        for r in ["R1", "R2"]:
            merged_file = out_path / f"{idx}_{r}.fastq.gz"
            assert merged_file.exists(), f"Merged file {merged_file} should exist"
            assert merged_file.stat().st_size > 0, f"Merged file {merged_file} should not be empty"


@pytest.fixture
def setup_filter_gene_type_data(tmp_path):
    """
    Fixture to create fake data for the filter_gene_type_matrix script.
    """
    # Create a temporary counts matrix file
    counts_matrix = tmp_path / "counts_matrix.tsv"
    counts_data = {
        "gene1": [10, 0],
        "gene2": [0, 0],
        "gene3": [5, 15],
    }
    counts_df = pd.DataFrame(counts_data, index=["1x1", "2x2"])
    counts_df.to_csv(counts_matrix, sep="\t")

    # Create a temporary annotation file
    annotation_file = tmp_path / "annotation.gtf"
    content = (
        "##gff-version 3\n"
        "chr1\tsource\tfeature\t100\t200\t.\t+\t.\tgene_id=gene1;gene_name=gene1;gene_type=protein_coding\n"
        "chr1\tsource\tfeature\t300\t400\t.\t-\t.\tgene_id=gene2;gene_name=gene2;gene_type=lincRNA\n"
        "chr2\tsource\tfeature\t500\t600\t.\t+\t.\tgene_id=gene3;gene_name=gene3;gene_type=protein_coding\n"
    )
    with open(annotation_file, "w") as f:
        f.write(content)

    # Output file path
    outfile = tmp_path / "filtered_counts_matrix.tsv"

    return counts_matrix, annotation_file, outfile


def test_filter_gene_type_matrix(setup_filter_gene_type_data, tmp_path):
    counts_matrix, annotation_file, outfile = setup_filter_gene_type_data

    # Run the script
    try:
        check_call(
            [
                "filter_gene_type_matrix",
                str(counts_matrix),
                "--annotation",
                str(annotation_file),
                "--gene-types-keep",
                "protein_coding",
                "--outfile",
                str(outfile),
            ]
        )
    except Exception as e:
        pytest.fail(f"filter_gene_type_matrix execution failed: {e}")

    # Validate the output file
    assert outfile.exists(), "Filtered counts matrix should exist"
    filtered_df = pd.read_csv(outfile, sep="\t", index_col=0)

    # Validate the filtered data
    assert "gene1" in filtered_df.columns.to_list(), "gene1 should be present in the filtered data"
    assert "gene3" in filtered_df.columns.to_list(), "gene3 should be present in the filtered data"
    assert "gene2" not in filtered_df.columns.to_list(), "gene2 should be filtered out"


@pytest.fixture
def setup_convert_ensembl_data(tmp_path):
    """
    Fixture to create fake data for the convertEnsemblToNames script.
    """
    # Create a temporary counts matrix file
    counts_matrix = tmp_path / "counts_matrix.tsv"
    counts_data = {
        "ENSG000001": [10, 0],
        "ENSG000002": [0, 0],
        "ENSG000003": [5, 15],
    }
    counts_df = pd.DataFrame(counts_data, index=["1x1", "2x2"])
    counts_df.to_csv(counts_matrix, sep="\t")

    # Create a temporary annotation file
    annotation_file = tmp_path / "annotation.gtf"
    content = (
        "##gff-version 3\n"
        "chr1\tsource\tfeature\t100\t200\t.\t+\t.\tgene_id=ENSG000001;gene_name=Gene1\n"
        "chr1\tsource\tfeature\t300\t400\t.\t-\t.\tgene_id=ENSG000002;gene_name=Gene2\n"
        "chr2\tsource\tfeature\t500\t600\t.\t+\t.\tgene_id=ENSG000003;gene_name=Gene3\n"
    )
    with open(annotation_file, "w") as f:
        f.write(content)

    # Output file path
    output_file = tmp_path / "output_counts_matrix.tsv"

    return counts_matrix, annotation_file, output_file


def test_convert_ensembl_to_names(setup_convert_ensembl_data):
    counts_matrix, annotation_file, output_file = setup_convert_ensembl_data

    # Run the script
    try:
        check_call(
            [
                "convertEnsemblToNames",
                str(counts_matrix),
                "--annotation",
                str(annotation_file),
                "--output",
                str(output_file),
            ]
        )
    except Exception as e:
        pytest.fail(f"convertEnsemblToNames execution failed: {e}")

    # Validate the output file
    assert output_file.exists(), "Converted counts matrix should exist"
    converted_df = pd.read_csv(output_file, sep="\t", index_col=0)

    # Validate the converted data
    assert "Gene1" in converted_df.columns.to_list(), "GeneA should be present in the converted data"
    assert "Gene2" in converted_df.columns.to_list(), "GeneB should be present in the converted data"
    assert "Gene3" in converted_df.columns.to_list(), "GeneC should be present in the converted data"
    assert "ENSG000001" not in converted_df.columns.to_list(), "ENSG000001 should not be present in the converted data"
    assert "ENSG000002" not in converted_df.columns.to_list(), "ENSG000002 should not be present in the converted data"
    assert "ENSG000003" not in converted_df.columns.to_list(), "ENSG000003 should not be present in the converted data"
