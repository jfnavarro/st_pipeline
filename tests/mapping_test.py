#! /usr/bin/env python
"""
Unit-test the package mapping
"""
import subprocess
from unittest.mock import patch, MagicMock, mock_open
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing


def test_alignReads():
    mock_log_content = """\
        Number of input reads 100000
        Average input read length 150
        Uniquely mapped reads number 90000
        Uniquely mapped reads 90
        Number of reads mapped to multiple loci 5000
        % of reads mapped to multiple loci 5
        % of reads unmapped: too short 5
    """

    with patch("subprocess.Popen") as mock_popen, patch("stpipeline.core.mapping.file_ok", return_value=True), patch(
        "stpipeline.core.mapping.shutil.move"
    ) as mock_shutil_move, patch(
        "stpipeline.core.mapping.open", mock_open(read_data=mock_log_content)
    ) as mock_open_file:
        # Mock the subprocess to simulate STAR execution
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"", b"")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        # Mock the log file content
        mock_open_file.return_value.__enter__.return_value.read.return_value = mock_log_content

        # Call the function
        total_reads = alignReads(
            reverse_reads="test.bam",
            ref_map="ref",
            outputFile="out.bam",
            annotation=None,
            outputFolder="output",
            trimReverse=10,
            invTrimReverse=10,
            cores=4,
            min_intron_size=20,
            max_intron_size=1000,
            disable_multimap=False,
            diable_softclipping=False,
            twopassMode=True,
            min_length=50,
            include_non_mapped=True,
            star_genome_loading="NoSharedMemory",
            star_sort_mem_limit=64000000,
        )

        # Assertions for subprocess call
        mock_popen.assert_called_once()
        mock_shutil_move.assert_called_once_with("output/Aligned.sortedByCoord.out.bam", "out.bam")

        expected_args = [
            "STAR",
            "--genomeDir",
            "ref",
            "--readFilesIn",
            "test.bam",
            "--outFileNamePrefix",
            "output/",
            "--clip3pNbases",
            "10",
            "--clip5pNbases",
            "10",
            "--runThreadN",
            "4",
            "--outFilterType",
            "Normal",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--alignEndsType",
            "Local",
            "--outSAMorder",
            "Paired",
            "--outSAMprimaryFlag",
            "OneBestScore",
            "--outFilterMultimapNmax",
            "20",
            "--alignIntronMin",
            "20",
            "--alignIntronMax",
            "1000",
            "--outFilterMatchNmin",
            "50",
            "--genomeLoad",
            "NoSharedMemory",
            "--limitBAMsortRAM",
            "64000000",
            "--readFilesType",
            "SAM",
            "SE",
            "--readFilesCommand",
            "samtools",
            "view",
            "-h",
            "--twopassMode",
            "Basic",
            "--outSAMunmapped",
            "Within",
        ]

        mock_popen.assert_called_once_with(
            expected_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True, shell=False
        )

        # Ensure the log file was read
        mock_open_file.assert_called_once_with("output/Log.final.out", "r")

        # Log file parsing validation
        assert total_reads == 5090


def test_barcodeDemultiplexing():
    with patch("subprocess.Popen") as mock_popen, patch("os.path.isfile", return_value=True), patch(
        "stpipeline.core.mapping.file_ok", return_value=True
    ):
        mock_process = MagicMock()
        mock_process.communicate.return_value = (b"Total reads: 100\nTotal reads written: 80", b"")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process

        total_reads = barcodeDemultiplexing(
            reads="reads.bam",
            idFile="barcodes.tsv",
            mismatches=1,
            kmer=8,
            over_hang=2,
            taggd_metric="Levenshtein",
            taggd_multiple_hits_keep_one=True,
            taggd_trim_sequences=[1, 2, 3],
            cores=4,
            outputFilePrefix="output/test",
            keep_discarded_files=False,
        )

        expected_args = [
            "taggd_demultiplex.py",
            "--max-edit-distance",
            "1",
            "--k",
            "8",
            "--barcode-tag",
            "B0",
            "--homopolymer-filter",
            "0",
            "--subprocesses",
            "4",
            "--metric",
            "Levenshtein",
            "--overhang",
            "2",
            "--trim-sequences",
            "1",
            "2",
            "3",
            "--multiple-hits-keep-one",
            "--no-unmatched-output",
            "--no-ambiguous-output",
            "--no-results-output",
            "barcodes.tsv",
            "reads.bam",
            "output/test",
        ]

        mock_popen.assert_called_once_with(
            expected_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True, shell=False
        )
        assert total_reads == 80
