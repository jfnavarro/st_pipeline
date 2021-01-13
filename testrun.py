#! /usr/bin/env python
""" 
Unit-test for run-tests, it just tests that the pipeline runs and produces correct results
"""

import unittest
import urllib
import tempfile
import multiprocessing
import pandas as pd
import numpy as np
from subprocess import check_call
from stpipeline.core.pipeline import Pipeline
import os
from shutil import copyfile


class TestPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        # Obtain paths and files.
        testdir = os.path.abspath("tests")
        self.infile_fw = os.path.join(testdir, 'input/arrayjet_1002/testdata_R1.fastq')
        self.infile_rv = os.path.join(testdir, 'input/arrayjet_1002/testdata_R2.fastq')
        self.annotfile = os.path.join(testdir, 'config/annotations/Homo_sapiens.GRCh38.79_chr19.gtf')
        self.chipfile = os.path.join(testdir, 'config/idfiles/150204_arrayjet_1000L2_probes.txt')
        self.expname = "test"

        # Obtain temp dir
        self.tmpdir = tempfile.mkdtemp(prefix="st_pipeline_test_temp")
        print("ST Pipeline Test Temporary directory {}".format(self.tmpdir))
        self.outdir = tempfile.mkdtemp(prefix="st_pipeline_test_output")
        print("ST Pipeline Test Temporary output {}".format(self.outdir))
        self.error_file = os.path.join(self.tmpdir, 'error.log')
        self.logFile = tempfile.mktemp(prefix="st_pipeline_test_log")
        print("ST Pipeline Test Log file {}".format(self.logFile))

        # Create genome index dirs.
        self.genomedir = os.path.join(self.tmpdir, 'config/genomes/mouse_grcm38')
        os.makedirs(self.genomedir)

        # STAR contaminant dir
        self.contamdir = os.path.join(self.tmpdir, 'config/contaminant_genomes/R45S5_R5S1')
        os.makedirs(self.contamdir)

        genomefasta = os.path.join(self.genomedir, "human_grcm38_chromosome19.fasta")
        genomefastagz = os.path.join(self.genomedir, "human_grcm38_chromosome19.fasta.gz")

        # Change dir to the temp folder
        os.chdir(self.tmpdir)

        # Download and unpack fasta files
        try:
            print("ST Pipeline Test Downloading genome files...")
            copyfile(os.path.join(testdir, "config/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz"),
                     genomefastagz)
            check_call(['gunzip', genomefastagz])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Downloading genome files failed \n")

        # Make genome indexes
        try:
            print("ST Pipeline Test Creating genome index...")
            check_call(["STAR",
                        "--runMode", "genomeGenerate",
                        "--genomeSAindexNbases", "11",
                        "--runThreadN", "4",
                        "--genomeDir", self.genomedir,
                        "--genomeFastaFiles", genomefasta])

            print("ST Pipeline Test Creating contaminant genome index...")
            contamfasta = os.path.join(testdir, "config/contaminant_genomes/R45S5_R5S1/Rn45s_Rn5s.fasta")
            check_call(["STAR",
                        "--runMode", "genomeGenerate",
                        "--genomeSAindexNbases", "8",
                        "--runThreadN", "4",
                        "--genomeDir", self.contamdir,
                        "--genomeFastaFiles", contamfasta])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Creating genome index failed \n")

        # Verify existence of input files
        assert (os.path.exists(self.infile_fw))
        assert (os.path.exists(self.infile_rv))
        assert (os.path.isdir(self.genomedir))
        assert (os.path.isdir(self.contamdir))
        assert (os.path.exists(self.annotfile))
        assert (os.path.exists(self.chipfile))
        assert (os.path.isdir(self.outdir))
        assert (os.path.isdir(self.tmpdir))

    @classmethod
    def tearDownClass(self):
        print("ST Pipeline Test Remove temporary output {}".format(self.outdir))
        for root, dirs, files in os.walk(self.outdir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        if os.path.exists(self.outdir):
            os.rmdir(self.outdir)

        print("ST Pipeline Test Remove temporary directory {}".format(self.tmpdir))
        for root, dirs, files in os.walk(self.tmpdir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        if os.path.exists(self.tmpdir):
            os.rmdir(self.tmpdir)

            # Remove STAR log files
        log_std = "Log.std.out"
        log = "Log.out"
        log_sj = "SJ.out.tab"
        log_final = "Log.final.out"
        log_progress = "Log.progress.out"
        if os.path.isfile(log_std):
            os.remove(log_std)
        if os.path.isfile(log):
            os.remove(log)
        if os.path.isfile(log_sj):
            os.remove(log_sj)
        if os.path.isfile(log_progress):
            os.remove(log_progress)
        if os.path.isfile(log_final):
            os.remove(log_final)

    def test_run(self):
        # Run st_pipeline_run.py
        try:
            print("Running st_pipeline_run.py")
            check_call(["st_pipeline_run.py",
                        "--verbose",
                        "--no-clean-up",
                        "--star-two-pass-mode",
                        "--htseq-no-ambiguous",
                        "--keep-discarded-files",
                        "--threads", "4",
                        "--log-file", self.logFile,
                        "--expName", self.expname,
                        "--ids", self.chipfile,
                        "--ref-map", self.genomedir,
                        "--contaminant-index", self.contamdir,
                        "--output-folder", self.outdir,
                        "--temp-folder", self.tmpdir,
                        "--ref-annotation", self.annotfile,
                        self.infile_fw, self.infile_rv])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running st_pipeline_run.py failed \n")

        # Verify existence of output files and temp files
        self.assertNotEqual(os.listdir(self.outdir), [], "Output folder is not empty")
        self.assertNotEqual(os.listdir(self.tmpdir), [], "Tmp folder is not empty")
        datafile = os.path.join(self.outdir, self.expname + "_stdata.tsv")
        readsfile = os.path.join(self.outdir, self.expname + "_reads.bed")
        # statsfile = os.path.join(self.outdir, self.expname + "_qa_stats.json")
        self.assertTrue(os.path.exists(datafile), "ST Data file exists")
        self.assertTrue(os.path.getsize(datafile) > 1024, "ST Data file is not empty")
        self.assertTrue(os.path.exists(readsfile), "ST Data BED file exists")
        self.assertTrue(os.path.getsize(readsfile) > 1024, "ST Data BED file is not empty")
        # self.assertTrue(os.path.exists(statsfile), "Stats JSON file exists")

        # Run st_qa.py
        try:
            print("Running st_qa.py")
            check_call(["st_qa.py", datafile])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running st_qa.py failed \n")
        # TODO verify all output files are present
        self.assertTrue(os.path.exists("{}_qa_stats.txt".format(
            os.path.basename(datafile).split(".")[0])),
            "Output of st_qa.py file exists")

        # Run convertEnsemblToNames.py
        try:
            print("Running convertEnsemblToNames.py")
            check_call(["convertEnsemblToNames.py",
                        "--annotation", self.annotfile,
                        datafile])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running convertEnsemblToNames.py failed \n")
        self.assertTrue(os.path.exists("output.tsv"),
                        "Output of convertEnsemblToNames.py file exists")

        # Run filter_gene_type_matrix.py
        try:
            print("Running filter_gene_type_matrix.py")
            check_call(["filter_gene_type_matrix.py",
                        "--annotation", self.annotfile,
                        "--gene-types-keep", "protein_coding",
                        "--outfile", "output2.tsv",
                        datafile])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running filter_gene_type_matrix.py failed \n")
        self.assertTrue(os.path.exists("output2.tsv"),
                        "Output of filter_gene_type_matrix.py file exists")

        # TODO complete run tests for adjust_matrix_coordinates.py, merge_fastq.py and multi_qa.py
        try:
            print("Running adjust_matrix_coordinates.py")
            check_call(["adjust_matrix_coordinates.py", "--help"])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running adjust_matrix_coordinates.py failed \n")

        try:
            print("Running merge_fastq.py")
            check_call(["merge_fastq.py", "--help"])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running merge_fastq.py failed \n")

        try:
            print("Running multi_qa.py")
            check_call(["multi_qa.py", "--help"])
        except Exception as e:
            print(str(e))
            self.assertTrue(0, "Running multi_qa.py failed \n")


if __name__ == '__main__':
    unittest.main()
