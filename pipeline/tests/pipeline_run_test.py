#! /usr/bin/env python
""" 
Unit-test for run-tests, it just tests that the pipeline runs and produces correct results
"""

import unittest
import urllib
import os
import tempfile
import time
import multiprocessing
from subprocess import check_call
from stpipeline.core.pipeline import *

class TestPipeline(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        # Obtain paths and files.
        testdir = str(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
        self.infile_fw = os.path.join(testdir, 'input/arrayjet_1002/testdata_R1.fastq')
        self.infile_rv = os.path.join(testdir, 'input/arrayjet_1002/testdata_R2.fastq')
        self.annotfile = os.path.join(testdir, 'config/annotations/Homo_sapiens.GRCh38.79_chr19.gtf')
        self.chipfile = os.path.join(testdir, 'config/idfiles/150204_arrayjet_1000L2_probes.txt')
        self.expname = "test"

        # Obtain temp dir
        self.tmpdir = tempfile.mkdtemp(prefix="st_pipeline_test_temp")
        print "ST Pipeline Test Temporary directory " + self.tmpdir
        self.outdir = tempfile.mkdtemp(prefix="st_pipeline_test_output")
        print "ST Pipeline Test Temporary output " + self.outdir
        self.error_file = os.path.join(self.tmpdir, 'error.log')
        self.logFile = tempfile.mktemp(prefix="st_pipeline_test_log")
        print "ST Pipeline Test Log file " + self.logFile

        # Create genome index dirs.
        self.genomedir = os.path.join(self.tmpdir, 'config/genomes/mouse_grcm38')
        os.makedirs(self.genomedir)

        # STAR contam dir
        self.contamdir = os.path.join(self.tmpdir, 'config/contaminant_genomes/R45S5_R5S1')
        os.makedirs(self.contamdir)

        # Download and unpack fasta files
        try:
            print "ST Pipeline Test Downloading genome files..."
            genomefasta = os.path.join(self.genomedir, "human_grcm38_chromosome19.fasta")
            genomefastagz = os.path.join(self.genomedir, "human_grcm38_chromosome19.fasta.gz")
            urllib.urlretrieve ("ftp://ftp.ensembl.org/pub/release-79/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz", 
                                genomefastagz)
            check_call(['gunzip', genomefastagz])
        except Exception as e:
            print e
            self.assertTrue(0, "Downloading genome files failed \n")

        # Make genome indexes
        try:
            print "ST Pipeline Test Creating genome index..."
            check_call(["STAR", "--runMode", "genomeGenerate",
                    "--runThreadN", str(multiprocessing.cpu_count() - 1),
                    "--genomeDir", self.genomedir,
                    "--genomeFastaFiles", genomefasta,
                    "--sjdbGTFfile", self.annotfile,
                    "--sjdbOverhang", "100"])

            print "ST Pipeline Test Creating contaminant genome index..."
            contamfasta = os.path.join(testdir, "config/contaminant_genomes/R45S5_R5S1/Rn45s_Rn5s.fasta")
            check_call(["STAR", "--runMode", "genomeGenerate",
                    "--runThreadN", str(multiprocessing.cpu_count() - 1),
                    "--genomeDir", self.contamdir,
                    "--genomeFastaFiles", contamfasta])
        except Exception as e:
            print e
            self.assertTrue(0, "Creating genome index failed \n")

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
               
        # Verify existence of input files
        assert(os.path.exists(self.infile_fw))
        assert(os.path.exists(self.infile_rv))
        assert(os.path.isdir(self.genomedir))
        assert(os.path.isdir(self.contamdir))
        assert(os.path.exists(self.annotfile))
        assert(os.path.exists(self.chipfile))
        assert(os.path.isdir(self.outdir))
        assert(os.path.isdir(self.tmpdir))

        # create a pipeline Instance
        self.pipeline = Pipeline()

        # init pipeline arguments
        self.pipeline.expName = self.expname
        self.pipeline.fastq_fw = self.infile_fw
        self.pipeline.fastq_rv = self.infile_rv
        self.pipeline.molecular_barcodes = True
        self.pipeline.mc_allowed_missmatches = 1
        self.pipeline.mc_start_position = 18
        self.pipeline.mc_end_position = 27
        self.pipeline.min_cluster_size = 2
        self.pipeline.trimming_fw = 31
        self.pipeline.keep_discarded_files = True
        self.pipeline.allowed_missed = 2
        self.pipeline.allowed_kimera = 6
        self.pipeline.min_length_trimming = 28
        self.pipeline.trimming_rv = 5
        self.pipeline.min_quality_trimming = 20
        self.pipeline.clean = False
        self.pipeline.barcode_start = 0
        self.pipeline.barcode_length = 18
        self.pipeline.threads = multiprocessing.cpu_count() - 1
        self.pipeline.verbose = True
        self.pipeline.ids = os.path.abspath(self.chipfile)
        self.pipeline.ref_map = os.path.abspath(self.genomedir)
        self.pipeline.ref_annotation = os.path.abspath(self.annotfile)
        self.pipeline.htseq_mode = "intersection-nonempty"
        self.pipeline.htseq_no_ambiguous = False
        self.pipeline.contaminant_index= os.path.abspath(self.contamdir)  
        self.pipeline.output_folder = os.path.abspath(self.outdir)
        self.pipeline.temp_folder = os.path.abspath(self.tmpdir)
        self.pipeline.logfile = self.logFile
        self.pipeline.remove_polyA_distance = 15
        self.pipeline.remove_polyT_distance = 15
        self.pipeline.remove_polyG_distance = 15
        self.pipeline.remove_polyC_distance = 15

    @classmethod
    def tearDownClass(self):
        return
        print "ST Pipeline Test Remove temporary output " + self.outdir
        for root, dirs, files in os.walk(self.outdir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))
        print "ST Pipeline Test Remove temporary directory " + self.tmpdir
        for root, dirs, files in os.walk(self.tmpdir, topdown=False):
            for name in files:
                os.remove(os.path.join(root, name))
            for name in dirs:
                os.rmdir(os.path.join(root, name))


    def validateOutputData(self, expName):
        # Verify existence of output files and temp files
        self.assertNotEqual(os.listdir(self.outdir), [], "Output folder is not empty")
        self.assertNotEqual(os.listdir(self.tmpdir), [], "Tmp folder is not empty")
        barcodesfile = os.path.join(self.outdir, "barcodes.json")
        readsfile = os.path.join(self.outdir, "reads.json")
        self.assertTrue(os.path.exists(barcodesfile), "Barcodes JSON file exists")
        self.assertTrue(os.path.getsize(barcodesfile) > 1024, "Barcordes JSON file is not empty")
        self.assertTrue(os.path.exists(readsfile), "Reads JSON file exists")
        self.assertTrue(os.path.getsize(readsfile) > 1024, "Reads JSON file is not empty")

    def test_normal_run(self):
        """
        Tests st_pipeline on a mouse data subset with normal fastq files
        """
        # Start the pipeline
        try:
            self.pipeline.createLogger()
            self.pipeline.sanityCheck()
            self.pipeline.run()
        except Exception as e:
            print e
            self.assertTrue(0, "Running Pipeline Test failed \n")

        self.validateOutputData(self.expname)

if __name__ == '__main__':
    unittest.main()