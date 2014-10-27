import pytest
import os
import subprocess
import json

class TestClass:
	
	def test_st_pipeline_run_on_mouse_data_subset(self):
		'''Tests st_pipeline on a mouse data subset'''
		
		# Obtain paths and files.
		testdir = str(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
		
		infile1 = testdir + '/input/miseqF1/testdata_R1.fastq'
		infile2 = testdir + '/input/miseqF1/testdata_R2.fastq'
		genomedir = testdir + '/config/genomes/mouse_grcm38'
		contamdir = testdir + '/config/contaminant_genomes/R45S5_R5S1'
		annotfile = testdir + '/config/annotations/mouse_grcm38_chromosome19.gtf'
		chipfile = testdir + '/config/idfiles/130307_Design3_27mer.txt'
		expname = "mouse_subset_test"
		outdir = testdir + '/output'
		tmpdir = testdir + '/tmp'
		
		assert os.path.exists(infile1)
		assert os.path.exists(infile2)
		assert os.path.isdir(genomedir)
		assert os.path.isdir(contamdir)
		assert os.path.exists(annotfile)
		assert os.path.exists(chipfile)
		assert os.path.isdir(outdir)
		assert os.path.isdir(tmpdir)
		
		
		# Start st_pipeline_run
		#cmd = 'st_pipeline_run.py %s %s --ref-map %s/chromosome19 --contaminant-bowtie2-index %s/rnagenome --ref-annotation %s --ids %s --expName %s --output-folder %s --temp-folder %s --allowed-missed 5 --allowed-kimer 7 --min-length-qual-trimming 28 --mapping-fw-trimming 51 --mapping-rv-trimming 5 --length-id 27 --min-quality-trimming 20 --verbose --no-clean-up --bowtie-threads 3 --error-id 0 --start-id 0 --htseq-mode intersection-nonempty' % (infile1, infile2, genomedir, contamdir, annotfile, chipfile, expname, outdir, tmpdir)
		discard_output = open(os.devnull,"w")		
		args = ["st_pipeline_run.py", infile1, infile2, \
			"--ref-map", genomedir + "/chromosome19", \
			"--contaminant-bowtie2-index", contamdir + "/rnagenome", \
			"--ref-annotation", annotfile, \
			"--ids", chipfile, \
			"--expName", expname, \
			"--output-folder", outdir, \
			"--temp-folder", tmpdir, \
			"--allowed-missed", "5", \
			"--allowed-kimer", "7", \
			"--min-length-qual-trimming", "28", \
			"--mapping-fw-trimming", "51", \
			"--mapping-rv-trimming", "5", \
			"--length-id", "27", \
			"--min-quality-trimming", "20", \
			"--verbose", "--no-clean-up", \
			"--bowtie-threads", "3", \
			"--error-id", "0", \
			"--start-id", "0", \
			"--htseq-mode", "intersection-nonempty"]
		subprocess.check_call(args, stdout=discard_output, stderr=subprocess.PIPE)
		
		
		# Verify existence of output files.
		outcnt = os.listdir(outdir)
		assert (outcnt != [])
		barcodesfile = outdir + "/" + expname + "_barcodes.json"
		assert os.path.exists(barcodesfile)
		assert (os.path.getsize(barcodesfile) > 1024)
		readsfile = outdir + "/" + expname + "_reads.json"
		assert os.path.exists(readsfile)
		assert (os.path.getsize(readsfile) > 1024)
		
		
		# Validate that barcodes is valid JSON
		with open(barcodesfile) as json_file:
			json_data = json.load(json_file)
		
		
		# Verify existence of temp files.
		tmpcnt = os.listdir(tmpdir)
		assert (tmpcnt != [])
		
		
		# Clean up.
		for cnt in outcnt:
			os.remove(outdir + "/" + cnt)
		for cnt in tmpcnt:
			os.remove(tmpdir + "/" + cnt)
		assert (os.listdir(outdir) == [])
		assert (os.listdir(tmpdir) == [])
		