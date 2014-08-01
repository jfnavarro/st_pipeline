#!/usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" This is the main API for the ST pipeline, it needs a bunch of files and parameters in order
to run the jobs, input files are fastq, output files are json. It logs status into a file.
"""
import sys
from main.common.utils import *
from main.core.mapping import *
from main.core.annotation import *
from main.common.json_utils import *
import os
from glob import glob
import logging
import subprocess
import random
import tempfile

class Pipeline():
    
    LogName = "STPipeline"
    DefaultLogLevel = 'DEBUG'
    
    def __init__(self):
        
        self.allowed_missed = 3
        self.allowed_kimera = 6
        self.min_length_trimming = 28
        self.trimming_fw_bowtie = 42
        self.trimming_rw_bowtie = 5 
        self.min_quality_trimming = 20 
        self.clean = True
        self.s = 0
        self.l = 18
        self.e = 0
        self.threads = 8
        self.verbose = False
        self.ids = None
        self.ref_map = None
        self.ref_annotation = None
        self.expName = None
        self.htseq_mode = "intersection-nonempty"
        self.htseq_no_ambiguous = False
        self.qual64 = False
        self.discard_fw = False
        self.discard_rv = False
        self.discordant = False
        self.contaminant_bt2_index = None
        self.Fastq_fw = None
        self.Fastq_rv = None
        self.path = None
        self.logger = None
        self.logfile = None
        self.output_folder = None
        
    def sanityCheck(self):
        
        conds = {"FW": fileOk(self.Fastq_fw), "RV": fileOk(self.Fastq_rv), 
                 "ids": fileOk(self.ids), "ref": fileOk(self.ref_annotation), 
                 "map": self.ref_map is not None, "Exp Name":  self.expName is not None}
        
        conds["htseq_gtf"] = self.ref_annotation.endswith("gtf")
        conds["htseq_mode"] = self.htseq_mode in ["union","intersection-nonempty","intersection-strict"]

        if not all(conds.values()):
            self.logger.error("Error: required file/s and or parameters not found or incorrect parameters :" + str(conds))
            raise RuntimeError("Error: required file/s and or parameters not found or incorrect parameters :" + str(conds))

        #test the presence of the scripts :
        required_scripts = set(['findIndexes','bowtie2'])

        unavailable_scripts = set()
        for script in required_scripts:
            if which(script) is None: 
                unavailable_scripts.add(script)
         
        if len(unavailable_scripts) == 0:
            self.logger.info("All tools present..starting the analysis")
        else:
            self.logger.error("Error, these programs not found:\t".join(unavailable_scripts))
            raise RuntimeError("Error, programs not found")
            
    def load_parameters(self):
      
        #TODO load the parameters here instead of forcing users to do so from outside
        
        # create a logger
        if self.logfile is not None and fileOk(self.logfile):
            logging.basicConfig(filename=self.logfile ,level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.DEBUG)
        self.logger = logging.getLogger(self.__class__.LogName)
        
        #load the given path into the system PATH
        if self.path is not None and os.path.isdir(self.path): 
            os.environ["PATH"] += os.pathsep + self.path
            
        #show parameters information and write them to stats
        parameters = "Parameters : m(" + str(self.allowed_missed) + ")" + \
                     "k(" + str(self.allowed_kimera) + ")" + "f(" + str(self.min_length_trimming) + ")" + \
                     "e(" + str(self.e) + ")" + "s(" + str(self.s) + ")" + "l(" + str(self.l) + ")" + \
                     "F(" + str(self.trimming_fw_bowtie) + ")" + "R(" + str(self.trimming_rw_bowtie) + ")"
    
        self.logger.info("Experiment : " + str(self.expName))
        self.logger.info("Forward reads file : " + str(self.Fastq_fw))
        self.logger.info("Reverse reads file : " + str(self.Fastq_rv))
        self.logger.info("Ids file : " + str(self.ids))
        self.logger.info("Reference mapping file : " + str(self.ref_map))
        self.logger.info("Reference annotation file : " + str(self.ref_annotation))
        if(self.contaminant_bt2_index):
            self.logger.info("Using bowtie2 contamination filter with " + str(self.contaminant_bt2_index))
        self.logger.info("Nodes : " + str(self.threads))
        self.logger.info(parameters)
        self.logger.info("Mapper : bowtie2")
        self.logger.info("Annotation Tool :  HTSeq")
  
    def run_pipeline(self,chunks):
        """ this function is called for Map Reduce jobs, when we want to run the pipeline,
        once all the streaming has been done
        in the input, it will iterate trough the chunks, create temp fastq files
        and call the pipeline on them, it will then parse the output and send
        the json formated features to the reducer
        """
        #TODO refactor and optimize this
        #TODO do mapping using gene as KEY
        for val in chunks:
            temp_name = tempfile.mktemp(prefix='stpipeline_temp_', suffix=str(random.random()), dir='') 
            new_filename = temp_name + "_1.fastq"
            new_filename2 = temp_name + "_2.fastq"
            outF = safeOpenFile(new_filename, 'w')
            outF_writer = writefq(outF)
            outF2 = safeOpenFile(new_filename2, 'w')
            outF_writer2 = writefq(outF2)
            for line in val.split("\n"):
                cols = line.replace("\t", " ").split(" ")
                if(len(cols) != 5):
                    continue
                header = cols[0]
                seq1 = cols[1]
                qual1 = cols[2]
                seq2 = cols[3]
                qual2 = cols[4]
                outF_writer.send((header, seq1, qual1))
                outF_writer2.send((header, seq2, qual2))
            outF.close()
            outF2.close()
            outF_writer.close()
            outF_writer2.close()
            #run pipeline with new files and name
            self.Fastq_fw = os.path.abspath(new_filename)
            self.Fastq_rv = os.path.abspath(new_filename2)
            self.expName = temp_name
            self.sanityCheck()
            #TODO wrap this into try catch
            self.run()
            # now we parse the output of pipeline (json) 
            outputPipeline = temp_name + "_barcodes.json"  # #this is very ugly (output file should be a parameter in the pipeline)
            outputPipelineReads = temp_name + "_reads.json"  # #this is very ugly (output file should be a parameter in the pipeline)
            it = json_iterator(outputPipeline)
            # send features to the reducer jobs (reads are ignored for now)
            for doc in it:
                feature_gene = (doc['y'], doc['x'], doc['gene'], doc['barcode'])
                hits = doc['hits']
                doc_json_formated = {}
                doc_json_formated['y'], doc_json_formated['x'], doc_json_formated['gene'], doc_json_formated['barcode'] = feature_gene
                #TODO would be nice to be able to send json objects trough
                yield doc_json_formated,hits
                
            # remove temp files
            safeRemove(new_filename)
            safeRemove(new_filename2)
            safeRemove(outputPipeline)
            safeRemove(outputPipelineReads)   
              
    def run(self):
        globaltime = TimeStamper()
        #starting time
        start_exe_time = globaltime.getTimestamp()
        self.logger.info("Starting the pipeline : " + str(start_exe_time))
        
        # add BC and PolyT from FW reads to the RW reads and apply quality filter
        Fastq_fw_trimmed, Fastq_rv_trimmed = reformatRawReads(self.Fastq_fw, self.Fastq_rv, 
                                                              self.trimming_fw_bowtie,
                                                              self.trimming_rw_bowtie, self.min_quality_trimming,
                                                              self.min_length_trimming, self.qual64)
        # First, do mapping against genome of both strands
        sam_mapped = bowtie2Map(Fastq_fw_trimmed, Fastq_rv_trimmed, self.ref_map, 
                                self.trimming_fw_bowtie, self.threads, self.qual64, self.discordant)
        
        ## filter unmapped and discordant reads
        sam_filtered = filterUnmapped(sam_mapped, self.discard_fw, self.discard_rv)
        if self.clean: safeRemove(sam_mapped)  
        
        ##annotate using htseq count
        annotatedFile = annotateReadsWithHTSeq(sam_filtered, self.ref_annotation, self.htseq_mode)
        if self.clean: safeRemove(sam_filtered)
    
        # get raw reads and quality from the forward and reverse reads
        withTr = getAnnotatedReadsFastq(annotatedFile, Fastq_fw_trimmed, 
                                        Fastq_rv_trimmed, self.htseq_no_ambiguous)
        if self.clean: safeRemove(annotatedFile)
        
        # Filter out contaminated reads with Bowtie2
        if self.contaminant_bt2_index: 
            oldWithTr = withTr
            withTr, contaminated_sam = bowtie2_contamination_map(withTr, self.contaminant_bt2_index,
                                                                 trim=self.trimming_fw_bowtie,
                                                                 cores=self.threads, qual64=self.qual64)
            if self.clean: safeRemove(contaminated_sam)
            if self.clean: safeRemove(oldWithTr)
    
        if self.clean: safeRemove(Fastq_fw_trimmed)
        if self.clean: safeRemove(Fastq_rv_trimmed)
        
        # Map against the barcodes
        mapFile = getTrToIdMap(withTr, self.ids, self.allowed_missed, self.allowed_kimera, 
                               self.s, self.l, self.e)
        if self.clean: safeRemove(withTr)
    
        # create json files with the results
        self.createDataset(mapFile, self.expName)
        if self.clean: safeRemove(mapFile)
        
        finish_exe_time = globaltime.getTimestamp()
        total_exe_time = finish_exe_time - start_exe_time
        self.logger.info("Total Execution Time : " + str(total_exe_time))

    def createDataset(self, mapFile, dbName):
        ''' parse annotated and mapped reads with the reads that contain barcodes to
            create json files with the barcodes and cordinates and json file with the raw reads
            and some useful stats and plots
        '''
        self.logger.info("Start Creating dataset")
        args = ['createDataset.py', '--input', str(mapFile), '--name', str(dbName)]
        if self.output_folder is not None: args += ['--output', str(self.output_folder)]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
        ##TODO should check for errors
        procOut = stdout.split("\n")
        self.logger.info('Create dataset stats :')
        for line in procOut: 
            self.logger.info(str(line))
        self.logger.info("Finish Creating dataset")
        return   