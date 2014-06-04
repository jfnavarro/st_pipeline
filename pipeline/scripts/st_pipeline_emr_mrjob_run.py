#!/usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""

''' This is a wrapper that uses MrJob to run the ST pipeline using Map-Reduce on Amazon'''

"""CURRENTLY NOT BEING USED"""

import sys
from mrjob.job import MRJob
from mrjob.util import bash_wrap
import mrjob.protocol
import random
import tempfile
from main.core.pipeline import *
from main.common.json_utils import json_iterator

class EMRPipeline(MRJob):
        
    INPUT_PROTOCOL = mrjob.protocol.RawValueProtocol
    INTERNAL_PROTOCOL = mrjob.protocol.JSONProtocol
    OUTPUT_PROTOCOL = mrjob.protocol.JSONProtocol
    
    def mapper_init(self):
        """ this function is called when the streaming starts
        we create the varibles and the pipeline object that we initialize
        with the parameters given in the input
        """
        self.batch = []
        self.chunks = []  
        self.pipeline = Pipeline()  
        self.pipeline.allowed_missed = self.options.allowed_missed
        self.pipeline.allowed_kimera = self.options.allowed_kimer
        self.pipeline.min_length_trimming = self.options.min_length_qual_trimming
        self.pipeline.trimming_fw_bowtie = self.options.mapping_fw_trimming
        self.pipeline.trimming_rw_bowtie = self.options.mapping_rv_trimming
        self.pipeline.min_quality_trimming = self.options.min_quality_trimming 
        self.pipeline.clean = self.options.no_clean_up
        self.pipeline.s = self.options.start_id
        self.pipeline.l = self.options.length_id
        self.pipeline.e = self.options.error_id
        self.pipeline.threads = self.options.bowtie_threads
        self.pipeline.verbose = self.options.verbose
        self.pipeline.ids = os.path.abspath(self.options.ids)
        self.pipeline.ref_map = os.path.abspath(self.options.ref_map)
        self.pipeline.ref_annotation = os.path.abspath(self.options.ref_annotation)
        self.pipeline.expName = self.options.expName
        self.pipeline.htseq_mode = self.options.htseq_mode
        self.pipeline.htseq_no_ambiguous = self.options.htseq_no_ambiguous
        self.pipeline.qual64 = self.options.qual_64
        self.pipeline.discard_fw = self.options.discard_fw
        self.pipeline.discard_rv = self.options.discard_rv
        self.pipeline.discordant = self.options.bowtie2_discordant
        self.pipeline.contaminant_bt2_index = self.options.contaminant_bowtie2_index
        self.pipeline.path = self.options.bin_path
        if(self.options.log_file != ""):
            self.pipeline.logfile = os.path.abspath(self.options.log_file)
        self.pipeline.load_parameters()
        
    def run_pipeline(self,file1,file2,expName):
        """ this function runs the pipelin with the given
        fastq files and the name of the experiment
        """
        self.pipeline.Fastq_fw = os.path.abspath(file1)
        self.pipeline.Fastq_rv = os.path.abspath(file2)
        self.pipeline.expName = expName
        #TODO should wrap this onto try catch
        self.pipeline.sanityCheck()
        self.pipeline.run()
        #run should return if everything is okay
        
    def mapper_final(self):
        """ this function is called once all the streaming has been done
        in the input, it will iterate trough the chunks, create temp fastq files
        and call the pipeline on them, it will then parse the output and send
        the json formated features to the reducer
        """
        #some lines left behind?
        if len(self.batch) > 0:
            self.chunks.append(''.join(self.batch))
        
        #call pipeline and output streams
        for key,value in self.pipeline.run_pipeline(self.chunks):
            yield key,value
               
        self.increment_counter('SkippingTaskCounters','MapProcessedRecords',1)
        
    def configure_options(self):
        """ creates parameter parser
        """
        super(EMRPipeline, self).configure_options()
        self.add_passthrough_option('--chunks', type='int', default=10000, help='Number of reads per chunk to split')
        self.add_passthrough_option('--output-file', type='str', default="merged.json", help='Name of the output file')
        self.add_passthrough_option('--ids', type='str', default="",help='the name of the file containing the barcodes and the coordinates')
        self.add_passthrough_option('--ref-map', type='str', default="", help= "<path_to_bowtie2_indexes>] = reference genome name for the genome that you want to use to align the reads")
        self.add_passthrough_option('--ref-annotation', type='str', default="",help="select the reference(htseq requires a GTF file annotation file that you want to use to annotate")
        self.add_passthrough_option('--expName', type='str', default="", help="the name of the experiment (outfile name)")
        self.add_passthrough_option('--allowed-missed', type='int', default=6, help="number of allowed mismatches when mapping against the barcodes")
        self.add_passthrough_option('--allowed-kimer', type='int', default=7, help="kMer length when mapping against the barcodes")
        self.add_passthrough_option('--min-length-qual-trimming', type='int', default=28, help="minimum lenght of the sequence for mapping after trimming, shorter reads will be discarded")
        self.add_passthrough_option('--mapping-fw-trimming', type='int',  default=42, help="he number of bases to trim in the forward reads for the Mapping [24 + ID_LENGTH]")
        self.add_passthrough_option('--mapping-rv-trimming', type='int', default=5, help="he number of bases to trim in the reverse reads for the Mapping")
        self.add_passthrough_option('--length-id', type='int', default=18, help="length of ID, a.k.a. the length of the barcodes")
        self.add_passthrough_option('--contaminant-bowtie2-index', type='str', default="", help="<path_to_bowtie2_indexes>] => When provided, reads will be filtered against the specified bowtie2 index, non-mapping reads will be saved and demultiplexed")
        self.add_passthrough_option('--qual-64', action="store_true", default=False, help="use phred-64 quality instead of phred-33(default)")
        self.add_passthrough_option('--htseq-mode', type='str', default="intersection-nonempty", help="Mode of Annotation when using HTSeq. Modes = {union,intersection-nonempty(default),intersection-strict}")
        self.add_passthrough_option('--htseq-no-ambiguous', action="store_true", help="When using htseq discard reads annotating ambiguous genes")
        self.add_passthrough_option('--start-id', type='int', default=0, help="start position of BARCODE ID [0]")
        self.add_passthrough_option('--error-id', type='int', default=0, help="Id positional error [0]")
        self.add_passthrough_option('--no-clean-up', action="store_true", default=False, help="do not remove temporary files at the end (useful for debugging)")
        self.add_passthrough_option('--bowtie-threads', type='int', default=8, help="bowtie-threads")
        self.add_passthrough_option('--min-quality-trimming', type='int', default=20, help="minimum quality for trimming")
        self.add_passthrough_option('--discard-fw', action="store_true", default=False, help="discard fw reads that maps uniquely")
        self.add_passthrough_option('--discard-rv', action="store_true", default=False, help="discard rw reads that maps uniquely")
        self.add_passthrough_option('--bowtie2-discordant', action="store_true", default=False, help="discard non-discordant alignments when mapping")
        self.add_passthrough_option('--bin-path', type='str', default="", help="path to folder where binary executables are present")
        self.add_passthrough_option('--log-file', type='str', default="", help="name of the file that we want to use to store the logs(default output to screen)")
     
    def load_options(self, args):
        """ load the parameters from the terminal
        """
        super(EMRPipeline, self).load_options(args=args)
                  
    def mapper(self, _, value): 
        """make a list of chunks from input file (key is empty)
        tell mrjob to make the chunks in the hadoop config file
        """
        self.batch.append(value + '\n')
        if len(self.batch) == int(self.options.chunks):
            self.chunks.append(''.join(self.batch))
            self.batch = []
 
    def combiner(self, key, values):
        """ combiner is called after the mapper jobs for the similar key so
        we can save some time and work to the reducers
        """ 
        yield key, sum(value for value in values)
           
    def reducer(self, key, values):   
        """ the reducer functions sums up the hits of similar features (key)
        and sends them to the output
        """       
        yield key, sum(value for value in values)  
        
        #tell mrjob I do not want to run task counters
        self.increment_counter('SkippingTaskCounters','ReduceProcessedRecords',1)
        
    def steps(self):
        return [self.mr(mapper_init=self.mapper_init,
                        mapper=self.mapper,
                        mapper_final=self.mapper_final,
                        combiner=self.combiner, 
                        reducer=self.reducer,),]


if __name__ == '__main__':
    EMRPipeline.run()
    ##parse output and write it to a file
    #with open(EMRPipeline.options.output_file, 'a') as out_file:
        #for line in runner.stream_output():
            #key, value = MRJob.parse_output_line(line)
            #out_file.write('{}\n'.format(EMRPipeline.serialize(key,value)))
