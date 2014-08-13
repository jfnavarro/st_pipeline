#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

import struct
import logging
import common
import monitor

from pydoop.pipes import Mapper, InputSplit
from pydoop.utils import jc_configure, jc_configure_int, jc_configure_bool

class Pipeline(object):
    
    def __init__(self):
        self.event_monitor = monitor.QuietMonitor()
        self.logger = monitor.QuietMonitor()   
        self.__batch = []
        ## I should put here all the parameters that will be set from the job runner
        self.qformat = "fastq-illumina"
        self.max_isize = 1000
        self.nthreads = 1
        self.trim_qual = 0
        self.reference = None
        self.annotation = None
        
        
        
        self.iterator = None ## would be smart to use an iterator to process the records
        #TODO rest of parameeters
    
    def load_pair_record(self, record):
        '''Append a tuple of the format (id, seq1, qual1, seq2, qual2) to this work batch.'''
        self.__batch.append(record)
        return len(self.__batch)

    def get_batch_size(self):
        return len(self.__batch)

    def clear_batch(self):
        self.__batch = []

    def write_batch_toFile(self):
        ''' write the current batch to a temp file 
        to be used as input for the pipeline'''
        return 
    
    def release_resources(self):
        
        self.__iterator = None
        
    def run_pipeline(self):
        
        #do checks and call the pipeline
        
        #if not self.reference:
        #    raise ValueError("You must set the reference path before calling run_alignment")
        #if not self.hit_visitor:
        #    raise ValueError("You must set the hit_visitor before calling run_alignment (else you'll lose the alignment results)")

        self.event_monitor.log_debug("Running pipeline in batch " + str(len(self.__batch)))
        self.event_monitor.count("reads processed", 2*len(self.__batch))

  
class mapper(Mapper):
    """
    Process sequences with STPipeline.

    @input-record: C{key} does not matter (standard LineRecordReader);
    C{value} is a tab-separated text line with 5 fields: ID, read_seq,
    read_qual, mate_seq, mate_qual.

    @output-record: protobuf-serialized mapped pairs (map-reduce job) or processsed records
    in json format (map-only job).

    @jobconf-param: C{mapred.reduce.tasks} number of Hadoop reduce tasks to launch.
    If the value of this property is set to 0, then the mapper will directly output
    the output json format.  If set to a value > 0 the mapper will output
    mappings in the protobuf serialized format for the rmdup reducer.

    @jobconf-param: C{mapred.batch.size}: how many
    sequences should be processed at a time by the pairing
    function. Status will be updated at each new batch: therefore,
    lowering this value can help avoid timeouts.

    @jobconf-param: C{mapred.create.symlink} must be set to 'yes'.

    @jobconf-param: C{stpipeline.log.level} logging level,
    specified as a logging module literal.

    @jobconf-param: C{mapred.cache.archives} distributed
    cache entry for the bowtie2 index archive, annotation files and so on. The entry
    is of the form HDFS_PATH#LINK_NAME..
    
    """
    
    SUPPORTED_FORMATS = "fastq-illumina"
    DEFAULT_FORMAT = "fastq-illumina"

    def __get_configuration(self, ctx):
        
        self.logger = logging.getLogger("STPipeline")
        jc = ctx.getJobConf()
        self.DeprecationMap = dict() ## deprecated variables
        jobconf = common.convert_job_conf(jc, self.DeprecationMap, self.logger)
        
        #use this to set parameters on the fly for the job comming from the runner
        jc_configure(self, jobconf, 'stpipeline.log.level', 'log_level', 'DEBUG')
        jc_configure(self, jobconf, "stpipeline.fastq-subformat", "format", self.DEFAULT_FORMAT)
        #jc_configure_int(self, jobconf, 'stpipeline.max.isize', 'max_isize', 1000)
        jc_configure_int(self, jobconf, 'stpipeline.batch.size', 'batch_size', 10000)
        #jc_configure_int(self, jobconf, 'stpipeline.min_hit_quality', 'min_hit_quality', 0)
        #jc_configure_bool(self, jobconf, 'stpipeline.remove_unmapped', 'remove_unmapped', False)
        #jc_configure_int(self, jobconf, 'stpipeline.nthreads', 'nthreads', 1)
        #jc_configure_int(self, jobconf, 'stpipeline.trim.qual', 'trim_qual', 0)

        #should do some checks here with the parameters/Configuration
        try:
            self.log_level = getattr(logging, self.log_level)
        except AttributeError:
            raise ValueError("Unsupported log level: %r" % self.log_level)
        
        if self.batch_size <= 0:
            raise ValueError("'stpipeline.batch.size' must be > 0, if specified [10000]")
          
        if jc.hasKey('mapred.reduce.tasks') and jc.getInt('mapred.reduce.tasks') > 0:
            self.__run_only = False
        else:
            self.__run_only = True

    def __is_last_record(self, k, v):
        return k + len(v) + 2 >= self.split_end

    def __init__(self, ctx):
        
        super(type(self), self).__init__(ctx)
        self.__get_configuration(ctx)
        logging.basicConfig(level=self.log_level)
        self.pipeline = Pipeline()
        
        ##assign variables and paths here (indexes and annotation can be shared on a virtual location)
        #self.pipeline.param = self.param
        
        # part of the code is a workaround for accumulating records, see #331
        isplit = InputSplit(ctx.getInputSplit())
        self.split_end = isplit.offset + isplit.length

    def map(self, ctx):
        # Accumulates reads in self.pairs, until batch size is reached or
        # until the input is finished. At that point it calls run_pipeline
        # and emits the output.
        k = struct.unpack(">q", ctx.getInputKey())[0]
        v = ctx.getInputValue()
        self.pipeline.load_pair_record(v.split("\t"))
        is_last_record = self.__is_last_record(k, v)
        if self.pipeline.get_batch_size() >= self.batch_size or is_last_record:
            self.pipeline.run_pipeline()
            self.pipeline.clear_batch()
        if is_last_record:
            self.pipeline.release_resources()
            
            
            
            
            
            
