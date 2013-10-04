#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

from pydoop.pipes import Reducer
from pydoop.utils import jc_configure, jc_configure_bool
import logging
import common
from hadoopMonitor import HadoopEventMonitor

class HitProcessorChainLink(object):
    
    def __init__(self, next_link = None):
        self.next_link = next_link

    def set_next(self, link):
        self.next_link = link
        return link

    def process(self, pair):
        """
        Identity action. Passes to the next link, if any.
        """
        if self.next_link:
            self.next_link.process(pair)
            
class EmitLink(HitProcessorChainLink):
    
    def __init__(self, context, event_monitor, next_link = None):
        super(type(self), self).__init__(next_link)
        self.ctx = context
        self.event_monitor = event_monitor

    def process(self, pair):
        for hit in pair:
            if hit:
                k, v = hit.split("\t", 1)
                self.ctx.emit(str(k), str(v))
                self.event_monitor.count("emitted records", 1)

        super(type(self), self).process(pair) # forward pair to next element in chain
        
class reducer(Reducer):

    COUNTER_CLASS = "STPipeline"
    DeprecationMap = dict()
    
    def __init__(self, ctx):
        
        super(reducer, self).__init__(ctx)
        logger = logging.getLogger("STPipeline")
        
        jc = ctx.getJobConf()  
        jobconf = common.convert_job_conf(jc, self.DeprecationMap, logger)
        jc_configure(self, jobconf, 'stpipeline.log.level', 'log_level', 'INFO')
        #jc_configure_bool(self, jobconf, 'seal.seqal.discard_duplicates', 'discard_duplicates', False)
        
        logging.basicConfig(level=self.log_level)
        self.event_monitor = HadoopEventMonitor(self.COUNTER_CLASS, logging.getLogger("reducer"), ctx)
        self.__output_sink = EmitLink(ctx, self.event_monitor)
        
    def reduce(self, ctx):
        # create the "workspace"
        self.__records = []

        # gather input
        key_values = ctx.getInputKey().split(':')
        self.event_monitor.log_debug("Reducer key values : " + key_values)
        
        # load mappings
        while ctx.nextValue():
            value = ctx.getInputValue()
            self.event_monitor.log_debug("Reducing..." + value)
            self.__output_sink.process( (value, None) )
            self.__records.append(value)
            
        # clean-up the workspace
        self.__records = None

