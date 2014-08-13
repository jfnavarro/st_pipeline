#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

import pydoop.hdfs as phdfs

import logging
import os
import random
import sys
import tempfile
import hadut as hadut
import config as config

""" this crases create the map/reduce jobs, it takes file inputs and create temp files for
    running the jobs using mapper/reducer classes """

class PipelineRun(object):

    DefaultReduceTasksPerNode = 6
    LogName = "STPipeline"
    DefaultLogLevel = 'DEBUG'
    ConfLogLevel = 'stpipeline.log.level'
    
    def __init__(self):
        
        self.parser = config.Config() ## load parameters
        
        # set default properties
        self.properties = {
            self.ConfLogLevel: self.DefaultLogLevel,
            'hadoop.pipes.java.recordreader': 'true',
            'hadoop.pipes.java.recordwriter': 'true',
            'mapred.create.symlink': 'yes',
            'mapred.compress.map.output': 'true',
            'bl.libhdfs.opts': '-Xmx48m'
        }
        self.hdfs = None
        self.options = None

    def parse_cmd_line(self, args):
        
        self.options, self.left_over_args = self.parser.load_config_and_cmd_line(args)
        
        # set the job name. Do it here so the user can override it
        self.properties['mapred.job.name'] = 'stpipeline_%s' % self.options.output

        # now collect the property values specified in the options and
        # copy them to properties
        for k,v in self.options.properties.iteritems():
            self.properties[k] = v

        #set properties for the pipeline here... (parameters read in Config object)
        # we need to make sure the files are present
        # reference
        #self.properties['mapred.cache.archives'] = '%s#reference' % self.options.reference
        #self.properties['stpipelien.batchsize] = ..
        
        # create a logger
        logging.basicConfig()
        self.logger = logging.getLogger(self.__class__.LogName)
        # temporarily set to a high logging level in case we have to print warnings
        # regarding deprecated properties
        self.logger.setLevel(logging.DEBUG)

        if self.options.num_reducers and self.options.num_reducers > 0:
            self.logger.warning("Number of reduce tasks must be 0 when doing --align-only.")
            self.logger.warning("Ignoring request for %d reduce tasks", self.options.num_reducers)
        elif self.options.num_reducers:
            n_red_tasks = self.options.num_reducers
        else:
            n_red_tasks = PipelineRun.DefaultReduceTasksPerNode * hadut.num_nodes()

        self.properties['mapred.reduce.tasks'] = n_red_tasks

    def __write_pipes_script(self, fd):
        
        ld_path = ":".join( filter(lambda x:x, [os.environ.get('LD_LIBRARY_PATH', None)]) )
        pypath = os.environ.get('PYTHONPATH', '')
        self.logger.debug("LD_LIBRARY_PATH for tasks: %s", ld_path)
        self.logger.debug("PYTHONPATH for tasks: %s", pypath)

        fd.write("#!/bin/bash\n")
        fd.write('""":"\n')
        # should we set HOME to ~?  Hadoop by default sets $HOME to /homes, unless the
        # cluster administrator sets mapreduce.admin.user.home.dir.  This kills local installations
        #fd.write('[ -d "${HOME}" ] || export HOME="$(echo ~)"\n')
        # which causes python not to add installations under ~/.local/ to the PYTHONPATH
        fd.write('export LD_LIBRARY_PATH="%s" # Pipeline dir + LD_LIBRARY_PATH copied from the env where you ran %s\n' % (ld_path, sys.argv[0]))
        fd.write('export PYTHONPATH="%s"\n' % pypath)
        if self.logger.isEnabledFor(logging.DEBUG):
            fd.write('env >&2\n') # write the environment to the stderr log
            fd.write('echo >&2; cat $0 >&2\n') # write the script to the stderr log
        fd.write('exec "%s" -u "$0" "$@"\n' % sys.executable)
        fd.write('":"""\n')
        script = """
import sys
try:
    from pydoop.pipes import runTask, Factory
    from mapper import mapper
    from reducer import reducer
    runTask(Factory(mapper, reducer))
except ImportError as e:
    sys.stderr.write(str(e) + "\\n")
    sys.stderr.write("Can't import a required module\\n")
    sys.stderr.write("Did you install the Pipeline to a system path on all the nodes?\\n")
    sys.stderr.write("If you installed to a non-system path (e.g. your home directory)\\n")
    sys.stderr.write("you'll have to set PYTHONPATH to point to it.\\n")
    sys.stderr.write("Current Python library paths:\\n")
    sys.stderr.write("  sys.path: %s:\\n" % ':'.join(map(str, sys.path)))
"""
        fd.write(script)
        
    def run(self):
        
        if self.options is None:
            raise RuntimeError("You must call parse_cmd_line before run")

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug("Running STPipeline")
            self.logger.debug("Properties:\n%s", "\n".join( sorted([ "%s = %s" % (str(k), str(v)) for k,v in self.properties.iteritems() ]) ))
        self.logger.info("Input: %s; Output: %s; reference: %s", self.options.input, self.options.output, self.options.reference)

        try:
            self.hdfs = phdfs.hdfs('default', 0)
            self.__validate()

            self.remote_bin_name = tempfile.mktemp(prefix='stpipeline_bin.', suffix=str(random.random()), dir='')
            try:
                with self.hdfs.open_file(self.remote_bin_name, 'w') as script:
                    self.__write_pipes_script(script)

                full_name = self.hdfs.get_path_info(self.remote_bin_name)['name']

                return hadut.run_pipes(full_name, self.options.input, self.options.output,
                    properties=self.properties, args_list=self.left_over_args)
            finally:
                try:
                    self.hdfs.delete(self.remote_bin_name) # delete the temporary pipes script from HDFS
                    self.logger.debug("pipes script %s deleted", self.remote_bin_name)
                except:
                    self.logger.error("Error deleting the temporary pipes script %s from HDFS", self.remote_bin_name)
                    ## don't re-raise the exception.  We're on our way out
        finally:
            if self.hdfs:
                tmp = self.hdfs
                self.hdfs = None
                tmp.close()
                self.logger.debug("HDFS closed")
                
    def __validate(self):
        
        if self.properties['mapred.reduce.tasks'] == 0:
            self.logger.warning("Running in single core mode (no rmdup).")

        # validate conditions and files
        
        if phdfs.path.exists(self.options.output):
            raise RuntimeError("Output directory %s already exists. Please delete it or specify a different output directory." % self.options.output)
        #if not phdfs.path.exists(self.options.reference):
            #raise RuntimeError("Can't read reference archive %s" % self.options.reference)
        
        
        
