#!/usr/bin/env python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

# Input should be the output directory of the PairReadsQSeq (PRQ) application.
#
# All paths are HDFS paths, and may be relative to your HDFS home
# (/user/<your username>);

# hadoop is expected to be either in $HADOOP_HOME/bin or in the PATH;
# if you use a non-standard Hadoop configuration directory, set
# HADOOP_CONF_DIR accordingly.

import logging
import sys
import hadut as hadut
import run
import config

def main(argv=None):
    
    print >>sys.stderr, "Using hadoop executable %s" % hadut.hadoop
    retcode = 0
    job = run.PipelineRun()

    try:
        job.parse_cmd_line(argv)
        retcode = job.run()
    except config.ConfigParser.Error as e:
        logger = logging.getLogger(run.PipelineRun.LogName)
        logger.critical("Error in STPipeline run configuration")
        logger.critical(">>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        logger.critical(e)
        logger.critical(">>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        retcode = 1

    if retcode != 0:
        print >>sys.stderr, "Error running STPipeline"
        
    return retcode

if __name__ == "__main__":
    main(sys.argv[1:])
