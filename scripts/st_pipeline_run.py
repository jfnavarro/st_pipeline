#! /usr/bin/env python
""" 
ST Pipeline is a tool to process the Spatial Transcriptomics raw data.
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers. 
The output contains the counts table, a stats file, a log file
and a BED file with all the transcripts.\n\n

The ST Pipeline requires two fastq files, an IDs files (BARCODE, X, Y)
,the path to a STAR genome index,
the path to a annotation file in GTF format an a dataset name.\n\n

The ST Pipeline has many parameters, you can see a description of them
by typing : st_pipeline_run.py --help\n\n

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import sys
from stpipeline.core.pipeline import *

def main(argv):
    
    # Create pipeline object
    pipeline = Pipeline()
    
    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__)
        
    # Parse parameters, sanity check and run the pipeline                  
    try:
        parser = pipeline.createParameters(parser)
        # Parse arguments
        options = parser.parse_args()
        pipeline.load_parameters(options)
        sys.stdout.write("ST Pipeline, parameters loaded")
        pipeline.createLogger()
        sys.stdout.write("ST Pipeline, logger created")
        pipeline.sanityCheck()
        sys.stdout.write("ST Pipeline, sanity check passed. Starting the run.")
        pipeline.run()
        sys.stdout.write("ST Pipeline, run completed!")
    except Exception as e:
        sys.stderr.write("Error running the pipeline\n")
        sys.stderr.write(str(e))
        sys.exit(1)
        
if __name__ == "__main__":
    main(sys.argv[1:])
    

