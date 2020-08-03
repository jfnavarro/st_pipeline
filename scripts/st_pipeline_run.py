#! /usr/bin/env python
""" 
ST Pipeline is a tool to process the Spatial Transcriptomics raw datasets.
The data is filtered, aligned to a genome, annotated to a reference,
demultiplexed by array coordinates and then aggregated by counts
that are not duplicates using the Unique Molecular Indentifiers. 
The output contains the counts matrix, a stats file, a log file
and a BED file with all the transcripts.

The ST Pipeline requires two fastq files, an IDs files (BARCODE, X, Y)
,the path to a STAR genome index,
the path to a annotation file in GTF format an a dataset name.

The ST Pipeline has many parameters, you can see a description of them
by typing : st_pipeline_run.py --help

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import sys
import argparse
from stpipeline.core.pipeline import Pipeline

def main(argv):
    
    # Create pipeline object
    pipeline = Pipeline()
    
    # Create a parser
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # Parse parameters, sanity check and run the pipeline                  
    try:
        parser = pipeline.createParameters(parser)
        # Parse arguments
        options = parser.parse_args()
        pipeline.load_parameters(options)
        sys.stdout.write("ST Pipeline, parameters loaded\n")
        pipeline.createLogger()
        sys.stdout.write("ST Pipeline, logger created\n")
        pipeline.sanityCheck()
        sys.stdout.write("ST Pipeline, sanity check passed. Starting the run.\n")
        pipeline.run()
        sys.stdout.write("ST Pipeline, run completed!\n")
    except Exception as e:
        sys.stderr.write("Error running the pipeline\n")
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)
    finally:
        pipeline.clean_filenames()
        
if __name__ == "__main__":
    main(sys.argv[1:])
    

