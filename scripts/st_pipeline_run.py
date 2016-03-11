#! /usr/bin/env python
""" This is a terminal based API to run the ST pipeline on a single node
"""

import sys
import argparse
from stpipeline.core.pipeline import *

def main(argv):
    
    pipeline = Pipeline()
    
    #create a parser
    parser = argparse.ArgumentParser(description=__doc__)
    parser = pipeline.createParameters(parser)
    
    #parse arguments
    options = parser.parse_args()
                          
    try:
        pipeline.load_parameters(options)
        pipeline.createLogger()
        pipeline.sanityCheck()
        pipeline.run()
    except Exception as e:
        print e
        sys.stderr.write("Error: Running the pipeline \n")
        sys.exit(-1)
        
if __name__ == "__main__":
    main(sys.argv[1:])
    

