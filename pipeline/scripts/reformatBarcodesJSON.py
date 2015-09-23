#! /usr/bin/env python
#@Author Jose Fernandez
""" 
Script to simply indents a JSON file containing barcodes.
It also has the option to add a floor value, so all the reads counts
above that value will be floored.
"""

import argparse
import json
from stpipeline.common.json_utils import json_iterator

def main(json_file, output_file="output.json", cutoff=None, noindent=False):

    adjustedList = list()
    it = json_iterator(json_file)
    indent = 2
    if noindent:
        indent = 0
        
    for doc in it:
        if cutoff is not None and cutoff >= 1 and int(doc['hits']) >= cutoff:
            doc['hits'] = cutoff
        adjustedList.append(doc)
        
    with open(output_file, "w") as filehandler:
        #write well formed json file
        json.dump(adjustedList, filehandler, indent, separators=(',', ': '))  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_file", help="JSON ST-data file")
    parser.add_argument("--output", default="output.json")
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("--noindent", default=False, action='store_true')
    
    args = parser.parse_args()

    main(args.json_file, args.output, args.cutoff, args.noindent)
