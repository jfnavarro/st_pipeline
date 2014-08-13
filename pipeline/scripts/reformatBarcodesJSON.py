#! /usr/bin/env python
@author Jose Fernandez <jose.fernandez.navarro@scilifelab.se>
""" Script to simple indent a JSON file containing barcodes.
It also has the option to add a floor value, so all the reads counts
above that value will be floored.
"""
import argparse
import os
import json

def json_iterator(json_file):
    """ Iterator over lines in an ST json file.
    """
    with open(json_file) as fh:
        for line in json.loads(fh.readlines()[0]):
            yield line

def main(json_file, output_file = "output.json", cutoff = None, noindent = False):

    adjustedList = list()
    it = json_iterator(json_file)
    indent = 2
    if noindent:
        indent = 0
        
    for doc in it:
        readsCount = int(doc['hits'])
        if cutoff is not None and cutoff >= 1 and int(doc['hits']) >= cutoff:
            doc['hits'] = cutoff
        adjustedList.append(doc)
        
    filehandler = open(output_file, "w")
    #write well formed json file
    json.dump(adjustedList, filehandler, indent, separators=(',', ': '))  
    filehandler.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_file", help="JSON ST-data file")
    parser.add_argument("--output", default="output.json")
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("--noindent", default=False, action='store_true')
    
    args = parser.parse_args()

    main(args.json_file, args.output, args.cutoff, args.noindent)
