#! /usr/bin/env python
""" 
Script that parses a ST features file JSON output
with ENSEMBL IDs as gene names and convert
the gene names to the real gene names using
a MAP file as input. The map file
must have two columns:

ENSEMBL_ID GENE_NAME

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import sys
import json
from stpipeline.common.json_utils import json_iterator
from stpipeline.common.utils import fileOk

def main(json_file, names_map, output_file):

    if not fileOk(json_file) or not fileOk(names_map):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)
        
    # loads a map with the ensembl -> gene name
    genes_map = dict()
    with open(names_map, "r") as map_file:
        for line in map_file.readlines():
            tokens = line.split()
            assert(len(tokens) == 2)
            genes_map[tokens[0]] = tokens[1]
            
    # iterates the JSON file to change the name
    adjustedList = list()
    it = json_iterator(json_file)
    for doc in it:
        try:
            doc['gene'] = genes_map[doc['gene']]
        except KeyError:
            sys.stdout.write("Warning, {} was not found in the MAP file\n".format(doc['gene']))
        adjustedList.append(doc)
        
    with open(output_file, "w") as filehandler:
        # write well formed json file
        json.dump(adjustedList, filehandler, indent=2, separators=(',', ': '))  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json_file", help="JSON ST-data file")
    parser.add_argument("--output", default="output.json", 
                        help="Name of the output file, default output.json")
    parser.add_argument("--names-map", required=True,
                        help="File containing the map of ensembl ID to gene \
                        name as a two columns tab delimited file")
    args = parser.parse_args()
    main(args.json_file, args.names_map, args.output)

