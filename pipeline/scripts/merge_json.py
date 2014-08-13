#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""
""" Script for merging json ST-data files.

The major point is to merge (feature, gene, expression) triples from several
json files so that (f, g, e_1) + (f, g, e_2 ) = (f, g, e_1 + e_2).
"""
import argparse
from collections import defaultdict

from main.common.json_utils import json_iterator
from main.common.json_utils import write_json

def main(files, out):
    
    hits = defaultdict(int)

    for f in files:
        it = json_iterator(f)
        for doc in it:
            feature_gene = (doc['y'], doc['x'], doc['gene'], doc['barcode'])
            hits[feature_gene] += doc['hits']
    
    write_json(out,hits)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('json_files', nargs='+', type=str,
                        help='Space seperated list of json ST files to merge.')
    parser.add_argument('-o', '--out', type=str,
                        help='Name of merged output file')

    args = parser.parse_args()
    main(args.json_files, args.out)
