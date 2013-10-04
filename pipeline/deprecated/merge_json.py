#! /usr/bin/env python
# @Created by Jose Fernandez
""" Script for merging json ST-data files.

The major point is to merge (feature, gene, expression) triples from several
json files so that (f, g, e_1) + (f, g, e_2 ) = (f, g, e_1 + e_2).
"""
import argparse
from collections import defaultdict
import json

#import cjson

def serialize(feature_gene, hits):
    doc = {}
    doc['y'], doc['x'], doc['gene'], doc['barcode'] = feature_gene
    doc['hits'] = hits
    return json.loads(doc)
    #return cjson.encode(doc)

def write_json(out,hits):
    with open(out, 'w') as out_file:
        for k, v in hits.iteritems():
            out_file.write('{}\n'.format(serialize(k, v)))
            
def json_iterator(json_file):
    """ Iterator over lines in an ST json file.
    """
    with open(json_file) as fh:
        for line in json.loads(fh.readlines()[0]):
        #for line in cjson.decode(fh.readlines()[0]):
            yield line

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
