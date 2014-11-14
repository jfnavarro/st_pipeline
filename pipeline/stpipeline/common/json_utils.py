#!/usr/bin/env python
""" Module for parsing and dealing with ST data formats.
"""

import json

def serialize(feature_gene, hits):
    doc = {}
    doc['y'], doc['x'], doc['gene'], doc['barcode'] = feature_gene
    doc['hits'] = hits
    return json.loads(doc)

def write_json(out,hits):
    with open(out, 'w') as out_file:
        for k, v in hits.iteritems():
            out_file.write('{}\n'.format(serialize(k, v)))
            
def json_iterator(json_file):
    """ Iterator over lines in an ST json file.
    """
    with open(json_file) as fh:
        for line in json.load(fh):
            yield line

def save_json(data, json_file):
    """ Save data in ST json format.
        data must be a list of dict types with ST barcodes
    """
    data = []
    with open(json_file, "w") as fh:
        for datum in data:
            data.append(datum)
        fh.write(json.dumps(data, indent=2, separators=(',', ': ')))

def load_id_map(id_file):
    """ Load a ids file in to a barcode -> coordinate dictionary.
    """
    id_map = {}
    with open(id_file, "r") as fh:
        for line in fh:
            bc, x, y = line.split("\t")
            id_map[bc] = (int(x), int(y))
    return id_map

def load_json(json_file):
    """ Load a json file with e.g. expression data.
    """
    data = []
    with open(json_file) as fh:
        for line in json.load(fh):
            data.append(json.loads(line))
    return data

