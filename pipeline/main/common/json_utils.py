#!/usr/bin/env python
# @Created by Jose Fernandez
""" Module for parsing and dealing with ST data formats.
"""
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

def save_json(data, json_file):
    """ Save data in ST json format.
    """
    data = []
    with open(json_file, "w") as fh:
        for datum in data:
            data.append(datum)
        fh.write(json.dumps(data))

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
        for line in json.loads(fh.readlines()[0]):
        #for line in cjson.decode(fh.readlines()[0]): ##faster but annoying dependency
            data.append(json.loads(line))
            #data.append(cjson.decode(line))
    return data

