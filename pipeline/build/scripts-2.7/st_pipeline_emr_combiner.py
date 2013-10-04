#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

''' This is the Amazon EMR wrapper for the combiner function to run the Map Reduce version of the ST pipeline'''

import sys
#import json
from collections import defaultdict
#from main.common.json_utils import json_iterator
#from main.common.json_utils import write_json

#TODO : write json file? 
#TODO : check the json features get merged properly

def main():
    
    hits = defaultdict(int)
    # input comes from STDIN
    for line in sys.stdin:
        # remove leading and trailing whitespace
        line = line.strip()
        word, count = line.split('\t', 1)
        try:
            count = int(count)
        except ValueError:
            continue
        #doc = json.loads(word)
        doc = eval(word)
        feature_gene = (doc['y'], doc['x'], doc['gene'], doc['barcode'])
        hits[feature_gene] += count
    
    for key,value in hits.iteritems():
        print "%s%s%d" % (key, '\t', value)

if __name__ == "__main__":
    main()