#!/usr/bin/env python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

''' This is the Amazon EMR wrapper for the reducer function to run the Map Reduce version of the ST pipeline'''

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

        doc = eval(word)
        feature_gene = (doc['x'], doc['y'], doc['gene'], doc['barcode'])
        hits[feature_gene] += count
    
    for key,value in hits.iteritems():
        print "%s%s%s%s%s%s%s%s%d" % (key[0],'\t', key[1], '\t', key[2], '\t', key[3], '\t', value)

if __name__ == "__main__":
    main()
