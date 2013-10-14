#!/usr/bin/env python
""" Super-ultra-extreme simple and basic stats generator for csv output files
for the ST pipeline in Amazon. 
Run it like python generateStats file.csv

Author : Jose Fernandez Navarro
"""

import os
import sys

if len(sys.argv) == 2:
    input = str(sys.argv[1])
else:
    print "Wrong number of parameters"
    sys.exit(1)

if not input.endswith("csv"):
    print "You must use csv files"
    sys.exit(1)

filehandler = open(input,"r")

unique_events = set()
unique_barcodes = set()
unique_genes = set()
count = 0

for line in filehandler.readlines():
    #['4', '13', 'Mycbp2', 'CGCTACCCTGATTCGACC', '1216']
    cols = line.split()
    if(len(cols) != 5):
        print "There was an error parsing the file"
        break
    unique_events.add( (str(cols[2]), str(cols[3]) ) )
    unique_barcodes.add(str(cols[3]))
    unique_genes.add(str(cols[2]))
    count += int(cols[4])
    
print "Stats for file " + str(input)
print "Unique Genes " + str(len(unique_genes))
print "Unique Barcodes " + str(len(unique_barcodes))
print "Unique Events " + str(len(unique_events))
print "Total Events " + str(count)
    