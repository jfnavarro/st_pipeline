#!/usr/bin/env python
""" Stupid script to add a fake molecular barcodes to normal reads"""
import sys
from stpipeline.common.fastq_utils import *
from stpipeline.common.utils import *
import random 

def getRandomMolecularBarcode(base_sequence):
    molecular_barcode = ""
    for ele in base_sequence:
        if ele == 'W':
            molecular_barcode += random.choice('AT')
        elif ele == 'S':
            molecular_barcode += random.choice('GC')
        elif ele == 'N':
            molecular_barcode += random.choice('ATGC')
        elif ele == 'V':
            molecular_barcode += random.choice('ACG')
        else:
            molecular_barcode += ele    
    return molecular_barcode

if len(sys.argv) == 2:
    input = str(sys.argv[1])
else:
    print "Wrong number of parameters, it should be a fastq file"
    sys.exit(1)

if not input.endswith("fastq"):
    print "You must use fastq files"
    sys.exit(1)

output_name = "reformated_" + getCleanFileName(input)
outF = safeOpenFile(output_name,'w')
outF_writer = writefq(outF)
fastq_file = safeOpenFile(input, "rU")
base_sequence = "WSNNNWSNNNV"
length = len(base_sequence)
bacode_length = 27

for line in readfq(fastq_file):
    #line[0], line[1], line[2]
    normal_barcode = line[1][:bacode_length]
    molecular_barcode = getRandomMolecularBarcode(base_sequence)
    clean_read = line[1][bacode_length:]
    new_line = ( line[0], normal_barcode + molecular_barcode + clean_read, line[2] + "".join("B" for i in xrange(0, length)) )
    outF_writer.send(new_line)
     
print "File reformated " + output_name          
outF_writer.close()
outF.close()
fastq_file.close()