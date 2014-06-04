#! /usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""

""" 
    Scripts that parses a tab delimited file with the format (Name,chromosome,gene,barcode,x,y,Qul,Read)
    into two json files one containing the features and another one containing the reads
"""
import sys
import os
import json
import argparse

def getNameToIdMap(NameToIdFile):
    ''' Parse the output from the barcode mapper, 
        VERY CAREFULL, the order of columns has to be like this  = Name,chromosome,gene,barcode,x,y,Qul,Read 
        TODO : refactor this, no need to rely on columns 
    '''
    nameToId = dict()
    inF = open(NameToIdFile,'r')
    
    for line in inF: 
        cols = line.split()
        if(len(cols) == 8):
            gene = str(cols[2].replace("Gene:",""))
            clear_name = str(cols[0].replace("@",""))
            chromosome = str(cols[1].replace("Chr:",""))
            barcode = str(cols[3])
            x = int(cols[4])
            y = int(cols[5])
            qual = str(cols[6])
            seq = str(cols[7])
            nameToId[clear_name] = [barcode,x,y,qual,seq,gene]
        else:
            sys.stderr.write("Error: parse Name to Id file, wrong number of columns " + str(len(cols)) + " \n")
            sys.exit(1)
            
    inF.close()
    return nameToId
        
def getIdMap(nameToId):
    """create a map of barcode -> gene -> transcripts 
    """
    idMap = dict()
    for name in nameToId.keys(): ## go trough the reads that have barcodes mapped  

        gene = nameToId[name][5] ## get the gene name
        Id = nameToId[name][0] ## get the barcode
        
        if Id in idMap: # if barcode is already added
            if gene in idMap[Id]: # if the gene is the same, increase hits and append reads and qualities and names
                idMap[Id][gene][0] += 1
                idMap[Id][gene][4].append(nameToId[name][3])
                idMap[Id][gene][5].append(nameToId[name][4])
                idMap[Id][gene][6].append(name)
            else: # create new record for the gene
                idMap[Id][gene] = [1,nameToId[name][0],nameToId[name][1],nameToId[name][2], 
                                  [nameToId[name][3]], [nameToId[name][4]], [name] ]
        else: # create new record for barcode and gene
            idMap[Id] = dict()
            idMap[Id][gene] = [1,nameToId[name][0],nameToId[name][1],nameToId[name][2], 
                              [nameToId[name][3]], [nameToId[name][4]], [name] ]
            
    return idMap

def main(NameToIdFile, dbName, output_folder):
    
    if  NameToIdFile is None or dbName is None or not os.path.isfile(NameToIdFile):
        sys.stderr.write("Error, one of the input file/s not present\n")
        sys.exit()

    if output_folder is None or not os.path.isdir(output_folder):
        output_folder = "."
        
    nameToId = getNameToIdMap(NameToIdFile)
    idMap = getIdMap(nameToId)
    
    total_record = 0
    json_barcodes = list()
    json_reads = list()
    unique_genes = set()
    unique_barcodes = set()
    total_barcodes = 0

    for Id in idMap.keys():  # for each barcode
        for g in idMap[Id].keys():  # for each gene
            x = int(idMap[Id][g][2])
            y = int(idMap[Id][g][3])
            hits = int(idMap[Id][g][0])
            #add json line with the barcode information
            json_barcodes.append({"barcode":Id,"gene":g,"x":x,"y":y,"hits":hits})
            #get the reads that mapped to the barcode
            for qula, read, name in zip(idMap[Id][g][4], idMap[Id][g][5], idMap[Id][g][6]):
                #add json line with the raw reads information
                json_reads.append({"name":str(name),"read":str(read),"quality":str(qula),"barcode":Id,"gene":g})
                
            unique_genes.add(str(g))
            unique_barcodes.add(str(Id))
            total_record += 1
            total_barcodes += int(hits)
    
    if(total_record == 0):
        sys.stderr.write("Error: the number of transcripts present is 0\n")
        sys.exit(1)
    
    print "Number of Transcripts with Barcode present : " + str(total_barcodes) 
    print "Number of unique events present : " + str(total_record) 
    print "Number of unique Barcodes present : " + str(len(unique_barcodes))
    print "Number of unique Genes present : " + str(len(unique_genes))
    
    filename = dbName + "_barcodes.json"
    filenameReads = dbName + "_reads.json"
    #jose (use generic folder joining)
    filehandler = open(os.path.join(output_folder,filename), "w")
    filehandlerReads = open(os.path.join(output_folder,filenameReads), "w")
    #write well formed json file
    json.dump(json_barcodes,filehandler,separators=(',', ': '))    
    json.dump(json_reads,filehandlerReads,separators=(',', ': '))    
    filehandler.close()
    filehandlerReads.close()    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i','--input', type=str,
                        help='Input files containing the records (Name,chromosome,gene,barcode,x,y,Qul,Read)')
    parser.add_argument('-o', '--out', type=str,
                        help='Path of the output folder (default is .)')
    parser.add_argument('-n', '--name', type=str,
                        help='Name of the output files')

    args = parser.parse_args()
    main(args.input, args.name, args.out)
                                    
