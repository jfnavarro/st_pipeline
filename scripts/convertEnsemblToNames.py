#! /usr/bin/env python
""" 
Script that parses a ST data file generated
with the pipeline in matrix (TSV) format where the genes are named
with ENSEMBL IDs and generates a new file
with the ENSEMBL IDs converted to gene names IDS.
For that the script also needs a tab delimited file like this:

ENSEMBL_ID GENE_NAME

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import sys
import json
from stpipeline.common.json_utils import json_iterator
from stpipeline.common.utils import fileOk

def main(st_data_file, names_map, output_file):

    if not fileOk(st_data_file) or not fileOk(names_map):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)
        
    # loads a map with the ensembl -> gene name
    genes_map = dict()
    with open(names_map, "r") as map_file:
        for line in map_file.readlines():
            tokens = line.split()
            assert(len(tokens) == 2)
            genes_map[tokens[0]] = tokens[1]
            
    # Iterates the genes IDs to get gene names
    st_data = pd.read_table(st_data_file, sep="\t", header=0, index_col=0)
    adjustedList = list()
    for gene in st_data.columns:
        try:
            gene = genes_map[gene]
        except KeyError:
            sys.stdout.write("Warning, {} was not found in the MAP file\n".format(gene))
        adjustedList.append(gene)
        
    # Update the table with the gene names
    st_data.columns = adjustedList
    
    # Write table to file
    st_data.to_csv(output_file, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("st_data_file", help="ST data file in TSV format")
    parser.add_argument("--output", default="output.tsv", 
                        help="Name of the output file, default output.tsv")
    parser.add_argument("--names-map", required=True,
                        help="File containing the map of ENSEMBL IDs to gene \
                        names as a two columns tab delimited file")
    args = parser.parse_args()
    main(args.st_data_file, args.names_map, args.output)

