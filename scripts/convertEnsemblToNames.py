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
import pandas as pd
import os
from stpipeline.common.gff_reader import *

def main(st_data_file, annotation, output_file):

    if not os.path.isfile(st_data_file) or not os.path.isfile(annotation):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)
        
    # loads a map with the ensembl -> gene name
    gene_map = dict()
    for line in gff_lines(annotation):
        gene_map[line["gene_id"]] = line["gene_name"]
    assert(len(gene_map) > 0)
        
    # Iterates the genes IDs to get gene names
    st_data = pd.read_table(st_data_file, sep="\t", header=0, index_col=0)
    adjustedList = list()
    for gene in st_data.columns:
        try:
            gene = gene_map[gene]
        except KeyError:
            sys.stdout.write("Warning, {} was not found in the annotation\n".format(gene))
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
    parser.add_argument("--annotation", required=True,
                        help="Path to the annotation file used to generate the data")
    args = parser.parse_args()
    main(args.st_data_file, args.annotation, args.output)

