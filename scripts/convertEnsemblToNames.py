#! /usr/bin/env python
"""
Script that parses a ST data file generated
with the pipeline in matrix (TSV) format where the genes are named
with ENSEMBL IDs and generates a new file
with the ENSEMBL IDs converted to gene names IDS.
For that the script also needs the annotation file (GFF format) used to create
the ST dataset.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import pandas as pd
import os
from stpipeline.common.gff_reader import *
from collections import Counter

def main(st_data_file, annotation, output_file):

    if not os.path.isfile(st_data_file) or not os.path.isfile(annotation):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)

    # loads a map with the ensembl ids -> gene name
    gene_map = dict()
    for line in gff_lines(annotation):
        gene_map[line["gene_id"]] = line["gene_name"]
    
    # Load the ST dataset
    st_data = pd.read_table(st_data_file, sep="\t", header=0, index_col=0)

    # Check that the annotation file given was valid
    if len(gene_map) < len(st_data.columns):
        sys.stdout.write("Error, the annotation file given is invalid\n")
        sys.exit(1)
            
    # Checks that there are no duplicated genes
    gene_ids_counter = Counter(st_data.columns)
    for gene_id, count in gene_ids_counter.most_common():
        if count > 1:
            sys.stdout.write("Error, Ensembl ID {} was found {} times in the input " \
                             "matrix.\n".format(gene_id, count))
            sys.exit(1)

    # Iterates the genes IDs to get gene names
    genes_replaced = set()
    adjustedList = list()
    for gene_id in st_data.columns:
        try:
            gene_name = gene_map[gene_id]
            # Check if the gene_name has been "used" before
            if gene_name not in genes_replaced:
                genes_replaced.add(gene_name)
            else:
                # This means the gene name would be duplicated in the output
                # so we keep the Ensmelb ID.
                # We assume input Ensemgl ids are unique as we checked this before
                gene_name = gene_id
                sys.stdout.write("Warning, gene name {} was already matched so the original " \
                                 "Ensembl ID {} will be kept\n".format(gene_name, gene_id))
        except KeyError:
            sys.stdout.write("Warning, {} was not found in the annotation, " \
                             "so the original Ensembl ID will be kept\n".format(gene_id))
            gene_name = gene_id
        adjustedList.append(gene_name)

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

