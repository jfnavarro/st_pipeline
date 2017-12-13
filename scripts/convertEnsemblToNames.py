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
from collections import Counter

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

    gene_ids_counter = Counter(st_data.columns)
    for gene_id, count in gene_ids_counter.most_common():
        if count == 1: break
        sys.stdout.write("Warning, gene_id {} was found {} times in the input matrix.\n".format(gene_id, count))

    genes_replaced = dict()

    adjustedList = list()
    for gene_id in st_data.columns:
        try:
            gene_name = gene_map[gene_id]

            # check if the gene_name has been "used" before
            if gene_name not in genes_replaced:
                genes_replaced[gene_name] = [gene_id]

            else: # this means the gene name will be duplicated in the output
                message = 'Warning, gene_name {} is already present in converted set.'.format(gene_name)
                genes_replaced[gene_name].append(gene_id)

                # check if the corresponding gene id was duplicated in the ST data input file,
                # in this case it will be dupluicated in output as well
                if gene_ids_counter[gene_id] > 1:
                    message += ' The corresponding gene_id was present {} times in the input ST data file.'.format( gene_ids_counter[gene_id] )

                # On the otherhand if the gene id was unique in the ST data input there must >1 ids mapping to the same name
                # in this case the original (unique) gene id will be used instead of the gene name
                else:
                    message += ' The gene_id {} was unique in the input ST data file.'.format( gene_id )
                    message += ' The following already converted gene_ids also map to the gene_name {}: {}.'.format( gene_name, ', '.join(genes_replaced[gene_name][:-1]) )
                    message += ' The original gene_id ({}) will be used in the output file for the entry.'.format(gene_id)
                    gene_name = gene_id

                sys.stdout.write('{}\n'.format(message))

        except KeyError:
            sys.stdout.write("Warning, {} was not found in the annotation\n".format(gene_id))
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

