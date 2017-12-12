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
    inverted_gene_map = dict()
    non_unique_gene_names = dict()

    for line in gff_lines(annotation):

        if line["gene_name"] in inverted_gene_map and line["gene_id"] != inverted_gene_map[line["gene_name"]]:
            # This will introduce a duplication of gene_name in the st_data_file if both the gene_ids are present within the inputfile
            try: non_unique_gene_names[line["gene_name"]][line["gene_id"]] = False
            except KeyError: non_unique_gene_names[line["gene_name"]]={line["gene_id"]:False}

        gene_map[line["gene_id"]] = line["gene_name"]
        inverted_gene_map[line["gene_name"]] = line["gene_id"]

    assert(len(gene_map) > 0)
    sys.stdout.write("Info, annotation file loaded.\n")

    # Iterates the genes IDs to get gene names
    st_data = pd.read_table(st_data_file, sep="\t", header=0, index_col=0)

    gene_ids_counter = Counter(st_data.columns)
    for gene_id, count in gene_ids_counter.most_common():
        if count == 1: break
        sys.stdout.write("Warning, gene_id {} was found {} times in the input matrix.\n".format(gene_id, count))

    adjustedList = list()
    for gene_id in st_data.columns:
        try:
            gene_name = gene_map[gene_id]

            if gene_name in non_unique_gene_names:
                non_unique_gene_names[gene_name][gene_id] = True
                if sum([1 for ensemblId in non_unique_gene_names[gene_name] if non_unique_gene_names[gene_name][ensemblId] ]) > 1:
                    sys.stdout.write(
                        "Warning, gene_name {} was found more than once in the annotation and will be duplicated in the st data file (the following EnsmblIds are affected: {}).\n".format(
                            gene_name,
                            ','.join(
                                [ensemblId for ensemblId in non_unique_gene_names[gene_name] if non_unique_gene_names[gene_name][ensemblId] ]
                                )
                            )
                        )

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

