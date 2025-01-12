#! /usr/bin/env python
"""
Script that parses a Spatial Transcriptomics (ST) data file generated
with the ST Pipeline in matrix (TSV) format where the genes are named
with ENSEMBL IDs and generates a new file with the ENSEMBL IDs converted to gene names.

The script needs the annotation file (GFF format) used to create the ST dataset with
the ST Pipeline.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import pandas as pd
import os
from stpipeline.common.gff_reader import gff_lines
from collections import Counter


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("counts_matrix", help="Matrix with gene counts (genes as columns) in TSV format")
    parser.add_argument("--output", default="output.tsv", help="Name of the output file, default output.tsv")
    parser.add_argument("--annotation", required=True, help="Path to the annotation file used to generate the data")
    args = parser.parse_args()
    sys.exit(run(args.counts_matrix, args.annotation, args.output))


def run(st_data_file: str, annotation: str, output_file: str) -> int:
    if not os.path.isfile(st_data_file) or not os.path.isfile(annotation):
        print("Error, input file(s) not present or invalid format")
        return 1

    # loads a map with the Ensembl ids -> gene name
    gene_map = {}
    for line in gff_lines(annotation):
        try:
            gene_map[line["gene_id"]] = line["gene_name"]
        except KeyError as e:
            print(f"Error, parsing annotation file, missing key {e}")
            return 1
    assert len(gene_map) > 0

    # Load the ST dataset
    st_data = pd.read_table(st_data_file, sep="\t", header=0, index_col=0)

    # Check that the annotation file given was valid
    if len(gene_map) < len(st_data.columns):
        print("Error, the annotation file given is invalid or does not match the ST data")
        return 1

    # Checks that there are no duplicated genes ids in the input data
    gene_ids_counter = Counter(st_data.columns)
    for gene_id, count in gene_ids_counter.most_common():
        if count > 1:
            print(f"Error, Ensembl ID {gene_id} was found {count} times in the input matrix.")
            return 1

    # Iterates the genes IDs to get gene names
    genes_replaced = set()
    adjustedList = []
    for gene_id in st_data.columns:
        try:
            gene_name = gene_map[gene_id]
            # Check if the gene_name has been "used" before
            if gene_name not in genes_replaced:
                genes_replaced.add(gene_name)
            else:
                # This means the gene name would be duplicated in the output
                # so we keep the Ensembl ID.
                # We assume input Ensembl ids are unique as we checked this before
                gene_name = gene_id
                print(
                    f"Warning, gene {gene_name} was already matched so the original Ensembl ID {gene_id} will be kept"
                )
        except KeyError:
            print(f"Warning, {gene_id} was not found in the annotation so the original Ensembl ID will be kept")
            gene_name = gene_id
        adjustedList.append(gene_name)

    # Update the table with the gene names
    st_data.columns = pd.Index(adjustedList)

    # Write table to file
    st_data.to_csv(output_file, sep="\t")

    return 0


if __name__ == "__main__":
    main()
