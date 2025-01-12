#! /usr/bin/env python
"""
Script that takes a matrix of counts
where the columns are genes and the rows
are spot coordinates

        gene    gene
XxY
XxY

And removes the columns of genes whose functional type is not in the
allowed types provided (Ensembl annotation gene type).

The script needs to be given an annotation file in GFF format.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import os
import pandas as pd
from stpipeline.common.gff_reader import gff_lines
from typing import List


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("counts_matrix", help="Matrix with gene counts (genes as columns) in TSV format")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument(
        "--gene-types-keep",
        required=True,
        nargs="+",
        type=str,
        help="List of Ensembl gene types to keep (E.x protein_coding lincRNA",
    )
    parser.add_argument("--annotation", help="The Ensembl annotation file", required=True, type=str)
    parser.add_argument(
        "--ensembl-ids",
        action="store_true",
        default=False,
        help="Pass this parameter if the genes in the matrix " "are named with Ensembl Ids instead of gene names",
    )
    args = parser.parse_args()
    sys.exit(run(args.counts_matrix, args.gene_types_keep, args.outfile, args.annotation, args.ensembl_ids))


def run(counts_matrix: str, gene_types_keep: List[str], outfile: str, annotation: str, ensembl_ids: bool) -> int:
    if not os.path.isfile(counts_matrix) or not os.path.isfile(annotation):
        print("Error, input file(s) not present or invalid format")
        return 1

    if not outfile:
        outfile = "filtered_{}".format(os.path.basename(counts_matrix))

    gene_types = {}
    for line in gff_lines(annotation):
        try:
            gene_name = line["gene_id"] if ensembl_ids else line["gene_name"]
            gene_types[gene_name] = line["gene_type"] if "gene_type" in line else line["gene_biotype"]
        except KeyError as e:
            print("Error, parsing annotation file, missing key {}".format(e))
    assert len(gene_types) > 0

    # Read the data frame (genes as columns)
    counts_table = pd.read_table(counts_matrix, sep="\t", header=0, index_col=0)
    genes = counts_table.columns

    # Filter out genes that match any of the allowed types
    genes_drop = []
    for gene in genes:
        try:
            if gene_types[gene] not in gene_types_keep:
                genes_drop.append(gene)
        except KeyError:
            print(f"Warning, {gene} was not found in the annotation")

    if len(genes_drop) > 0:
        counts_table.drop(genes_drop, axis=1, inplace=True)
    else:
        print("Not a single gene could be discarded...")

    # Write filtered table
    counts_table.to_csv(outfile, sep="\t")

    return 0


if __name__ == "__main__":
    main()
