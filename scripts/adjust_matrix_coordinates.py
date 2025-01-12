#! /usr/bin/env python
"""
Script that takes a matrix of counts
where the columns are genes and the rows
are the spot coordinates in a format like this:

    gene    gene
XxY
XxY

And then removes the spots that are not present in
the spot coordinates file or if the under_tissue flag is 0
The format of the spot coordinates file can be like this:

x y new_x new_y

or

x y new_x new_y pixel_x pixel_y

or

x y new_x new_y pixel_x pixel_y under_tissue

Optionally, the coordinates of the spots in the matrix
can be changed to the adjusted new coordinates (array).

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import sys
import os
import pandas as pd


def main(counts_matrix: str, coordinates_file: str, update_coordinates: bool, outfile: str) -> int:
    if not os.path.isfile(counts_matrix) or not os.path.isfile(coordinates_file):
        print("Error, input file(s) not present or invalid format")
        return 1

    if not outfile:
        outfile = f"adjusted_{os.path.basename(counts_matrix)}"

    # Get a map of the new coordinates
    new_coordinates = {}
    with open(coordinates_file, "r") as filehandler:
        for line in filehandler.readlines():
            tokens = line.split()
            assert len(tokens) == 6 or len(tokens) == 4 or len(tokens) == 7
            if tokens[0] != "x":
                old_x = int(tokens[0])
                old_y = int(tokens[1])
                new_x = round(float(tokens[2]), 2)
                new_y = round(float(tokens[3]), 2)
                if len(tokens) == 7 and not bool(tokens[6]):
                    continue
                new_coordinates[(old_x, old_y)] = (new_x, new_y)

    # Read the data frame (spots as rows)
    counts_table = pd.read_table(counts_matrix, sep="\t", header=0, index_col=0)
    new_index_values = []

    # Replace spot coordinates and remove row if not present
    for index in counts_table.index:
        tokens = index.split("x")
        x = int(tokens[0])
        y = int(tokens[1])
        try:
            new_x, new_y = new_coordinates[(x, y)]
            if not update_coordinates:
                new_x, new_y = x, y
            new_index_values.append(f"{new_x}x{new_y}")
        except KeyError:
            counts_table.drop(index, inplace=True)

    # Assign the new indexes
    counts_table.index = pd.Index(new_index_values)

    # Remove genes that have now a total count of zero
    counts_table = counts_table.transpose()[counts_table.sum(axis=0) > 0].transpose()

    # Write table again
    counts_table.to_csv(outfile, sep="\t")

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("counts_matrix", help="Matrix with gene counts (genes as columns) in TSV format")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument(
        "--update-coordinates",
        action="store_true",
        default=False,
        help="Updates the spot coordinates in the output matrix with the\n"
        "new coordinates present in the coordinates file",
    )
    parser.add_argument("--coordinates-file", required=True, help="New coordinates in a tab delimited file")
    args = parser.parse_args()

    sys.exit(main(args.counts_matrix, args.coordinates_file, args.update_coordinates, args.outfile))
