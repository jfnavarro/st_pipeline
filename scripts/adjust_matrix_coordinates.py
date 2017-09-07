#! /usr/bin/env python
""" 
Script that takes a matrix of counts
where the columns are genes and the rows
are the spot coordinates like: 
 
    gene    gene
XxY
XxY

And then removes the spots that are not present in 
a tab delimited coordinates file that has a least 4 columns:

x y new_x new_y pixel_x pixel_y

Optionally, the coordinates of the spots in the matrix
can be changed to the new coordinates (pixel or array).


@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import sys
import os
import pandas as pd

def main(counts_matrix, coordinates_file, update_coordinates, outfile, outformat):

    if not os.path.isfile(counts_matrix) or not os.path.isfile(coordinates_file):
        sys.stderr.write("Error, input file not present or invalid format\n")
        sys.exit(1)
     
    if not outfile:
        outfile = "adjusted_{}".format(os.path.basename(counts_matrix))
           
    # Get a map of the new coordinates
    new_coordinates = dict()
    with open(coordinates_file, "r") as filehandler:
        for line in filehandler.readlines():
            tokens = line.split()
            assert(len(tokens) == 6 or len(tokens) == 4)
            if tokens[0] != "x":
                old_x = int(tokens[0])
                old_y = int(tokens[1])
                new_x = round(float(tokens[2]),2)
                new_y = round(float(tokens[3]),2)
                if outformat == "array":
                    new_coordinates[(old_x, old_y)] = (new_x,new_y)
                elif len(tokens) == 6:
                    pixel_x = float(tokens[4])
                    pixel_y = float(tokens[5])
                    new_coordinates[(old_x, old_y)] = (pixel_x, pixel_y)
                else:
                    sys.stderr.write("Error, output format is pixel coordinates but\n "
                                     "the coordinates file only contains 4 columns\n")
                    sys.exit(1)                    
    # Read the data frame (spots as rows)
    counts_table = pd.read_table(counts_matrix, sep="\t", header=0, index_col=0)
    new_index_values = list()

    # Replace spot coordinates and remove row if not present
    for index in counts_table.index:
        tokens = index.split("x")
        x = int(tokens[0])
        y = int(tokens[1])
        try:
            new_x, new_y = new_coordinates[(x,y)] 
            if not update_coordinates:
                new_x, new_y = x,y
            new_index_values.append("{0}x{1}".format(new_x,new_y))
        except KeyError:
            counts_table.drop(index, inplace=True)

    # Write table again
    counts_table.index = new_index_values
    counts_table.to_csv(outfile, sep='\t')
               
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts-matrix", required=True,
                        help="Matrix with gene counts (genes as columns)")
    parser.add_argument("--outfile", help="Name of the output file")
    parser.add_argument("--update-coordinates", action="store_true", default=False,
                        help="Updates the spot coordinates in the output matrix with the\n"
                        "new coordinates present in the coordinates file")
    parser.add_argument("--coordinates-file",  required=True,
                        help="New coordinates in a tab delimited file")
    parser.add_argument("--outformat", default="array",
                        help="Output array positions or pixel positions"
                        " array or pixel", choices=["array", "pixel"])
    args = parser.parse_args()

    main(args.counts_matrix, args.coordinates_file, args.update_coordinates, args.outfile, args.outformat)


