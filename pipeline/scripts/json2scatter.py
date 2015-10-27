#! /usr/bin/env python
#@Author Jose Fernandez
""" Script for creating a quality scatter plot from a json ST-data file.
The output will be a .png file with the same name as the input file.

It allows to plot clusters generated in stclust of the following form

CLUSTER_NUMBER BARCODE X Y

It allows to choose transparency for the data points

It allows to give the high-res image 

If a regular expression for a gene symbol to highlight is provided, the output
image will be stored in a file called *.0.png. It is possible to give several
regular expressions to highlight, by adding another --highlight parameter.
The images from these queries will be stored in files with the number
increasing: *.0.png, *.1.png, *.2.png, etc.
"""

import argparse
import os
import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import transforms
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from stpipeline.common.json_utils import json_iterator
import numpy as np

blues = LinearSegmentedColormap.from_list('blues', [[0.02, 0.44, 0.69, 0.3],
                                                    [0.02, 0.44, 0.69, 1.0]])


def parseJSON(input_data, cutoff, highlight_regexes, highlights, expression):
    
    it = json_iterator(input_data)
    for doc in it:   
        
        if cutoff is not None and int(doc['hits']) < cutoff:
            continue
        
        expression[doc['x'], doc['y']] += doc['hits']
        
        if not highlight_regexes:
            continue
        
        for i, regex in enumerate(highlight_regexes):
            if re.search(regex, doc["gene"]):
                highlights[i].add((doc['x'], doc['y']))
                
    return highlights, expression
    
def parseCSV(input_data, cutoff, highlight_regexes, highlights, expression):
    
    with open(input_data,"r") as filehandler_read:
        for line in filehandler_read.readlines():
            if line.find("#") != -1:
                continue
            tokens = line.split()
            hits = int(tokens[4])
            gene = tokens[0]
            x = int(tokens[2])
            y = int(tokens[3])
            bc = tokens[1]
            
            if cutoff is not None and hits < cutoff:
                continue
            
            expression[x, y] += hits
            
            if not highlight_regexes:
                continue
            
            for i, regex in enumerate(highlight_regexes):
                if re.search(regex, gene):
                    highlights[i].add((x, y))
                    
    return highlights, expression

def main(input_data, highlight_regexes, image, only_highlight, cutoff, 
         highlight_barcodes, alignment, data_alpha, dot_size):
    fig = []
    ax = []
    highlights = []
    expression = np.zeros((33, 33), dtype=np.int)
    colors = np.zeros((33,33), dtype=np.int)
    alignment_matrix = np.zeros((3,3), dtype=np.float)

    if alignment is not None:
        alignment_matrix[0,0] = alignment[0]
        alignment_matrix[0,1] = alignment[1]
        alignment_matrix[0,2] = alignment[2]
        alignment_matrix[1,0] = alignment[3]
        alignment_matrix[1,1] = alignment[4]
        alignment_matrix[1,2] = alignment[5]
        alignment_matrix[2,0] = alignment[6]
        alignment_matrix[2,1] = alignment[7]
        alignment_matrix[2,2] = alignment[8]
        
    # Parse the clusters colors if needed
    if highlight_barcodes is not None:        
        with open(highlight_barcodes, "r") as filehandler_read:
            for line in filehandler_read.readlines():
                tokens = line.split()
                barcode = tokens[1]
                cluster = int(tokens[0])
                x = int(tokens[2])
                y = int(tokens[3])
                colors[x,y] = cluster

    for _ in highlight_regexes if highlight_regexes else [0]:
        f = plt.figure()
        a = f.add_subplot(111, aspect='equal')
        fig.append(f)
        ax.append(a)
        highlights.append(set([]))

    if input_data.endswith(".json"):
        highlights,expression = parseJSON(input_data, cutoff,
                                          highlight_regexes, highlights, expression)
    else:
        highlights,expression = parseCSV(input_data, cutoff, 
                                         highlight_regexes, highlights, expression)
     
    x, y = expression.nonzero()

    if image:
        img = plt.imread(image)

    for a in ax:
        if not only_highlight:
            base_trans = a.transData 
            tr = transforms.Affine2D(matrix = alignment_matrix) + base_trans
            
            a.scatter(x, y, c=expression[x, y],
                      edgecolor="none",
                      s=dot_size,
                      label="Expression",
                      transform = tr,
                      alpha=data_alpha) 
                                     
            if highlight_barcodes is not None:
                x2, y2 = colors.nonzero()
                sc = a.scatter(x2, y2, c=colors[x2, y2],
                               edgecolor="none",
                               s=dot_size,
                               transform = tr
                               )
        if image:
            a.imshow(img)

    for i, regex in enumerate(highlight_regexes if highlight_regexes else [0]):
        
        if len(highlights[i]) == 0:
            continue

        x, y = zip(*highlights[i])
        ax[i].scatter(x, y, c="#CA0020",
                      edgecolor="#CA0020",
                      s=dot_size + 10,
                      label=highlight_regexes[i])

    for i, a in enumerate(ax):

        a.set_xlabel("X")
        a.set_ylabel("Y")
        a.legend()
        a.set_title("Scatter", size=20)

        if highlight_regexes:
            ending = ".{0}.png".format(i)
        else:
            ending = ".png"

        img_file = "data_plot" + ending
        if image:
            fig[i].set_size_inches(16, 16)
        else:
            fig[i].set_size_inches(10, 8)

        fig[i].savefig(img_file, dpi=900)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_data", help="JSON ST-data file or CSV ST-data file")
    parser.add_argument("--highlight", help="Regular expression for \
                        gene symbols to highlight in the quality \
                        scatter plot. Can be given several times.",
                        default=None,
                        type=str,
                        action='append')
    parser.add_argument("--image", default=None)
    parser.add_argument("--only-highlight", help="Only shows highlighted genes",
                        default=False, action='store_true')
    parser.add_argument("--cutoff", help="Do not include genes below this reads cut off",
                        type=int, default=None)
    parser.add_argument("--highlight-barcodes", default=None,
                        help="File with a list of barcodes in a column and a list of clusters in the other column")
    parser.add_argument("--alignment", help="Aligment matrix needed when using the image", nargs="+", type=float, default=None)
    parser.add_argument("--data-alpha", type=float, default=1.0, 
                        help="The brightness level for the data points, 0 min and 1 max")
    parser.add_argument("--dot-size", type=int, default=30,
                        help="The size of the dots")
    args = parser.parse_args()

    main(args.input_data, args.highlight, args.image, 
         args.only_highlight, args.cutoff, 
         args.highlight_barcodes, args.alignment, float(args.data_alpha), int(args.dot_size))
