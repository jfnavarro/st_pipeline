#! /usr/bin/env python
#@Author Jose Fernandez
""" Script for creating a quality scatter plot from a json ST-data file.
The output will be a .png file with the same name as the input file.

It allows to plot clusters generated in stclust of the following form

CLUSTER_NUMBER BARCODE X Y

It allows to choose transparency for the data points

It allows to give the high-res image 

It allows to give the normalized counts from stclust

It allows to filter out by counts or gene names

If a regular expression for a gene symbol to highlight is provided, the output
image will be stored in a file called *.0.png. It is possible to give several
regular expressions to highlight, by adding another --highlight parameter.
The images from these queries will be stored in files with the number
increasing: *.0.png, *.1.png, *.2.png, etc.
"""

import argparse
import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import transforms
from matplotlib import pyplot as plt
from stpipeline.common.json_utils import json_iterator
import numpy as np
import pandas as pd
import math

## This must be the same as R colors from numbers
color_map = ["white", "black", "red", "green", "blue", "cyan", "pink", "yellow", "grey"]

def parseJSON(input_data, cutoff, highlight_regexes, highlights, 
              expression, filter_genes, norm_counts_table, use_norm=False):
    
    it = json_iterator(input_data)
    for doc in it:   
        
        if cutoff and int(doc['hits']) < cutoff:
            continue
        
            filter_passed = True
            for regex in filter_genes if filter_genes else []:
                if not re.match(regex, doc["gene"]):
                    filter_passed = False
                    continue
              
            if not filter_passed:
                continue
                
        if not use_norm:
            expression[doc['x'], doc['y']] += doc['hits']
        else:
            barcode = doc['barcode']
            if barcode in norm_counts_table.columns:
                ##TODO this is done several times for same barcode-gene
                if doc["gene"] in norm_counts_table[barcode]:
                    norm_exp = norm_counts_table[barcode][doc["gene"]]
                    expression[doc['x'], doc['y']] = math.log(norm_exp,2)
                else:
                    continue
            else:
                continue
            
        if not highlight_regexes:
            continue
        
        for i, regex in enumerate(highlight_regexes):
            if re.search(regex, doc["gene"]):
                highlights[i].add((doc['x'], doc['y']))
                
    return highlights, expression
    
def parseCSV(input_data, cutoff, highlight_regexes, highlights, 
             expression, filter_genes, norm_counts_table, use_norm=False):
    
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
            
            if cutoff and hits < cutoff:
                continue
            
            filter_passed = True
            for regex in filter_genes if filter_genes else []:
                if not re.match(regex, gene):
                    filter_passed = False
                    continue
              
            if not filter_passed:
                continue
                               
            if not use_norm:  
                expression[x, y] += hits
            else:
                ##TODO this is done several times for same barcode-gene
                if bc in norm_counts_table.columns:
                    if gene in norm_counts_table[bc]:
                        norm_exp = norm_counts_table[bc][gene]
                        expression[x, y] = math.log(norm_exp,2)
                    else:
                        continue
                else:
                    continue
            
            if not highlight_regexes:
                continue
            
            for i, regex in enumerate(highlight_regexes):
                if re.search(regex, gene):
                    highlights[i].add((x, y))
                    
    return highlights, expression

def main(input_data, 
         highlight_regexes, 
         image,
         cutoff, 
         highlight_barcodes, 
         alignment, 
         data_alpha,
         highlight_alpha,
         dot_size,
         normalized_counts,
         filter_genes):
    
    #TODO add checks for parameters
    
    fig = []
    ax = []
    highlights = []
    expression = np.zeros((33, 33), dtype=np.float)
    colors = np.zeros((33,33), dtype=np.int)
    alignment_matrix = np.zeros((3,3), dtype=np.float)
    # load normalized counts
    use_norm = False
    norm_counts_table = None
    if normalized_counts is not None:
        norm_counts_table = pd.read_table(normalized_counts, sep="\t", header=0)
        use_norm = True
    
    # Create alignment matrix 
    alignment_matrix[0,0] = 1
    alignment_matrix[0,1] = 0
    alignment_matrix[0,2] = 0
    alignment_matrix[1,0] = 0
    alignment_matrix[1,1] = 1
    alignment_matrix[1,2] = 0
    alignment_matrix[2,0] = 0
    alignment_matrix[2,1] = 0
    alignment_matrix[2,2] = 1
    if alignment:
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
    if highlight_barcodes:        
        with open(highlight_barcodes, "r") as filehandler_read:
            for line in filehandler_read.readlines():
                tokens = line.split()
                barcode = tokens[1]
                cluster = int(tokens[0])
                x = int(tokens[2])
                y = int(tokens[3])
                colors[x,y] = cluster

    # Create figures
    for _ in highlight_regexes if highlight_regexes else [0]:
        f = plt.figure()
        a = f.add_subplot(111, aspect='equal')
        fig.append(f)
        ax.append(a)
        highlights.append(set([]))

    # Parse input data
    if input_data.endswith(".json"):
        highlights,expression = parseJSON(input_data, cutoff,
                                          highlight_regexes, 
                                          highlights, expression, 
                                          filter_genes, norm_counts_table, use_norm)
    else:
        highlights,expression = parseCSV(input_data, cutoff, 
                                         highlight_regexes, 
                                         highlights, expression, 
                                         filter_genes, norm_counts_table, use_norm)
     
    x, y = expression.nonzero()

    if image:
        img = plt.imread(image)

    for i, a in enumerate(ax):
        base_trans = a.transData
        tr = transforms.Affine2D(matrix = alignment_matrix) + base_trans
        
        a.scatter(x, y,
                  c=expression[x, y],
                  edgecolor="none",
                  cmap=plt.get_cmap("YlOrBr"),
                  s=dot_size,
                  transform=tr,
                  alpha=data_alpha) 
                                     
        if highlight_barcodes:
            x2, y2 = colors.nonzero()
            color_list = set(colors[x2,y2].tolist())
            cmap = color_map[min(color_list):max(color_list)+1]
            sc = a.scatter(x2, y2,
                           c=colors[x2, y2],
                           cmap=matplotlib.colors.ListedColormap(cmap),
                           edgecolor="none",
                           s=dot_size,
                           transform=tr,
                           alpha=highlight_alpha)
            #Show legend with clusters
                
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

        fig[i].savefig(img_file, dpi=300)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_data", help="JSON ST-data file or CSV ST-data file")
    parser.add_argument("--highlight", help="Regular expression for \
                        gene symbols to highlight in the quality \
                        scatter plot. Can be given several times.",
                        default=None,
                        type=str,
                        action='append')
    parser.add_argument("--image", default=None, 
                        help="When given the data will plotted on top of the image, \
                        if the alignment matrix is given the data will be aligned")
    parser.add_argument("--cutoff", help="Do not include genes below this reads cut off",
                        type=int, default=None)
    parser.add_argument("--highlight-barcodes", default=None,
                        help="File with a list of barcodes in a column and a list of clusters in the other column")
    parser.add_argument("--alignment", 
                        help="Alignment matrix needed when using the image", 
                        nargs="+", type=float, default=None)
    parser.add_argument("--data-alpha", type=float, default=1.0, 
                        help="The transparency level for the data points, 0 min and 1 max")
    parser.add_argument("--highlight-alpha", type=float, default=1.0, 
                        help="The transparency level for the highlighted barcodes, 0 min and 1 max")
    parser.add_argument("--dot-size", type=int, default=30,
                        help="The size of the dots")
    parser.add_argument("--normalized-counts", default=None, type=str,
                        help="A table with the normalized counts for each barcode-gene combination to be used \
                        instead of the reads counts. Note that if a gene-barcode is not present in the table it will \
                        not be shown")
    parser.add_argument("--filter-genes", help="Regular expression for \
                        gene symbols to filter out. Can be given several times.",
                        default=None,
                        type=str,
                        action='append')
    args = parser.parse_args()

    main(args.input_data, 
         args.highlight, 
         args.image,
         args.cutoff, 
         args.highlight_barcodes, 
         args.alignment, 
         float(args.data_alpha),
         float(args.highlight_alpha),
         int(args.dot_size),
         args.normalized_counts,
         args.filter_genes)
