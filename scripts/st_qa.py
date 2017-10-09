#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script that performs a basic Quality Control analysis 
of a ST dataset (matrix in TSV) format.

The scripts writes stats and generates
some plots in the folder that is run. 

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""
import pandas as pd
import numpy as np
import os.path
import argparse
import matplotlib.pyplot as plt

def scatter_plot(x_points, y_points, output, colors,
                 title="Scatter", xlabel="X", ylabel="Y",
                 xlim=[1,33], ylim=[1,35]):
    """ 
    This function makes a scatter plot of a set of points (x,y).
    and a given list of color values for each point.
    The plot will be written to a file.
    :param x_points: a list of x coordinates
    :param y_points: a list of y coordinates
    :param output: the name/path of the output file
    :param colors: a color value for each point
    :param title: the title for the plot
    :param xlabel: the name of the X label
    :param ylabel: the name of the Y label
    :raises: RuntimeError
    """
    # Plot spots with the color class in the tissue image
    fig = plt.figure()
    plt.scatter(x_points, 
              y_points,  
              c=colors, 
              cmap=plt.get_cmap("YlOrBr"), 
              edgecolor="none", 
              s=50)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.colorbar()
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    fig.savefig(output, format='pdf', dpi=300)
    
def histogram(x_points, output, title="Histogram", xlabel="X",
              ylabel="Y", nbins=50, color="blue"):

    """ This function generates a simple histogram
    with the points given as input.
    :param x_points: a list of x coordinates
    :param title: the title for the plot
    :param xlabel: the name of the X label
    :param ylabel: the name of the X label
    :param output: the name/path of the output file
    :param nbins: the number of bings for the histogram
    :param color: the color for the histogram
    """
    # Create the plot
    fig = plt.figure()
    plt.hist(x_points, bins=nbins, facecolor=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    fig.savefig(output, format='pdf', dpi=300)

def main(input_data):
    # Parse the data
    counts_table = pd.read_table(input_data, sep="\t", header=0, index_col=0)
    # Get the basename
    input_name = os.path.basename(input_data)
    # Compute some statistics
    total_barcodes = len(counts_table.index)
    total_transcripts = np.sum(counts_table.values, dtype=np.int32)
    number_genes = len(counts_table.columns)
    max_count = counts_table.values.max()
    min_count = counts_table.values.min()
    aggregated_spot_counts = counts_table.sum(axis=1).values
    aggregated_gene_counts = (counts_table != 0).sum(axis=1).values
    max_genes_feature = aggregated_gene_counts.max()
    min_genes_feature = aggregated_gene_counts.min()
    max_reads_feature = aggregated_spot_counts.max()
    min_reads_feature = aggregated_spot_counts.min()
    average_reads_feature = np.mean(aggregated_spot_counts)
    average_genes_feature = np.mean(aggregated_gene_counts)
    std_reads_feature = np.std(aggregated_spot_counts)
    std_genes_feature = np.std(aggregated_gene_counts)
    # Generate heatmap plots
    histogram(aggregated_spot_counts, nbins=20, xlabel="#Reads", ylabel="#Spots",
              output=input_name + "_hist_counts.pdf", title="Reads per spot")
    histogram(aggregated_gene_counts, nbins=20, xlabel="#Genes", ylabel="#Spots", 
              output=input_name + "_hist_genes.pdf", title="Genes per spot")
    # Get the spot coordinates
    x_points = list()
    y_points = list()
    for spot in counts_table.index:
        tokens = spot.split("x")
        assert(len(tokens) == 2)
        y_points.append(float(tokens[1]))
        x_points.append(float(tokens[0]))
    # Generate scater plots
    scatter_plot(x_points, y_points, colors=aggregated_spot_counts, 
                 xlabel="X", ylabel="Y", output=input_name + "_heatmap_counts.pdf", 
                 title="Heatmap expression")
    scatter_plot(x_points, y_points, colors=aggregated_gene_counts, 
                 xlabel="X", ylabel="Y", output=input_name + "_heatmap_genes.pdf", 
                 title="Heatmap genes")

    qa_stats = [
    ("Number of features: {}".format(total_barcodes)+"\n"),
    ("Number of unique molecules present: {}".format(total_transcripts)+"\n"),
    ("Number of unique genes present: {}".format(number_genes)+"\n"),
    ("Max number of genes over all spots: {}".format(max_genes_feature)+"\n"),
    ("Min number of genes over all spots: {}".format(min_genes_feature)+"\n"),
    ("Max number of unique molecules over all spots: {}".format(max_reads_feature)+"\n"),
    ("Min number of unique molecules over all spots: {}".format(min_reads_feature)+"\n"),
    ("Average number genes per spots: {}".format(average_genes_feature)+"\n"),
    ("Average number unique molecules per spot: {}".format(average_reads_feature)+"\n"),
    ("Std number genes per spot: {}".format(std_genes_feature)+"\n"),
    ("Std number unique molecules per spot: {}".format(std_reads_feature)+"\n"),
    ("Max number of unique molecules over all unique events: {}".format(max_count)+"\n"),
    ("Min number of unique molecules over all unique events: {}".format(min_count)+"\n")
    ]
    # Print stats to stdout and a file
    print("".join(qa_stats))
    with open("{}_qa_stats.txt".format(input_name), "a") as outfile:
        outfile.write("".join(qa_stats))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-data", required=True, type=str,
                        help="The matrix of counts (spots as row names and genes as column names)")
    args = parser.parse_args()
    main(args.input_data)
