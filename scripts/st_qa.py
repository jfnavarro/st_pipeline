#! /usr/bin/env python
"""
Script that performs a basic Quality Control analysis
of a Spatial Transcriptomics dataset (matrix in TSV) format.

The script writes stats and generates some plots in the folder
where it is run.

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from typing import List
import pandas as pd
import numpy as np
import numpy.typing as npt
import os.path
import argparse
import matplotlib.pyplot as plt
import seaborn as sns  # type: ignore
import sys


def scatter_plot(
    x_points: List[float],
    y_points: List[float],
    output: str,
    colors: npt.NDArray[np.int32],
    title: str = "Scatter",
    xlabel: str = "X",
    ylabel: str = "Y",
) -> None:
    """
    Creates a scatter plot of a set of points (x, y) with corresponding color values
    and saves it as a PDF.

    Args:
        x_points: A 1D array of x coordinates.
        y_points: A 1D array of y coordinates.
        output: The file path where the plot will be saved.
        colors: A 1D array of color values for each point.
        title: The title of the plot. Defaults to "Scatter".
        xlabel: The label for the x-axis. Defaults to "X".
        ylabel: The label for the y-axis. Defaults to "Y".

    Returns:
        None

    Raises:
        RuntimeError: If an error occurs during plot creation or saving.
    """
    try:
        fig = plt.figure()
        plt.scatter(x_points, y_points, c=colors, cmap=plt.get_cmap("YlOrBr"), edgecolor="none", s=50)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.gca().invert_yaxis()
        plt.title(title)
        plt.colorbar()
        plt.subplots_adjust(left=0.15)
        fig.savefig(output, format="pdf", dpi=300)
    except Exception as e:
        raise RuntimeError("Failed to create scatter plot") from e


def histogram(
    x_points: npt.NDArray[np.int32],
    output: str,
    title: str = "Histogram",
    xlabel: str = "X",
    ylabel: str = "Y",
    nbins: int = 50,
    color: str = "blue",
) -> None:
    """
    Generates a histogram with the given data points and saves it as a PDF.

    Args:
        x_points: A list of x coordinates to be plotted.
        output: The file path where the plot will be saved.
        title: The title of the plot. Defaults to "Histogram".
        xlabel: The label for the x-axis. Defaults to "X".
        ylabel: The label for the y-axis. Defaults to "Y".
        nbins: The number of bins for the histogram. Defaults to 50.
        color: The color of the histogram. Defaults to "blue".

    Returns:
        None

    Raises:
        RuntimeError: If an error occurs during plot creation or saving.
    """
    try:
        fig = plt.figure()
        plt.hist(x_points, bins=nbins, facecolor=color)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.subplots_adjust(left=0.15)
        fig.savefig(output, format="pdf", dpi=300)
    except Exception as e:
        raise RuntimeError("Failed to create histogram") from e


def main(input_data: str) -> int:
    # Parse the data
    counts_table = pd.read_table(input_data, sep="\t", header=0, index_col=0)

    # Get the basename
    input_name = os.path.basename(input_data).split(".")[0]

    # Compute some statistics
    total_barcodes = len(counts_table.index)
    total_transcripts = np.sum(counts_table.values, dtype=np.int32)
    number_genes = len(counts_table.columns)
    max_count = counts_table.max()
    min_count = counts_table.min()
    aggregated_spot_counts = counts_table.sum(axis=1).to_numpy()
    aggregated_gene_counts = (counts_table > 0).sum(axis=1).to_numpy()
    aggregated_gene_counts_1 = (counts_table > 1).sum(axis=1).to_numpy()
    aggregated_gene_counts_2 = (counts_table > 2).sum(axis=1).to_numpy()
    aggregated_gene_gene_counts = (counts_table > 0).sum(axis=0).to_numpy()
    aggregated_gene_gene_counts_1 = (counts_table > 1).sum(axis=0).to_numpy()
    aggregated_gene_gene_counts_2 = (counts_table > 2).sum(axis=0).to_numpy()
    max_genes_feature = aggregated_gene_counts.max()
    min_genes_feature = aggregated_gene_counts.min()
    max_reads_feature = aggregated_spot_counts.max()
    min_reads_feature = aggregated_spot_counts.min()
    average_reads_feature = np.mean(aggregated_spot_counts)
    average_genes_feature = np.mean(aggregated_gene_counts)
    std_reads_feature = np.std(aggregated_spot_counts)
    std_genes_feature = np.std(aggregated_gene_counts)

    # Generate heatmap plots
    histogram(
        aggregated_spot_counts,
        nbins=20,
        xlabel="#Reads",
        ylabel="#Spots",
        output=input_name + "_hist_reads_spot.pdf",
        title="Reads per spot",
    )
    histogram(
        aggregated_gene_counts,
        nbins=20,
        xlabel="#Genes",
        ylabel="#Spots",
        output=input_name + "_hist_genes_spot.pdf",
        title="Genes per spot (>0)",
    )
    histogram(
        aggregated_gene_counts_1,
        nbins=20,
        xlabel="#Genes",
        ylabel="#Spots",
        output=input_name + "_hist_genes_spots_1.pdf",
        title="Genes per spot (>1)",
    )
    histogram(
        aggregated_gene_counts_2,
        nbins=20,
        xlabel="#Genes",
        ylabel="#Spots",
        output=input_name + "_hist_genes_spots_2.pdf",
        title="Genes per spot (>2)",
    )
    histogram(
        aggregated_gene_gene_counts,
        nbins=20,
        xlabel="#Spots",
        ylabel="#Genes",
        output=input_name + "_hist_spots_gene.pdf",
        title="Spots per gene (>0)",
    )
    histogram(
        aggregated_gene_gene_counts_1,
        nbins=20,
        xlabel="#Spots",
        ylabel="#Genes",
        output=input_name + "_hist_spots_gene_1.pdf",
        title="Spots per gene (>1)",
    )
    histogram(
        aggregated_gene_gene_counts_2,
        nbins=20,
        xlabel="#Spots",
        ylabel="#Genes",
        output=input_name + "_hist_spots_gene_2.pdf",
        title="Spots per gene (>2)",
    )
    plt.clf()

    # Generate density plots
    sns.displot(aggregated_gene_counts, kind="kde", label="Counts > 0")
    sns.displot(aggregated_gene_counts_1, kind="kde", label="Counts > 1")
    sns_plot = sns.displot(aggregated_gene_counts_2, axlabel="#Genes", kind="kde", label="Counts > 2")
    fig = sns_plot.get_figure()
    fig.savefig(input_name + "_density_genes_by_spot.pdf")
    plt.clf()

    sns.displot(aggregated_gene_gene_counts, kind="kde", label="Counts > 0")
    sns.displot(aggregated_gene_gene_counts_1, kind="kde", label="Counts > 1")
    sns_plot = sns.displot(aggregated_gene_gene_counts_2, axlabel="#Spots", kind="kde", label="Counts > 2")
    fig = sns_plot.get_figure()
    fig.savefig(input_name + "_density_spots_by_gene.pdf")
    plt.clf()

    sns.scatterplot(x=aggregated_spot_counts, y=aggregated_gene_counts, label="Gene counts >0")
    sns.scatterplot(x=aggregated_spot_counts, y=aggregated_gene_counts_1, label="Gene counts >1")
    sns_plot = sns.scatterplot(x=aggregated_spot_counts, y=aggregated_gene_counts_2, label="Gene counts >2")
    fig = sns_plot.get_figure()
    fig.savefig(input_name + "_scatter_reads_vs_genes.pdf")
    plt.clf()

    # sns_plot = sns.jointplot(x=aggregated_spot_counts,
    #                         y=aggregated_gene_counts, kind='kde', color="skyblue")
    # fig = sns_plot.get_figure()
    # fig.savefig(input_name + "_join_density_reads_vs_genes.pdf")
    # plt.clf()

    qa_stats = [
        ("Number of spots: {}".format(total_barcodes) + "\n"),
        ("Number of reads present: {}".format(total_transcripts) + "\n"),
        ("Number of unique genes present: {}".format(number_genes) + "\n"),
        ("Max number of genes over all spots: {}".format(max_genes_feature) + "\n"),
        ("Min number of genes over all spots: {}".format(min_genes_feature) + "\n"),
        ("Max number of reads over all spots: {}".format(max_reads_feature) + "\n"),
        ("Min number of reads over all spots: {}".format(min_reads_feature) + "\n"),
        ("Average number genes per spots: {}".format(average_genes_feature) + "\n"),
        ("Average number reads per spot: {}".format(average_reads_feature) + "\n"),
        ("Std number genes per spot: {}".format(std_genes_feature) + "\n"),
        ("Std number reads per spot: {}".format(std_reads_feature) + "\n"),
        ("Max number of reads over all spots/genes: {}".format(max_count) + "\n"),
        ("Min number of reads over all spots/genes: {}".format(min_count) + "\n"),
    ]
    # Print stats to stdout and a file
    print("".join(qa_stats))
    with open("{}_qa_stats.txt".format(input_name), "a") as outfile:
        outfile.write("".join(qa_stats))

    # Generate scatter plots
    # Get the spot coordinates
    x_points = []
    y_points = []
    for spot in counts_table.index:
        tokens = spot.split("x")
        assert len(tokens) == 2
        y_points.append(float(tokens[1]))
        x_points.append(float(tokens[0]))
    scatter_plot(
        x_points,
        y_points,
        colors=aggregated_spot_counts,
        xlabel="X",
        ylabel="Y",
        output=input_name + "_heatmap_counts.pdf",
        title="Heatmap expression",
    )
    scatter_plot(
        x_points,
        y_points,
        colors=aggregated_gene_counts,
        xlabel="X",
        ylabel="Y",
        output=input_name + "_heatmap_genes.pdf",
        title="Heatmap genes",
    )

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("counts_matrix", help="Matrix with gene counts (genes as columns) in TSV format")
    args = parser.parse_args()
    sys.exit(main(args.counts_matrix))
