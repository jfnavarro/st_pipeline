#! /usr/bin/env python
"""
Script that creates multiple QC plots and stats
from multiple Spatial Transcriptomics datasets in TSV format.

This is useful when one have several consecutive sections
or sections from the same model.

The tool generates:

- Violin plots (genes and counts)
- Genes shared % pair-wise matrix
- Correlation pair-wise matrix
- Correlation pair-wise scatter plots
- PCA plot (one dot per dataset)

@Author Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr  # type: ignore
from sklearn.decomposition import PCA  # type: ignore
import os
import sys
from typing import List, Any

color_map = [
    "red",
    "green",
    "blue",
    "orange",
    "cyan",
    "yellow",
    "orchid",
    "saddlebrown",
    "darkcyan",
    "gray",
    "darkred",
    "darkgreen",
    "darkblue",
    "antiquewhite",
    "bisque",
    "black",
    "slategray",
    "gold",
    "floralwhite",
    "aliceblue",
    "plum",
    "cadetblue",
    "coral",
    "olive",
    "khaki",
    "lightsalmon",
]


def create_violin_plot(data: List[List[float]], pos: List[int], title: str, outfile: str) -> None:
    """
    Creates a violin plot and saves it as a PDF.

    Args:
        data: The data to plot, where each sublist represents a dataset.
        pos: Positions of the datasets on the x-axis.
        title: The title of the plot.
        outfile: The file path where the plot will be saved.

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    ax.violinplot(data, pos, showmeans=True, showextrema=True, showmedians=True)
    ax.set_axisbelow(True)
    ax.set_title(title)
    fig.savefig(outfile, format="pdf", dpi=90)


def create_pca_plot(data: Any, labels: List[str], title: str, outfile: str) -> None:
    """
    Creates a PCA scatter plot and saves it as a PDF.

    Args:
        data: A 2D array where each row represents a data point, and the columns are the PCA components.
        labels: Labels corresponding to the data points for the legend.
        title: The title of the plot.
        outfile: The file path where the plot will be saved.

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    class_colours = [color_map[i] for i in range(len(data))]
    ax.scatter(data[:, 0], data[:, 1], s=20, c=class_colours, edgecolor="none")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(title)
    recs = [mpatches.Rectangle((0, 0), 1, 1, fc=color) for color in class_colours]
    ax.legend(recs, labels, loc=4, prop={"size": 6})
    fig.savefig(outfile, format="pdf", dpi=90)


def main(counts_table_files: List[str], outdir: str, use_log: bool) -> int:
    if len(counts_table_files) == 0 or any(not os.path.isfile(f) for f in counts_table_files):
        print("Error, input file(s) not present or invalid format")
        return 1

    if len(counts_table_files) < 2:
        print("Error, minimum number of input datasets is 2")
        return 1

    if outdir is None or not os.path.isdir(outdir):
        outdir = os.getcwd()
    outdir = os.path.abspath(outdir)

    print(f"Output directory {outdir}")
    print(f"Input datasets {' '.join(counts_table_files)}")

    # Parse datasets and sort them by column
    datasets = [
        (pd.read_table(x, sep="\t", header=0, index_col=0).sort_index(axis=1), os.path.splitext(os.path.basename(x))[0])
        for x in counts_table_files
    ]

    # Common genes
    common_genes = set(datasets[0][0].columns)
    all_genes = set(datasets[0][0].columns)
    for dataset, _ in datasets[1:]:
        common_genes = common_genes.intersection(set(dataset.columns))
        all_genes = all_genes.union(set(dataset.columns))

    # Compute violin plots data
    violin_data_reads = []
    violin_data_genes = []
    violin_data_pos = []
    genes_counts = pd.DataFrame(index=pd.Index([x[1] for x in datasets]), columns=pd.Index(all_genes))
    genes_counts.fillna(0, inplace=True)
    for i, (dataset, name) in enumerate(datasets):
        violin_data_reads.append(dataset.sum(axis=1).to_list())
        violin_data_genes.append((dataset > 0).sum(axis=1).to_list())
        violin_data_pos.append(i + 1)
        genes_counts.loc[name, dataset.columns] = dataset.sum(axis=0)
    genes_counts.fillna(0, inplace=True)

    # Create the violin plots
    create_violin_plot(
        violin_data_reads,
        violin_data_pos,
        "Total reads",
        os.path.join(outdir, "violin_plot_reads.pdf"),
    )
    create_violin_plot(
        violin_data_genes,
        violin_data_pos,
        "Total genes",
        os.path.join(outdir, "violin_plot_genes.pdf"),
    )

    # Compute and plot PCA (sum gene counts)
    decomp_model = PCA(n_components=2, whiten=True, copy=True)
    reduced_data = decomp_model.fit_transform(np.log1p(genes_counts))
    create_pca_plot(
        reduced_data, genes_counts.index.to_list(), "PCA (sum gene counts)", os.path.join(outdir, "pca.pdf")
    )

    # Measure percentage of genes and correlations
    genes_similarities = pd.DataFrame(index=counts_table_files, columns=counts_table_files)
    genes_similarities.fillna(0, inplace=True)

    genes_correlations = pd.DataFrame(index=counts_table_files, columns=counts_table_files)
    genes_correlations.fillna(0, inplace=True)

    # Compute and create gene correlation plots
    n_col = len(counts_table_files)
    n_row = len(counts_table_files)
    plt.rcParams.update({"font.size": 6})
    fig, ax = plt.subplots(n_col, n_row, sharex="col", sharey="row", figsize=(3 * n_col, 3 * n_row))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i, (d1, n1) in enumerate(datasets):
        genes_1 = set(d1.columns)
        n_genes_1 = float(len(genes_1))
        common_to_all = (n_genes_1 - len(genes_1 - common_genes)) / n_genes_1
        print("{} shares {} genes with all the rest".format(n1, common_to_all))
        for j, (d2, n2) in enumerate(datasets):
            genes_2 = set(d2.columns)
            common_to_d2 = (n_genes_1 - len(genes_1 - genes_2)) / n_genes_1
            common_d1_d2 = genes_1.intersection(genes_2)
            sum_counts_1 = d1.loc[:, common_d1_d2].sum(axis=0)  # type: ignore
            sum_counts_1 = np.log(sum_counts_1) if use_log else sum_counts_1
            sum_counts_2 = d2.loc[:, common_d1_d2].sum(axis=0)  # type: ignore
            sum_counts_2 = np.log(sum_counts_2) if use_log else sum_counts_2
            genes_similarities.loc[n1, n2] = common_to_d2
            genes_correlations.loc[n1, n2] = pearsonr(sum_counts_1, sum_counts_2)[0]
            ax[i, j].scatter(sum_counts_1, sum_counts_2, s=5, c="blue", edgecolor="none")
            ax[i, j].set_xlabel(n1)
            ax[i, j].set_ylabel(n2)

    fig.savefig(os.path.join(outdir, "gene_correlations.png"), format="png", dpi=180)
    genes_similarities.to_csv(os.path.join(outdir, "gene_similarities.tsv"), sep="\t")
    genes_correlations.to_csv(os.path.join(outdir, "gene_correlations.tsv"), sep="\t")

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "counts_matrix_files", nargs="+", help="One or more matrices with gene counts (genes as columns) in TSV format"
    )
    parser.add_argument("--outdir", default=None, help="Path to the output directory")
    parser.add_argument(
        "--use-log-scale", action="store_true", default=False, help="Convert counts to log space for the correlation"
    )
    args = parser.parse_args()

    sys.exit(main(args.counts_matrix_files, args.outdir, args.use_log_scale))
