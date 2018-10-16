#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" 
Script that creates multiple QC plots and stats
from multiple Spatial Transcriptomics datasets. 
Useful when one have several consecutive sections
or sections from the same model. 

The tool generates:

- Violin plots (genes and counts)
- Genes shared % pair-wise matrix
- Correlation pair-wise matrix
- Correlation pair-wise scatter plots
- PCA plot (one dot per dataset)

@Author Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>
"""

import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
from sklearn.decomposition import PCA
import os
import sys

color_map = ["red", "green", "blue", "orange", "cyan", "yellow", "orchid", 
             "saddlebrown", "darkcyan", "gray", "darkred", "darkgreen", "darkblue", 
             "antiquewhite", "bisque", "black", "slategray", "gold", "floralwhite",
             "aliceblue", "plum", "cadetblue", "coral", "olive", "khaki", "lightsalmon"]

def main(counts_table_files, outdir, use_log, n_col):
    
    if outdir is None or not os.path.isdir(outdir): 
        outdir = os.getcwd()
    outdir = os.path.abspath(outdir)

    print("Output directory {}".format(outdir))
    print("Input datasets {}".format(" ".join(counts_table_files))) 
    
    # Parse datasets and sort them by column
    datasets = [(pd.read_table(x, 
                               sep="\t", 
                               header=0, 
                               index_col=0).sort_index(axis=1), x) for x in counts_table_files]
    
    # Common genes
    common_genes = set(datasets[0][0].columns)
    for (dataset,_) in datasets[1:]:
        common_genes = common_genes.intersection(dataset.columns)
  
    # Compute violin plots
    violin_data_reads = list()
    violin_data_genes = list()
    violin_data_names = list()
    violin_data_pos = list()
    genes_counts = pd.DataFrame(index=counts_table_files, columns=common_genes)
    genes_counts.fillna(0)
    for i,(dataset, name) in enumerate(datasets):
        violin_data_reads.append(dataset.sum(axis=1).values)
        violin_data_genes.append((dataset > 0).sum(axis=1).values)
        violin_data_names.append(name)
        violin_data_pos.append(i)
        genes_counts.loc[name,:] = dataset.loc[:,common_genes].sum(axis=0)
        
    fig, axes = plt.subplots(nrows=2, ncols=1)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    axes[0].violinplot(violin_data_reads, violin_data_pos, points=50, widths=0.5,
                       showmeans=True, showextrema=True, showmedians=True)
    axes[0].set_title('Violin plot (reads)', fontsize=6)
    axes[1].violinplot(violin_data_genes, violin_data_pos, points=50, widths=0.5,
                       showmeans=True, showextrema=True, showmedians=True)
    axes[1].set_title('Violin plot (genes)', fontsize=6)
    for ax in axes.flatten():
        ax.set_xticklabels(violin_data_names)
    fig.savefig(os.path.join(outdir,"violin_plots.pdf"), format='pdf', dpi=90)
          
    decomp_model = PCA(n_components=2, whiten=True, copy=True)
    reduced_data = decomp_model.fit_transform(genes_counts)
    
    fig, ax = plt.subplots()
    class_colours = [color_map[i] for i in range(0,len(reduced_data))]
    ax.scatter(reduced_data[:,0], reduced_data[:,1], s=10, c=class_colours,
               edgecolor="none")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA scatter (sum. gene counts)")
    recs = list()
    for i in range(0,len(class_colours)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[i]))
    ax.legend(recs,genes_counts.index,loc=4)
    fig.savefig(os.path.join(outdir,"pca.pdf"), format='pdf', dpi=90)
 
    # Measssure percentage of genes and correlations
    genes_similarities = pd.DataFrame(index=counts_table_files, columns=counts_table_files)
    genes_similarities.fillna(0)
    
    genes_correlations = pd.DataFrame(index=counts_table_files, columns=counts_table_files)
    genes_correlations.fillna(0)

    n_plots = len(counts_table_files) * len(counts_table_files)
    n_col = min(n_col, n_plots) 
    n_row = max(int(n_plots / n_col), 1) 
    plt.rcParams.update({'font.size': 6})
    fig, ax = plt.subplots(n_col, n_row, sharex='col', sharey='row')
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i,(d1,n1) in enumerate(datasets):
        genes_1 = set(d1.columns)
        n_genes_1 = float(len(genes_1))
        common_to_all = (n_genes_1 - len(genes_1 - common_genes)) / n_genes_1
        print("{} shares {}% genes with all the rest".format(n1, common_to_all))
        for j,(d2,n2) in enumerate(datasets):
            genes_2 = set(d2.columns)
            common_to_d2 = (n_genes_1 - len(genes_1 - genes_2)) / n_genes_1
            common_d1_d2 = genes_1.intersection(genes_2)
            sum_counts_1 = d1.loc[:,common_d1_d2].sum(axis=0)
            sum_counts_1 = np.log(sum_counts_1) if use_log else sum_counts_1
            sum_counts_2 = d2.loc[:,common_d1_d2].sum(axis=0)
            sum_counts_2 = np.log(sum_counts_2) if use_log else sum_counts_2
            genes_similarities.loc[n1,n2] = common_to_d2
            genes_correlations.loc[n1,n2] = pearsonr(sum_counts_1, sum_counts_2)[0]
            ax[i,j].scatter(sum_counts_1, sum_counts_2, s=5, c="blue", edgecolor="none")
            ax[i,j].set_xlabel(n1)
            ax[i,j].set_ylabel(n2)
    
    fig.savefig(os.path.join(outdir,"gene_correlations.pdf"), format='pdf', dpi=90)
    genes_similarities.to_csv(os.path.join(outdir,"gene_similarities.tsv"), sep='\t')
    genes_correlations.to_csv(os.path.join(outdir,"gene_correlations.tsv"), sep='\t')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--counts-table-files", required=True, nargs='+', type=str,
                        help="One or more matrices with gene counts per feature/spot (genes as columns)")
    parser.add_argument("--outdir", default=None, help="Path to output dir")
    parser.add_argument("--use-log-scale", action="store_true", default=False, 
                        help="Plot expression in log space (log2)")
    parser.add_argument("--num-columns", default=1, type=int, metavar="[INT]",
                        help="The number of columns when using --joint-plot (default: %(default)s)")
    args = parser.parse_args()

    main(args.counts_table_files,
         args.outdir,
         args.use_log_scale,
         args.num_columns)
