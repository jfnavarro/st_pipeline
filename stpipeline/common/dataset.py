"""
This module contains routines to create
a ST dataset and some statistics. The dataset
will contain several files with the ST data in different
formats
"""

import logging
import os
import random
from collections import defaultdict
from typing import Any, Callable, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from stpipeline.common.clustering import dedup_adj, dedup_dir_adj, dedup_hierarchical
from stpipeline.common.transcript import Transcript
from stpipeline.common.unique_events_parser import parse_unique_events

logger = logging.getLogger("STPipeline")


def compute_unique_umis(
    transcripts: List[Transcript],
    umi_counting_offset: int,
    umi_allowed_mismatches: int,
    group_umi_func: Callable[[List[str], int], List[str]],
) -> List[Transcript]:
    """
    Computes unique UMIs from a list of transcripts, grouping them by genomic coordinates and strand.

    The function groups UMIs by strand and start position (with an offset), clusters the UMIs within
    each group using a grouping function, and computes unique UMIs based on Hamming distance.

    Args:
        transcripts: A list of transcript data. Each transcript should include positional and strand information.
        umi_counting_offset: The maximum offset allowed when grouping transcripts.
        umi_allowed_mismatches: The maximum allowed mismatches for UMIs in a group.
        group_umi_func: A function to group UMIs, accepting a list of UMIs and mismatch threshold.

    Returns:
        A list of unique transcripts, one for each unique UMI.
    """
    # Sort transcripts by strand and start position
    sorted_transcripts = sorted(transcripts, key=lambda x: (x.strand, x.start))
    grouped_transcripts = defaultdict(list)
    unique_transcripts = []
    num_transcripts = len(transcripts)
    for i in range(num_transcripts - 1):
        current = sorted_transcripts[i]
        nextone = sorted_transcripts[i + 1]
        grouped_transcripts[current.umi].append(current)
        if abs(current.start - nextone.start) > umi_counting_offset or current.strand != nextone.strand:
            # A new group has been reached (strand, start-pos, offset)
            unique_umis = group_umi_func(list(grouped_transcripts.keys()), umi_allowed_mismatches)
            unique_transcripts += [random.choice(grouped_transcripts[u_umi]) for u_umi in unique_umis]
            grouped_transcripts = defaultdict(list)

    # Process the last group
    lastone = sorted_transcripts[num_transcripts - 1]
    grouped_transcripts[lastone.umi].append(lastone)
    unique_umis = group_umi_func(list(grouped_transcripts.keys()), umi_allowed_mismatches)
    unique_transcripts += [random.choice(grouped_transcripts[u_umi]) for u_umi in unique_umis]

    return unique_transcripts


def createDataset(
    input_file: str,
    output_folder: str,
    gff_filename: Optional[str] = None,
    umi_cluster_algorithm: str = "AdjacentBi",
    umi_allowed_mismatches: int = 1,
    umi_counting_offset: int = 250,
    disable_umi: bool = False,
    output_template: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Parses an annotated and demultiplexed BAM file with reads, it groups them by gene-barcode
    to count unique transcripts, and removes duplicates using UMIs. It creates a data frame of unique
    transcripts by spots and genes. It returns a dictionary of basic stats.

    Args:
        input_file: Path to the BAM file containing annotated-demultiplexed records.
        output_folder: Directory for output files.
        gff_filename: Annotation reference file. Defaults to None.
        umi_cluster_algorithm : Algorithm for clustering UMIs. Defaults to "hierarchical".
        umi_allowed_mismatches: Allowed mismatches for UMI deduplication. Defaults to 1.
        umi_counting_offset: Offset for grouping transcripts by position. Defaults to 250.
        disable_umi: Disables UMI deduplication if True. Defaults to False.
        output_template: Template for output file names. Defaults to None.
        verbose: Enables verbose logging if True. Defaults to True.

    Returns:
        A dictionary with basic stats.

    Raises:
        RuntimeError: If input file is missing or errors occur during processing.
    """
    if not os.path.isfile(input_file):
        error = f"Error creating dataset, input file not present {input_file}"
        logger.error(error)
        raise RuntimeError(error)

    # Set default filenames for the output files
    filenameDataFrame = f"{output_template}_stdata.tsv" if output_template else "stdata.tsv"
    filenameReadsBED = f"{output_template}_reads.bed" if output_template else "reads.bed"

    # Initialize counters
    total_record = 0
    discarded_reads = 0

    # Obtain the appropriate UMI clustering function
    group_umi_func = {"hierarchical": dedup_hierarchical, "Adjacent": dedup_adj, "AdjacentBi": dedup_dir_adj}.get(
        umi_cluster_algorithm
    )

    if not group_umi_func:
        error = f"Error creating dataset. Incorrect clustering algorithm {umi_cluster_algorithm}"
        logger.error(error)
        raise RuntimeError(error)

    # Containers to store data for creating the DataFrame
    list_row_values = []
    list_indexes = []

    # Parse unique events to generate the unique counts (reads) dataframe and a BED file
    unique_events = parse_unique_events(input_file, gff_filename)
    with open(os.path.join(output_folder, filenameReadsBED), "w") as reads_handler:
        # unique_events is a list of tuples (gene, spots)
        # where gene is a str and spots is a dictionary of transcripts per spot [spot] -> List[Transcript]
        # this loop is to make the transcripts unique (deduplicate them)
        for gene, spots in unique_events:
            unique_transcripts_by_spot = {}
            for spot_coordinates, transcripts in spots.items():
                x, y = spot_coordinates
                transcripts_count = len(transcripts)

                # Compute unique transcripts based on UMI, strand, and start position
                unique_transcripts = (
                    compute_unique_umis(transcripts, umi_counting_offset, umi_allowed_mismatches, group_umi_func)  # type: ignore
                    if not disable_umi
                    else transcripts
                )

                unique_transcripts_count = len(unique_transcripts)
                assert 0 < unique_transcripts_count <= transcripts_count
                discarded_reads += transcripts_count - unique_transcripts_count
                unique_transcripts_by_spot[f"{x}x{y}"] = unique_transcripts_count

                # Write unique transcripts to the BED file
                for t in unique_transcripts:
                    reads_handler.write(
                        f"{t.chrom}\t{t.start}\t{t.end}\t{t.clear_name}\t{t.mapping_quality}\t{t.strand}\t{gene}\t{x}\t{y}\n"
                    )

                total_record += 1

            # Add data for the DataFrame
            list_indexes.append(gene)
            list_row_values.append(unique_transcripts_by_spot)

    if total_record == 0:
        error = "Error creating dataset, input file did not contain any transcript"
        logger.error(error)
        raise RuntimeError(error)

    # Create the counts DataFrame
    counts_table = pd.DataFrame(list_row_values, index=list_indexes).fillna(0).T

    # Write the counts DataFrame to a TSV file
    counts_table.to_csv(os.path.join(output_folder, filenameDataFrame), sep="\t", na_rep=0)  # type: ignore

    # Compute statistics for the dataset
    total_spots, number_genes = counts_table.shape
    total_reads = np.sum(counts_table.values, dtype=np.int32)
    aggregated_spot_counts = counts_table.sum(axis=1)
    aggregated_gene_counts = (counts_table != 0).sum(axis=1)

    stats_dict: Dict[str, Union[int, float]] = {}
    stats_dict["max_genes_feature"] = int(aggregated_gene_counts.max())
    stats_dict["min_genes_feature"] = int(aggregated_gene_counts.min())
    stats_dict["max_reads_feature"] = float(aggregated_spot_counts.max())
    stats_dict["min_reads_feature"] = float(aggregated_spot_counts.min())
    stats_dict["average_reads_feature"] = float(np.mean(aggregated_spot_counts))
    stats_dict["average_genes_feature"] = float(np.mean(aggregated_gene_counts))
    stats_dict["std_reads_feature"] = float(np.std(aggregated_spot_counts))
    stats_dict["std_genes_feature"] = float(np.std(aggregated_gene_counts))
    stats_dict["reads_after_duplicates_removal"] = int(total_reads)
    stats_dict["barcodes_found"] = total_spots
    stats_dict["genes_found"] = number_genes
    stats_dict["duplicates_found"] = discarded_reads

    # Log statistics if verbose mode is enabled
    if verbose:
        logger.info(f"Number of reads present: {total_reads}")
        logger.info(f"Number of unique events (gene-spot) present: {total_record}")
        logger.info(f"Number of unique genes present: {number_genes}")
        logger.info(f"Max number of genes over all spots: {stats_dict['max_genes_feature']}")
        logger.info(f"Min number of genes over all spots: {stats_dict['min_genes_feature']}")
        logger.info(f"Max number of reads over all spots: {stats_dict['max_reads_feature']}")
        logger.info(f"Min number of reads over all spots: {stats_dict['min_reads_feature']}")
        logger.info(f"Average number genes per spot: {stats_dict['average_genes_feature']}")
        logger.info(f"Average number reads per spot: {stats_dict['average_reads_feature']}")
        logger.info(f"Std. number genes per spot: {stats_dict['std_genes_feature']}")
        logger.info(f"Std. number reads per spot: {stats_dict['std_reads_feature']}")
        logger.info(f"Number of discarded reads (possible duplicates): {discarded_reads}")

    return stats_dict
