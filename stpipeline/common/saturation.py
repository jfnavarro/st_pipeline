"""
This module contains routines
to compute saturation points on a
set of annotated reads in BAM/SAM format
"""
import pysam
import random
import math
import os
from collections import defaultdict
import logging
from stpipeline.common.dataset import createDataset
from stpipeline.common.utils import safe_remove
from typing import List, Dict, Optional, Tuple

logger = logging.getLogger("STPipeline")


def compute_saturation(
    nreads: int,
    annotated_reads: str,
    gff_filename: str,
    umi_cluster_algorithm: str,
    umi_allowed_mismatches: int,
    umi_counting_offset: int,
    disable_umi: bool,
    expName: str,
    temp_folder: str,
    saturation_points: Optional[List[int]] = None,
) -> None:
    """
    Computes saturation points from annotated reads and logs the results.

    Args:
        nreads: Total number of reads in the annotated_reads file.
        annotated_reads: Path to a BAM file with the annotated reads.
        gff_filename: Path to the GFF file.
        umi_cluster_algorithm: Clustering algorithm for UMIs.
        umi_allowed_mismatches: Number of allowed mismatches.
        umi_counting_offset: Number of bases allowed as offset when counting UMIs.
        disable_umi: If True, disables UMI filtering.
        expName: Experiment name for logging and file organization.
        temp_folder: Path to temporary folder for intermediate files.
        saturation_points: List of saturation points to be used.

    Raises:
        RuntimeError: If the input file is missing or invalid.
    """
    if not os.path.isfile(annotated_reads):
        msg = f"Error, input file not present: {annotated_reads}"
        logger.error(msg)
        raise RuntimeError(msg)

    saturation_points = _determine_saturation_points(nreads, saturation_points)
    files, file_names, subsampling = _generate_subsamples(nreads, annotated_reads, saturation_points, temp_folder)

    _write_subsamples_to_files(files, subsampling, annotated_reads, saturation_points)

    results = _compute_saturation_metrics(
        file_names,
        saturation_points,
        gff_filename,
        umi_cluster_algorithm,
        umi_allowed_mismatches,
        umi_counting_offset,
        disable_umi,
        temp_folder,
        expName,
    )

    _cleanup_files(file_names)

    # TODO write this to a CSV file
    logger.info("Saturation points: %s", ", ".join(map(str, saturation_points)))
    logger.info("Reads per saturation point: %s", ", ".join(map(str, results["reads"])))
    logger.info("Genes per saturation point: %s", ", ".join(map(str, results["genes"])))
    logger.info("Average genes/spot per saturation point: %s", ", ".join(map(str, results["avg_genes"])))
    logger.info("Average reads/spot per saturation point: %s", ", ".join(map(str, results["avg_reads"])))


def _determine_saturation_points(nreads: int, saturation_points: Optional[List[int]]) -> List[int]:
    """
    Creates a list of saturations points using the total number of reads. If the points are
    provided they will be filtered if they are not then they will be extracted from an
    exponential distribution.
    """
    if saturation_points is None:
        saturation_points = [int(math.floor(1e3 + (math.exp(x) * 1e3))) for x in range(15)]
    points = [p for p in sorted(saturation_points) if p < nreads]
    if not points:
        msg = "All saturation points are greater than the number of reads."
        logger.error(msg)
        raise RuntimeError(msg)
    return points


def _generate_subsamples(
    nreads: int, annotated_reads: str, saturation_points: List[int], temp_folder: Optional[str]
) -> Tuple[Dict[int, pysam.AlignmentFile], Dict[int, str], Dict[int, List[int]]]:
    """
    Creates the saturation files for each saturation point.
    """
    file_ext = os.path.splitext(annotated_reads)[1].lower()
    flag_read, flag_write = ("rb", "wb") if file_ext != ".sam" else ("r", "wh")

    annotated_sam = pysam.AlignmentFile(annotated_reads, flag_read)  # type: ignore
    files, file_names, subsampling = {}, {}, {}

    for spoint in saturation_points:
        file_name = f"subsample_{spoint}{file_ext}"
        file_name = os.path.join(temp_folder, file_name) if temp_folder else file_name
        files[spoint] = pysam.AlignmentFile(file_name, flag_write, template=annotated_sam)  # type: ignore
        file_names[spoint] = file_name
        indices = list(range(nreads))
        random.shuffle(indices)
        subsampling[spoint] = sorted(indices[:spoint])

    annotated_sam.close()
    return files, file_names, subsampling


def _write_subsamples_to_files(
    files: Dict[int, pysam.AlignmentFile],
    subsampling: Dict[int, List[int]],
    annotated_reads: str,
    saturation_points: List[int],
) -> None:
    """
    Fill the content of each saturation point (file) using the annodated reads.
    """
    annotated_sam = pysam.AlignmentFile(annotated_reads, "rb")
    index = 0
    sub_indexes = defaultdict(int)  # type: ignore

    for read in annotated_sam.fetch(until_eof=True):
        for spoint in saturation_points:
            sub_index = sub_indexes[spoint]
            if sub_index < len(subsampling[spoint]) and subsampling[spoint][sub_index] == index:
                files[spoint].write(read)
                sub_indexes[spoint] += 1
        index += 1

    annotated_sam.close()
    for file_sam in files.values():
        file_sam.close()


def _compute_saturation_metrics(
    file_names: Dict[int, str],
    saturation_points: List[int],
    gff_filename: str,
    umi_cluster_algorithm: str,
    umi_allowed_mismatches: int,
    umi_counting_offset: int,
    disable_umi: bool,
    temp_folder: str,
    expName: str,
) -> Dict[str, List[int]]:
    """
    Generates the dataset for each saturation point (file) and fetch stats.
    """
    results = {"reads": [], "genes": [], "avg_genes": [], "avg_reads": []}  # type: ignore

    for spoint in saturation_points:
        input_file = file_names[spoint]
        try:
            stats = createDataset(
                input_file,
                temp_folder,
                gff_filename,
                umi_cluster_algorithm,
                umi_allowed_mismatches,
                umi_counting_offset,
                disable_umi,
                expName,
                verbose=False,
            )
        except Exception as e:
            logger.error("Error computing saturation curve: createDataset failed.")
            raise e

        results["reads"].append(stats["reads_after_duplicates_removal"])
        results["genes"].append(stats["genes_found"])
        results["avg_genes"].append(stats["average_genes_feature"])
        results["avg_reads"].append(stats["average_reads_feature"])

    return results


def _cleanup_files(file_names: Dict[int, str]) -> None:
    """
    Remove the files created during the saturation curve.
    """
    for file_name in file_names.values():
        safe_remove(file_name)
