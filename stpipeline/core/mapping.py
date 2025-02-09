"""
This module contains functions related to sequence alignment and barcode
demultiplexing in the ST pipeline
"""

import logging
import os
import shutil
import subprocess
from typing import List, Optional

from stpipeline.common.utils import file_ok

logger = logging.getLogger("STPipeline")


def alignReads(
    reverse_reads: str,
    ref_map: str,
    outputFile: str,
    annotation: Optional[str],
    outputFolder: str,
    trimReverse: int,
    invTrimReverse: int,
    cores: int,
    min_intron_size: int,
    max_intron_size: int,
    disable_multimap: bool,
    diable_softclipping: bool,
    twopassMode: bool,
    min_length: int,
    include_non_mapped: bool,
    star_genome_loading: str,
    star_sort_mem_limit: int,
) -> int:
    """
    Perform sequence alignment using STAR.

    Args:
        reverse_reads: File containing reverse reads in BAM format.
        ref_map: Path to the STAR genome/transcriptome index.
        outputFile: Name of the SAM/BAM output file for alignments.
        annotation: GTF annotation file path (optional).
        outputFolder: Path to the output folder.
        trimReverse: Number of bases to trim from the 5' end of reverse reads.
        invTrimReverse: Number of bases to trim from the 3' end of reverse reads.
        cores: Number of cores for alignment.
        min_intron_size: Minimum allowed intron size for splice junctions.
        max_intron_size: Maximum allowed intron size for splice junctions.
        disable_multimap: If True, disallow multiple alignments.
        diable_softclipping: If True, disable local alignment.
        twopassMode: If True, enable 2-pass mode.
        min_length: Minimum allowed read length after trimming.
        include_non_mapped: If True, include unaligned reads in the output.
        star_genome_loading: Type of genome loading for STAR.
        star_sort_mem_limit: Memory limit for BAM sorting by STAR.

    Returns:
        The total number of reads mapped.

    Raises:
        RuntimeError: If input files are missing or output file creation fails.
        ValueError: For invalid input arguments.
        OSError: If STAR executable is not found.
    """

    if not file_ok(reverse_reads):
        error = f"Error mapping with STAR, input file not present {reverse_reads}\n"
        logger.error(error)
        raise RuntimeError(error)

    # STAR has predefined output names for the files
    tmpOutputFile = "Aligned.sortedByCoord.out.bam"
    log_final = "Log.final.out"

    if outputFolder is not None:
        tmpOutputFile = os.path.join(outputFolder, tmpOutputFile)
        log_final = os.path.join(outputFolder, log_final)

    multi_map_number = 1 if disable_multimap else 20  # 10 is the STAR default
    alignment_mode = "EndToEnd" if diable_softclipping else "Local"

    # Prepare STAR command arguments
    args = [
        "STAR",
        "--genomeDir",
        ref_map,
        "--readFilesIn",
        reverse_reads,
        "--outFileNamePrefix",
        f"{outputFolder}{os.sep}",
        "--clip3pNbases",
        str(invTrimReverse),
        "--clip5pNbases",
        str(trimReverse),
        "--runThreadN",
        str(max(cores, 1)),
        "--outFilterType",
        "Normal",
        "--outSAMtype",
        "BAM",
        "SortedByCoordinate",
        "--alignEndsType",
        alignment_mode,
        "--outSAMorder",
        "Paired",
        "--outSAMprimaryFlag",
        "OneBestScore",
        "--outFilterMultimapNmax",
        str(multi_map_number),
        "--alignIntronMin",
        str(min_intron_size),
        "--alignIntronMax",
        str(max_intron_size),
        "--outFilterMatchNmin",
        str(min_length),
        "--genomeLoad",
        star_genome_loading,
        "--limitBAMsortRAM",
        str(star_sort_mem_limit),
        "--readFilesType",
        "SAM",
        "SE",  # Input in BAM format
        "--readFilesCommand",
        "samtools",
        "view",
        "-h",
    ]

    if twopassMode:
        args += ["--twopassMode", "Basic"]

    if annotation:
        args += ["--sjdbGTFfile", annotation]

    if include_non_mapped:
        args += ["--outSAMunmapped", "Within"]
    else:
        args += ["--outSAMunmapped", "None"]

    try:
        logger.debug(f"STAR mapping, running command: {' '.join(args)}")
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True, shell=False)
        _, errmsg = proc.communicate()

        if proc.returncode != 0:
            logger.error(f"Error mapping with STAR: {errmsg.decode()}")
            raise RuntimeError(f"STAR mapping failed with error: {errmsg.decode()}")
    except OSError as e:
        logger.error("Error mapping with STAR: Executable not found.")
        raise e

    if not file_ok(tmpOutputFile):
        error = f"Error mapping with STAR: Output file not present {tmpOutputFile}"
        logger.error(error)
        raise RuntimeError(error)

    # Rename output file
    shutil.move(tmpOutputFile, outputFile)

    uniquely_mapped = 0
    multiple_mapped = 0
    if not file_ok(log_final):
        logger.warning("Log output file from STAR is not present")
    else:
        logger.info("Mapping stats:")
        logger.info("Mapping stats are computed from all the pair reads present in the raw files")
        with open(log_final, "r") as star_log:
            for line in star_log:
                if "Uniquely mapped reads number" in line:
                    uniquely_mapped = int(line.strip().split()[-1])
                    logger.info(line.strip())
                elif "Number of reads mapped to multiple loci" in line:
                    multiple_mapped += int(line.strip().split()[-1])
                    logger.info(line.strip())
                elif "% of reads mapped to multiple loci" in line or "% of reads unmapped: too short" in line:
                    logger.info(line.strip())
            logger.info("Total mapped reads: {}".format(uniquely_mapped + multiple_mapped))

    return uniquely_mapped + multiple_mapped


def barcodeDemultiplexing(
    reads: str,
    idFile: str,
    mismatches: int,
    kmer: int,
    over_hang: int,
    taggd_metric: str,
    taggd_multiple_hits_keep_one: bool,
    taggd_trim_sequences: Optional[List[int]],
    cores: int,
    outputFilePrefix: str,
    keep_discarded_files: bool = False,
    taggd_chunk_size: int = 10000,
) -> int:
    """
    Perform demultiplexing using Taggd.

    Args:
        reads: Path to the FASTQ/BAM file containing barcoded reads.
        idFile: Path to the tab-delimited barcode file (BARCODE - X - Y).
        mismatches: Allowed mismatches for barcode matching.
        kmer: K-mer length for barcode search.
        over_hang: Allowed flanking bases around barcodes.
        taggd_metric: Distance metric algorithm ('Hamming', 'Levenshtein', 'Subglobal').
        taggd_multiple_hits_keep_one: If True, keep one random hit for multiple candidates.
        taggd_trim_sequences: List of coordinates to trim in the barcode (optional).
        cores: Number of subprocesses for Taggd.
        outputFilePrefix: Prefix for output files.
        keep_discarded_files: If True, generate files for unmatched reads.
        taggd_chunk_size: Number of reads to process in each thread.

    Returns:
        The total number of reads demultiplexed.

    Raises:
        RuntimeError: If input files are missing or output file creation fails.
        ValueError: For invalid input arguments.
        OSError: If Taggd executable is not found.
    """

    if not os.path.isfile(reads):
        error = f"Error, input file not present {reads}"
        logger.error(error)
        raise RuntimeError(error)

    # Taggd options
    # --metric (subglobal (default) , Levenshtein or Hamming)
    # --slider-increment (space between kmer searches, 0 is default = kmer length)
    # --seed
    # --overhang additional flanking bases around read barcode to allow
    # --estimate-min-edit-distance is set estimate the min edit distance among true barcodes
    # --no-offset-speedup turns off speed up, it might yield more hits (exactly as findIndexes)
    # --homopolymer-filter if set excludes reads where the barcode
    # --chunk-size number of reads to use in each thread for parallel processing
    #  contains a homolopymer of the given length (0 no filter), default 8

    args = [
        "taggd_demultiplex",
        "--max-edit-distance",
        str(mismatches),
        "--k",
        str(kmer),
        "--barcode-tag",
        "B0",  # if input is BAM we tell taggd which tag contains the barcode
        "--homopolymer-filter",
        "0",
        "--subprocesses",
        str(cores),
        "--metric",
        taggd_metric,
        "--chunk-size",
        str(taggd_chunk_size),
        "--overhang",
        str(over_hang if taggd_metric != "Hamming" else 0),
    ]

    if taggd_trim_sequences:
        args += ["--trim-sequences"] + list(map(str, taggd_trim_sequences))

    if taggd_multiple_hits_keep_one:
        args.append("--multiple-hits-keep-one")

    if not keep_discarded_files:
        args += ["--no-unmatched-output", "--no-ambiguous-output", "--no-results-output"]

    args += [idFile, reads, outputFilePrefix]

    try:
        logger.debug(f"Taggd demultiplexing, running command: {' '.join(args)}")
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True, shell=False)
        stdout, errmsg = proc.communicate()

        if proc.returncode != 0:
            logger.error(f"Error demultiplexing with Taggd: {errmsg.decode()}")
            raise RuntimeError(f"Taggd demultiplexing failed with error: {errmsg.decode()}")
    except OSError as e:
        logger.error("Error demultiplexing with Taggd: Executable not found.")
        raise e

    outputFile = f"{outputFilePrefix}_matched{os.path.splitext(reads)[1].lower()}"
    if not file_ok(outputFile):
        error = f"Error demultiplexing with Taggd: Output file not present {outputFile}"
        logger.error(error)
        raise RuntimeError(error)

    # Write log file and collect stats
    outputLog = f"{outputFilePrefix}_log.txt"
    reads_after_demultiplexing = 0
    with open(outputLog, "w") as fwrite:
        for line in stdout.decode().splitlines():
            fwrite.write(line + "\n")
            tokens = [
                "Total reads:",
                "Perfect Matches:",
                "Imperfect Matches",
                "Ambiguous matches:",
                "Non-unique ambiguous matches:",
                "Unmatched:",
            ]
            if any(x in line for x in tokens):
                logger.info(line.strip())
            if "Total reads written:" in line:
                logger.info(line.strip())
                reads_after_demultiplexing = int(line.strip().split()[-1])

    return reads_after_demultiplexing
