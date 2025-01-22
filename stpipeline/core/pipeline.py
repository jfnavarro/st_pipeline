# coding=utf-8
"""
This is the main object for the ST pipeline.
It contains the Pipeline object that has methods
to parse the input parameters, generate the input parameters,
do sanity check and ultimately run the pipeline.
"""

import argparse
import inspect
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile

import pysam
from taggd.io.barcode_utils import read_barcode_file  # type: ignore

from stpipeline.common.dataset import createDataset
from stpipeline.common.filter import filter_input_data
from stpipeline.common.saturation import compute_saturation
from stpipeline.common.stats import Stats
from stpipeline.common.utils import (
    TimeStamper,
    get_htseq_count_version,
    get_star_version,
    get_taggd_count_version,
    safe_remove,
    which_program,
)
from stpipeline.core.annotation import annotateReads
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing
from stpipeline.version import version_number

FILENAMES = {
    "mapped": "mapped.bam",
    "annotated": "annotated.bam",
    "contaminated_clean": "contaminated_clean.bam",
    "demultiplexed_prefix": "demultiplexed",
    "demultiplexed_matched": "demultiplexed_matched.bam",
    "quality_trimmed_R2": "R2_quality_trimmed.bam",
}

FILENAMES_DISCARDED = {
    "mapped_discarded": "mapping_discarded.bam",
    "contaminated_discarded": "contaminated.bam",
    "demultiplexed_ambiguos": "demultiplexed_ambiguos.bam",
    "demultiplexed_unmatched": "demultiplexed_unmatched.bam",
    "demultiplexed_results": "demultiplexed_results.tsv",
    "quality_trimmed_discarded": "R2_quality_trimmed_discarded.fastq",
    "annotated_discarded": "annotated_discarded.bam",
}

logger = logging.getLogger("STPipeline")


class Pipeline:
    """
    This class contains all the ST pipeline
    attributes and a bunch of methods to parse
    the input parameters, do sanity check and
    run the pipeline steps.
    """

    def __init__(self) -> None:
        self.allowed_missed = 2
        self.allowed_kmer = 6
        self.overhang = 0
        self.min_length_trimming = 20
        self.trimming_rv = 0
        self.min_quality_trimming = 20
        self.clean = True
        self.barcode_start = 0
        self.threads = 8
        self.verbose = False
        self.ids = ""
        self.ref_map = ""
        self.ref_annotation = None
        self.expName = ""
        self.htseq_mode = "intersection-nonempty"
        self.htseq_no_ambiguous = False
        self.htseq_features = ["exon"]
        self.qual64 = False
        self.contaminant_index = None
        self.fastq_fw = ""
        self.fastq_rv = ""
        self.logfile = None
        self.output_folder = None
        self.temp_folder = None
        self.umi_allowed_mismatches = 1
        self.umi_start_position = 18
        self.umi_end_position = 27
        self.keep_discarded_files = False
        self.remove_polyA_distance = 10
        self.remove_polyT_distance = 10
        self.remove_polyG_distance = 10
        self.remove_polyC_distance = 10
        self.remove_polyN_distance = 10
        self.filter_AT_content = 90
        self.filter_GC_content = 90
        self.disable_clipping = False
        self.disable_multimap = False
        self.umi_cluster_algorithm = "AdjacentBi"
        self.min_intron_size = 1
        self.max_intron_size = 1
        self.umi_filter = False
        self.umi_filter_template = "WSNNWSNNV"
        self.compute_saturation = False
        self.include_non_annotated = False
        self.inverse_trimming_rv = 0
        self.two_pass_mode = False
        self.strandness = "yes"
        self.umi_quality_bases = 6
        self.umi_counting_offset = 250
        self.taggd_metric = "Subglobal"
        self.taggd_multiple_hits_keep_one = False
        self.taggd_trim_sequences = None
        self.taggd_chunk_size = 10000
        self.adaptor_missmatches = 0
        self.star_genome_loading = "NoSharedMemory"
        self.star_sort_mem_limit = 0
        self.disable_trimming = False
        self.disable_mapping = False
        self.disable_annotation = False
        self.disable_umi = False
        self.disable_barcode = False
        self.transcriptome = False
        self.saturation_points = None
        self.qa_stats = Stats()

    def clean_filenames(self) -> None:
        """
        Just makes sure to remove all temp files
        """
        if self.clean:
            for file_name in list(FILENAMES.values()):
                safe_remove(file_name)
        if not self.keep_discarded_files:
            for file_name in list(FILENAMES_DISCARDED.values()):
                safe_remove(file_name)
        if self.temp_folder is not None and os.path.isdir(self.temp_folder):
            star_temp1 = os.path.join(self.temp_folder, "_STARgenome")
            star_temp2 = os.path.join(self.temp_folder, "_STARpass1")
            if os.path.isdir(star_temp1):
                shutil.rmtree(star_temp1)
            if os.path.isdir(star_temp2):
                shutil.rmtree(star_temp2)
            if self.clean and not self.keep_discarded_files and self.temp_folder != self.output_folder:
                shutil.rmtree(self.temp_folder)

    def sanityCheck(self) -> None:
        """
        Performs some basic sanity checks on the input parameters
        """

        if self.ref_annotation is not None and (
            not os.path.isfile(self.ref_annotation)
            or (
                not self.ref_annotation.endswith(".gtf")
                and not self.ref_annotation.endswith(".gff3")
                and not self.ref_annotation.endswith(".gff")
            )
        ):
            error = f"Error parsing parameters.\nInvalid annotation file {self.ref_annotation}."
            logger.error(error)
            raise RuntimeError(error)

        if self.ref_annotation is None and not self.transcriptome:
            error = "Error, annotation file is missing but the transcriptome option is disabled."
            logger.error(error)
            raise RuntimeError(error)

        if self.ref_map is None and not self.disable_mapping:
            error = "Error, genome reference is missing but the disable_mapping option is disabled."
            logger.error(error)
            raise RuntimeError(error)

        if not os.path.isfile(self.fastq_fw) or not os.path.isfile(self.fastq_rv):
            error = f"Error parsing parameters.\nInvalid input files {self.fastq_fw} {self.fastq_rv}"
            logger.error(error)
            raise RuntimeError(error)

        if (
            not self.fastq_fw.endswith(".fastq")
            and not self.fastq_fw.endswith(".fq")
            and not self.fastq_fw.endswith(".gz")
            and not self.fastq_rv.endswith(".bz2")
        ) or (
            not self.fastq_rv.endswith(".fastq")
            and not self.fastq_rv.endswith(".fq")
            and not self.fastq_rv.endswith(".gz")
            and not self.fastq_rv.endswith(".bz2")
        ):
            error = f"Error parsing parameters.\nIncorrect format for input files {self.fastq_fw} {self.fastq_rv}"
            logger.error(error)
            raise RuntimeError(error)

        if not self.disable_barcode and not os.path.isfile(self.ids):
            error = f"Error parsing parameters.\nInvalid IDs file {self.ids}"
            logger.error(error)
            raise RuntimeError(error)

        if not self.disable_barcode and self.ids is None:
            error = "Error IDs file is missing but the option to disable the " "demultiplexing step is not activated\n"
            logger.error(error)
            raise RuntimeError(error)

        if self.saturation_points is not None and not self.compute_saturation:
            logger.warning("Saturation points are provided but the option" "to compute saturation is disabled.")

        if not self.disable_umi and self.umi_filter:
            # Check template validity
            regex = "[^ACGTURYSWKMBDHVN]"
            if re.search(regex, self.umi_filter_template) is not None:
                error = f"Error invalid UMI template given {self.umi_filter_template}."
                logger.error(error)
                raise RuntimeError(error)
            # Check template length
            if len(self.umi_filter_template) != (self.umi_end_position - self.umi_start_position):
                error = f"Error the UMI template given does not have the correct length {self.umi_filter_template}."
                logger.error(error)
                raise RuntimeError(error)
            # Convert the template into a reg-exp
            temp_reg_exp = ""
            for ele in self.umi_filter_template:
                if ele == "W":
                    temp_reg_exp += "[AT]"
                elif ele == "S":
                    temp_reg_exp += "[CG]"
                elif ele == "N":
                    temp_reg_exp += "[ATCG]"
                elif ele == "V":
                    temp_reg_exp += "[ACG]"
                elif ele == "A":
                    temp_reg_exp += "[A]"
                elif ele == "C":
                    temp_reg_exp += "[C]"
                elif ele == "G":
                    temp_reg_exp += "[G]"
                elif ele == "T":
                    temp_reg_exp += "[T]"
                elif ele == "U":
                    temp_reg_exp += "[U]"
                elif ele == "R":
                    temp_reg_exp += "[AG]"
                elif ele == "Y":
                    temp_reg_exp += "[CT]"
                elif ele == "K":
                    temp_reg_exp += "[GT]"
                elif ele == "M":
                    temp_reg_exp += "[AC]"
                elif ele == "B":
                    temp_reg_exp += "[CGT]"
                elif ele == "D":
                    temp_reg_exp += "[AGT]"
                elif ele == "H":
                    temp_reg_exp += "[ACT]"
            self.umi_filter_template = temp_reg_exp

        # Add checks for trimming parameters, demultiplex parameters and UMI parameters
        if self.allowed_missed > self.allowed_kmer and not self.disable_barcode:
            error = (
                "Error starting the pipeline.\n" "Taggd allowed mismatches is bigger or equal than the Taggd k-mer size"
            )
            logger.error(error)
            raise RuntimeError(error)

        if (
            self.umi_start_position < self.barcode_start < self.umi_end_position
            and not self.disable_barcode
            and not self.disable_umi
        ):
            error = (
                "Error starting the pipeline.\n"
                "The start position of the barcodes is between the UMIs start-end position"
            )
            logger.error(error)
            raise RuntimeError(error)

        if self.umi_allowed_mismatches > (self.umi_end_position - self.umi_start_position) and not self.disable_umi:
            error = "Error starting the pipeline.\n" "The allowed UMI mismatches is bigger than the UMI size"
            logger.error(error)
            raise RuntimeError(error)

        if self.taggd_chunk_size < 100:
            error = "Error starting the pipeline.\n" "The chunk size for the demultiplexing step is too small"
            logger.error(error)
            raise RuntimeError(error)

        # Test the presence of the required tools
        required_scripts = ["STAR"]
        unavailable_scripts = set()
        for script in required_scripts:
            if not which_program(script):
                unavailable_scripts.add(script)
        if len(unavailable_scripts) != 0:
            error = f"Error starting the pipeline. Required software not found: {' '.join(required_scripts)}"
            logger.error(error)
            raise RuntimeError(error)

    def createParameters(self, parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
        """
        Adds the pipeline"s parameters to a given
        Argparse object and returns it.
        """

        class readable_dir(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None) -> None:  # type: ignore
                prospective_dir = values
                if not os.path.isdir(prospective_dir):
                    raise argparse.ArgumentTypeError(f"{prospective_dir} is not a valid path")
                if os.access(prospective_dir, os.R_OK):
                    setattr(namespace, self.dest, prospective_dir)
                else:
                    raise argparse.ArgumentTypeError(f"{prospective_dir} is not a readable dir")

        parser.add_argument("fastq_files", nargs=2)
        parser.add_argument(
            "--ids",
            metavar="[FILE]",
            required=False,
            help="Path to the file containing the map of barcodes to the array coordinates",
        )
        parser.add_argument(
            "--ref-map",
            metavar="[FOLDER]",
            action=readable_dir,
            required=False,
            help="Path to the folder with the STAR index " "for the genome that you want to use as reference",
        )
        parser.add_argument(
            "--ref-annotation",
            metavar="[FILE]",
            required=False,
            help="Path to the reference annotation file "
            "(GTF or GFF format is required) to be used to annotated the mapped reads",
        )
        parser.add_argument(
            "--expName",
            type=str,
            metavar="[STRING]",
            required=True,
            help="Name of the dataset (The output files will prepend this name)",
        )
        parser.add_argument(
            "--contaminant-index",
            metavar="[FOLDER]",
            action=readable_dir,
            default=None,
            help="Path to the folder with a STAR index with a contaminant genome reference.\n"
            "Reads will be filtered using the specified genome and mapping reads will be discarded",
        )
        parser.add_argument(
            "--no-clean-up",
            action="store_false",
            default=True,
            help="Do not remove temporary/intermediary files (useful for debugging)",
        )
        parser.add_argument(
            "--verbose", action="store_true", default=False, help="Show extra information on the log file"
        )
        parser.add_argument(
            "--threads",
            default=8,
            metavar="[INT]",
            type=int,
            choices=range(1, 81),
            help="Number of threads to use (default: %(default)s)",
        )
        parser.add_argument(
            "--bin-path",
            metavar="[FOLDER]",
            action=readable_dir,
            default=None,
            help="Path to folder where binary executables are present (system path by default)",
        )
        parser.add_argument(
            "--log-file",
            metavar="[STR]",
            default=None,
            help="Name of the file that we want to use to store the logs (default output to screen)",
        )
        parser.add_argument(
            "--output-folder", metavar="[FOLDER]", action=readable_dir, default=None, help="Path of the output folder"
        )
        parser.add_argument(
            "--temp-folder",
            metavar="[FOLDER]",
            action=readable_dir,
            default=None,
            help="Path of the location for temporary files",
        )
        parser.add_argument(
            "--keep-discarded-files",
            action="store_true",
            default=False,
            help="Keep files with discarded reads in every step",
        )
        parser.add_argument(
            "--qual-64",
            action="store_true",
            default=False,
            help="Use phred-64 quality instead of phred-33(default) in the quality trimming step",
        )
        parser.add_argument(
            "--min-length-qual-trimming",
            default=20,
            metavar="[INT]",
            type=int,
            choices=range(5, 151),
            help="Minimum length of the reads after trimming, "
            "shorter reads will be discarded (default: %(default)s)",
        )
        parser.add_argument(
            "--min-quality-trimming",
            default=20,
            metavar="[INT]",
            type=int,
            choices=range(1, 61),
            help="Minimum phred quality a base must have in order to be kept "
            "in the quality trimming step (default: %(default)s)",
        )
        parser.add_argument(
            "--remove-polyA",
            default=10,
            metavar="[INT]",
            type=int,
            choices=range(0, 35),
            help="Remove PolyA stretches of the given length from R2 " "(Use 0 to disable it) (default: %(default)s)",
        )
        parser.add_argument(
            "--remove-polyT",
            default=10,
            metavar="[INT]",
            type=int,
            choices=range(0, 35),
            help="Remove PolyT stretches of the given length from R2 " "(Use 0 to disable it) (default: %(default)s)",
        )
        parser.add_argument(
            "--remove-polyG",
            default=10,
            metavar="[INT]",
            type=int,
            choices=range(0, 35),
            help="Remove PolyG stretches of the given length from R2 " "(Use 0 to disable it) (default: %(default)s)",
        )
        parser.add_argument(
            "--remove-polyC",
            default=10,
            metavar="[INT]",
            type=int,
            choices=range(0, 35),
            help="Remove PolyC stretches of the given length from R2 " "(Use 0 to disable it) (default: %(default)s)",
        )
        parser.add_argument(
            "--remove-polyN",
            default=10,
            metavar="[INT]",
            type=int,
            choices=range(0, 35),
            help="Remove PolyN stretches of the given length from R2 " "(Use 0 to disable it) (default: %(default)s)",
        )
        parser.add_argument(
            "--homopolymer-mismatches",
            default=0,
            metavar="[INT]",
            type=int,
            choices=range(0, 9),
            help="Number of mismatches allowed when removing " "homopolymers (A, T, G, C or N) (default: %(default)s)",
        )
        parser.add_argument(
            "--filter-AT-content",
            default=90,
            metavar="[INT%]",
            type=int,
            choices=range(0, 101),
            help="Discards reads whose number of A and T bases in total are more "
            "or equal than the percentage given as input (0-100) (default: %(default)s)",
        )
        parser.add_argument(
            "--filter-GC-content",
            default=90,
            metavar="[INT%]",
            type=int,
            choices=range(0, 100),
            help="Discards reads whose number of G and C bases in total are more "
            "or equal the percentage given as input (0-100) (default: %(default)s)",
        )
        parser.add_argument(
            "--mapping-rv-trimming",
            default=0,
            metavar="[INT]",
            type=int,
            choices=range(0, 51),
            help="Number of bases to trim in the reverse reads (R2) for "
            "the mapping step (5' end) (default: %(default)s)",
        )
        parser.add_argument(
            "--inverse-mapping-rv-trimming",
            default=0,
            type=int,
            metavar="[INT]",
            choices=range(0, 51),
            help="Number of bases to trim in the reverse reads (R2) for "
            "the mapping step (3' end) (default: %(default)s)",
        )
        parser.add_argument(
            "--disable-multimap",
            action="store_true",
            default=False,
            help="If activated, multiple aligned reads obtained during mapping will be all discarded. "
            "Otherwise the highest scored one will be kept",
        )
        parser.add_argument(
            "--disable-clipping",
            action="store_true",
            default=False,
            help="If activated, disable soft-clipping (local alignment) in the mapping step",
        )
        parser.add_argument(
            "--min-intron-size",
            default=1,
            metavar="[INT]",
            type=int,
            choices=range(1, 1000),
            help="Minimum allowed intron size when searching for splice variants with STAR\n"
            "Splices alignments are disabled by default (=1) but to turn it on set this"
            "parameter\n to a bigger number, for example 10 or 20. (default: %(default)s)",
        )
        parser.add_argument(
            "--max-intron-size",
            default=1,
            metavar="[INT]",
            type=int,
            choices=range(1, 1000000),
            help="Maximum allowed intron size when searching for splice variants with STAR\n"
            "Splices alignments are disabled by default (=1) but to turn it on set this"
            "parameter\n to a big number, for example 10000 or 100000. (default: %(default)s)",
        )
        parser.add_argument(
            "--star-two-pass-mode",
            default=False,
            action="store_true",
            help="Activates the 2-pass mode in STAR to improve mapping accuracy",
        )
        parser.add_argument(
            "--star-genome-loading",
            default="NoSharedMemory",
            metavar="[STRING]",
            type=str,
            choices=["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"],
            help="Similar to the STAR option --genomeLoad. It allows to load the genome index\n"
            "into memory so it can easily be shared by other jobs to save loading time.\n"
            "Read the STAR manual for more info on this. (default: %(default)s)",
        )
        parser.add_argument(
            "--star-sort-mem-limit",
            default=0,
            type=int,
            help="The maximum available RAM for sorting BAM during mapping with STAR."
            "\nDefault is 0 which means that it will be set to the genome index size",
        )
        parser.add_argument(
            "--demultiplexing-mismatches",
            default=2,
            metavar="[INT]",
            type=int,
            choices=range(0, 31),
            help="Number of allowed mismatches when demultiplexing the reads "
            "against the barcodes with TaggD (default: %(default)s)",
        )
        parser.add_argument(
            "--demultiplexing-kmer",
            default=6,
            metavar="[INT]",
            type=int,
            choices=range(1, 51),
            help="KMer size to use when demultiplexing against the " "barcodes with TaggD (default: %(default)s)",
        )
        parser.add_argument(
            "--demultiplexing-overhang",
            default=0,
            metavar="[INT]",
            type=int,
            choices=range(0, 11),
            help="Extra flanking bases added on each side of the barcode when demultiplexing against "
            "the barcodes with TaggD (default: %(default)s)",
        )
        parser.add_argument(
            "--demultiplexing-start",
            default=0,
            metavar="[INT]",
            type=int,
            help="Start position of the IDs (Barcodes) in R1 (counting from 0) (default: %(default)s)",
        )
        parser.add_argument(
            "--demultiplexing-metric",
            default="Subglobal",
            metavar="[STRING]",
            type=str,
            choices=["Subglobal", "Levenshtein", "Hamming"],
            help="Distance metric to use for TaggD demultiplexing:\n"
            "Options:\n"
            "  Subglobal, Levenshtein or Hamming (default: Subglobal)",
        )
        parser.add_argument(
            "--demultiplexing-multiple-hits-keep-one",
            default=False,
            action="store_true",
            help="When multiple ambiguous hits with same score are "
            "found in the demultiplexing step, keep only one (random).",
        )
        parser.add_argument(
            "--demultiplexing-trim-sequences",
            nargs="+",
            type=int,
            default=None,
            help="Trim the barcodes in the input file when doing demultiplexing.\n"
            "The input given is a list of tuples: START END START END .. where\n"
            "START is the integer position of the first base (0 based) and END is the integer\n"
            "position of the last base (1 based).\n"
            "The final barcode will be obtained by combining all the sequences given in the input.\n"
            "This is useful when having a barcode composed of multiple sequences in the read"
            "or when the barcode needs to be trimmed out.\n"
            "Trimmng sequences can be given several times.",
        )
        parser.add_argument(
            "--demultiplexing-chunk-size",
            default=10000,
            metavar="[INT]",
            type=int,
            help="Chunk size for parallel processing (number of reads assigned to each thread) (default: %(default)s)",
        )
        parser.add_argument(
            "--htseq-mode",
            default="intersection-nonempty",
            type=str,
            metavar="[STRING]",
            choices=["union", "intersection-nonempty", "intersection-strict"],
            help="Mode of annotation when using htseq-count. "
            "Modes = {union, intersection-nonempty(default), intersection-strict}",
        )
        parser.add_argument(
            "--htseq-no-ambiguous",
            action="store_true",
            default=False,
            help="When using htseq-count discard reads annotating ambiguous genes (default False)",
        )
        parser.add_argument(
            "--htseq-features",
            nargs="+",
            default=["exon"],
            type=str,
            help="Which feature types to use from the GTF/GFF file in the annotation.\n "
            "Can be given more than one type (default exon)",
        )
        parser.add_argument(
            "--strandness",
            default="yes",
            type=str,
            metavar="[STRING]",
            choices=["no", "yes", "reverse"],
            help="What strandness mode to use when annotating " "with htseq-count [no, yes(default), reverse]",
        )
        parser.add_argument(
            "--include-non-annotated",
            action="store_true",
            default=False,
            help="Do not discard un-annotated reads (they will be labeled __no_feature)",
        )
        parser.add_argument(
            "--umi-cluster-algorithm",
            default="AdjacentBi",
            metavar="[STRING]",
            type=str,
            choices=["hierarchical", "Adjacent", "AdjacentBi"],
            help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n"
            "Options = {hierarchical, Adjacent and AdjacentBi(default)}.",
        )
        parser.add_argument(
            "--umi-allowed-mismatches",
            default=1,
            metavar="[INT]",
            type=int,
            choices=range(0, 9),
            help="Number of allowed mismatches (hamming distance) "
            "that UMIs of the same gene-spot must have in order to "
            "cluster together (default: %(default)s)",
        )
        parser.add_argument(
            "--umi-start-position",
            default=18,
            metavar="[INT]",
            type=int,
            help="Position in R1 (base wise) of the first base of the " "UMI (starting by 0) (default: %(default)s)",
        )
        parser.add_argument(
            "--umi-end-position",
            default=27,
            metavar="[INT]",
            type=int,
            help="Position in R1 (base wise) of the last base of the " "UMI (starting by 1) (default: %(default)s)",
        )
        parser.add_argument(
            "--umi-filter",
            action="store_true",
            default=False,
            help="Enables the UMI quality filter based on the template given in --umi-filter-template",
        )
        parser.add_argument(
            "--umi-filter-template",
            default="WSNNWSNNV",
            type=str,
            metavar="[STRING]",
            help="UMI template (IUPAC nucleotide code) for the UMI filter, default = WSNNWSNNV",
        )
        parser.add_argument(
            "--umi-quality-bases",
            default=6,
            metavar="[INT]",
            type=int,
            choices=range(0, 13),
            help="Maximum number of low quality bases allowed in an UMI (default: %(default)s)",
        )
        parser.add_argument(
            "--umi-counting-offset",
            default=250,
            metavar="[INT]",
            type=int,
            choices=range(0, 1001),
            help="UMI count for each gene-spot combination is computed "
            "as the number of unique UMIs in each strand/start position. However "
            "some reads might have slightly different start positions due to "
            "amplification artifacts. This parameters allows to define an "
            "offset window from where to count unique UMIs. You can set it to a very "
            "high value +9999 to count unique UMIs for the whole gene (default: %(default)s)",
        )
        parser.add_argument(
            "--compute-saturation",
            action="store_true",
            default=False,
            help="Performs a saturation curve computation by sub-sampling the annotated reads, computing "
            "unique UMIs and adding the stats to the log file (this can be used to plot saturation curves)",
        )
        parser.add_argument(
            "--saturation-points",
            default=None,
            nargs="+",
            type=int,
            help="Saturation points for the saturation curve computation can be "
            "provided instead of using default values.\n"
            "Provide a list of values like this for example: 10000 20000 50000 100000",
        )
        parser.add_argument(
            "--disable-trimming",
            default=False,
            action="store_true",
            help="Use this flag if you want to skip the trimming step",
        )
        parser.add_argument(
            "--disable-mapping",
            default=False,
            action="store_true",
            help="Use this flag if you want to skip the mapping step",
        )
        parser.add_argument(
            "--disable-annotation",
            default=False,
            action="store_true",
            help="Use this flag if you want to skip the annotation",
        )
        parser.add_argument(
            "--disable-barcode",
            default=False,
            action="store_true",
            help="Use this flag if you want to skip the barcode demultiplexing step",
        )
        parser.add_argument(
            "--disable-umi",
            default=False,
            action="store_true",
            help="Use this flag if you want to skip the UMI filtering step",
        )
        parser.add_argument(
            "--transcriptome",
            default=False,
            action="store_true",
            help="Use this flag if you want to use transcriptome instead of a genome, the gene tag will be "
            "obtained from the transcriptome file",
        )
        parser.add_argument("--version", action="version", version="%(prog)s " + str(version_number))
        return parser

    def load_parameters(self, options: argparse.Namespace) -> None:
        """
        Load the input parameters from the argparse object given as parameter
        """
        self.allowed_missed = options.demultiplexing_mismatches
        self.allowed_kmer = options.demultiplexing_kmer
        self.overhang = options.demultiplexing_overhang
        self.min_length_trimming = options.min_length_qual_trimming
        self.trimming_rv = options.mapping_rv_trimming
        self.min_quality_trimming = options.min_quality_trimming
        self.clean = options.no_clean_up
        self.barcode_start = options.demultiplexing_start
        self.threads = options.threads
        self.verbose = options.verbose
        self.ids = os.path.abspath(options.ids)
        if options.ref_map is not None:
            self.ref_map = os.path.abspath(options.ref_map)
        if options.ref_annotation is not None:
            self.ref_annotation = os.path.abspath(options.ref_annotation)
        self.expName = options.expName
        self.htseq_mode = options.htseq_mode
        self.htseq_no_ambiguous = options.htseq_no_ambiguous
        self.htseq_features = options.htseq_features
        self.qual64 = options.qual_64
        if options.contaminant_index is not None:
            self.contaminant_index = os.path.abspath(options.contaminant_index)
        # Load the given path into the system PATH
        if options.bin_path is not None and os.path.isdir(options.bin_path):
            os.environ["PATH"] += os.pathsep + options.bin_path
        if options.log_file is not None:
            self.logfile = os.path.abspath(options.log_file)
        self.fastq_fw = os.path.abspath(options.fastq_files[0])
        self.fastq_rv = os.path.abspath(options.fastq_files[1])
        if options.output_folder is not None and os.path.isdir(options.output_folder):
            self.output_folder = os.path.abspath(options.output_folder)
        else:
            self.output_folder = os.path.abspath(os.getcwd())  # type: ignore
        if options.temp_folder is not None and os.path.isdir(options.temp_folder):
            self.temp_folder = os.path.abspath(options.temp_folder)
        else:
            self.temp_folder = tempfile.mkdtemp(prefix="st_pipeline_temp")  # type: ignore
        self.umi_allowed_mismatches = options.umi_allowed_mismatches
        self.umi_start_position = options.umi_start_position
        self.umi_end_position = options.umi_end_position
        self.keep_discarded_files = options.keep_discarded_files
        self.remove_polyA_distance = options.remove_polyA
        self.remove_polyT_distance = options.remove_polyT
        self.remove_polyG_distance = options.remove_polyG
        self.remove_polyC_distance = options.remove_polyC
        self.remove_polyN_distance = options.remove_polyN
        self.filter_AT_content = options.filter_AT_content
        self.filter_GC_content = options.filter_GC_content
        self.disable_multimap = options.disable_multimap
        self.disable_clipping = options.disable_clipping
        self.umi_cluster_algorithm = options.umi_cluster_algorithm
        self.min_intron_size = options.min_intron_size
        self.max_intron_size = options.max_intron_size
        self.umi_filter = options.umi_filter
        self.umi_filter_template = options.umi_filter_template.upper()
        self.compute_saturation = options.compute_saturation
        self.include_non_annotated = options.include_non_annotated
        self.inverse_trimming_rv = options.inverse_mapping_rv_trimming
        self.two_pass_mode = options.star_two_pass_mode
        self.strandness = options.strandness
        self.umi_quality_bases = options.umi_quality_bases
        self.umi_counting_offset = options.umi_counting_offset
        self.taggd_metric = options.demultiplexing_metric
        self.taggd_multiple_hits_keep_one = options.demultiplexing_multiple_hits_keep_one
        self.taggd_trim_sequences = options.demultiplexing_trim_sequences
        self.taggd_chunk_size = options.demultiplexing_chunk_size
        self.adaptor_missmatches = options.homopolymer_mismatches
        self.star_genome_loading = options.star_genome_loading
        self.star_sort_mem_limit = options.star_sort_mem_limit
        self.disable_barcode = options.disable_barcode
        self.disable_trimming = options.disable_trimming
        self.disable_mapping = options.disable_mapping
        self.disable_annotation = options.disable_annotation
        self.transcriptome = options.transcriptome
        self.disable_umi = options.disable_umi
        self.transcriptome = options.transcriptome
        if options.saturation_points is not None:
            self.saturation_points = [int(p) for p in options.saturation_points]  # type: ignore
        # Assign class parameters to the QA stats object
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        attributes_filtered = {a[0]: a[1] for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))}
        # Assign general parameters to the qa_stats object
        self.qa_stats.input_parameters = attributes_filtered  # type: ignore
        self.qa_stats.annotation_tool = f"htseq-count {get_htseq_count_version()}"
        self.qa_stats.demultiplex_tool = f"Taggd {get_taggd_count_version()}"
        self.qa_stats.pipeline_version = version_number
        self.qa_stats.mapper_tool = f"STAR {get_star_version()}"

    def createLogger(self) -> None:
        """
        Creates a logging object and logs some information from the input parameters.
        """
        # Create a logger
        if self.logfile is not None:
            logging.basicConfig(filename=self.logfile, level=logging.DEBUG)
        else:
            logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

        logger.info(f"ST Pipeline {version_number}")

        # Log general information
        logger.info(f"Output directory: {self.output_folder}")
        logger.info(f"Temporary directory: {self.temp_folder}")
        logger.info(f"Dataset name: {self.expName}")
        logger.info(f"Forward(R1) input file: {self.fastq_fw}")
        logger.info(f"Reverse(R2) input file: {self.fastq_rv}")
        logger.info(f"Reference mapping STAR index folder: {self.ref_map}")
        if self.ref_annotation is not None:
            logger.info(f"Reference annotation file: {self.ref_annotation}")
        if self.contaminant_index is not None:
            logger.info(f"Using contamination filter STAR index: {self.contaminant_index}")
        logger.info(f"CPU Nodes: {self.threads}")

        if not self.disable_trimming:
            logger.info("Quality and trimming settings")
            if self.remove_polyA_distance > 0:
                logger.info(f"Removing polyA sequences of a length of at least: {self.remove_polyA_distance}")
            if self.remove_polyT_distance > 0:
                logger.info(f"Removing polyT sequences of a length of at least: {self.remove_polyT_distance}")
            if self.remove_polyG_distance > 0:
                logger.info(f"Removing polyG sequences of a length of at least: {self.remove_polyG_distance}")
            if self.remove_polyC_distance > 0:
                logger.info(f"Removing polyC sequences of a length of at least: {self.remove_polyC_distance}")
            if self.remove_polyN_distance > 0:
                logger.info(f"Removing polyN sequences of a length of at least: {self.remove_polyN_distance}")
            logger.info(f"Allowing {self.adaptor_missmatches} mismatches when removing homopolymers")
            logger.info(f"Remove reads whose AT content is {self.filter_AT_content}%")
            logger.info(f"Remove reads whose GC content is {self.filter_GC_content}%")

        if not self.disable_mapping:
            logger.info("Mapping settings")
            logger.info(f"Mapping reverse trimming: {self.trimming_rv}")
            logger.info(f"Mapping inverse reverse trimming: {self.inverse_trimming_rv}")
            logger.info("Mapping tool: STAR")
            logger.info(f"Mapping minimum intron size allowed (splice alignments) with STAR: {self.min_intron_size}")
            logger.info(f"Mapping maximum intron size allowed (splice alignments) with STAR: {self.max_intron_size}")
            logger.info(f"STAR genome loading strategy: {self.star_genome_loading}")
            if self.disable_clipping:
                logger.info("Not allowing soft clipping when mapping with STAR")
            if self.disable_multimap:
                logger.info("Not allowing multiple alignments when mapping with STAR")
            if self.two_pass_mode:
                logger.info("Using the STAR 2-pass mode for the mapping step")

        if not self.disable_barcode:
            logger.info("Demultiplexing settings")
            logger.info(f"Ids(barcodes) file: {self.ids}")
            logger.info(f"TaggD allowed mismatches: {self.allowed_missed}")
            logger.info(f"TaggD kmer size: {self.allowed_kmer}")
            logger.info(f"TaggD overhang: {self.overhang}")
            logger.info(f"TaggD metric: {self.taggd_metric}")
            if self.taggd_multiple_hits_keep_one:
                logger.info("TaggD multiple hits keep one (random) is enabled")
            if self.taggd_trim_sequences is not None:
                logger.info(f"TaggD trimming from the barcodes: {'-'.join(str(x) for x in self.taggd_trim_sequences)}")
            logger.info(f"TaggD chunk size: {self.taggd_chunk_size }")

        if not self.disable_annotation:
            logger.info("Annotation settings")
            logger.info("Annotation tool: HTSeq")
            logger.info(f"Annotation mode: {self.htseq_mode}")
            logger.info(f"Annotation strandness: {self.strandness}")
            logger.info(f"Annotation feature types: {','.join(self.htseq_features)}")
            if self.include_non_annotated:
                logger.info("Including non annotated reads in the output")

        if self.compute_saturation:
            logger.info("Computing saturation curve with several sub-samples...")
            if self.saturation_points is not None:
                logger.info(f"Using the following points: {' '.join(str(p) for p in self.saturation_points)}")

        if not self.disable_umi:
            logger.info("UMI settings")
            logger.info(f"UMIs start position: {self.umi_start_position}")
            logger.info(f"UMIs end position: {self.umi_end_position}")
            logger.info(f"UMIs allowed mismatches: {self.umi_allowed_mismatches}")
            logger.info(f"UMIs clustering algorithm: {self.umi_cluster_algorithm}")
            logger.info(
                f"Allowing an offset of {self.umi_counting_offset} when clustering UMIs by strand-start in a gene-spot"
            )
            logger.info(f"Allowing {self.umi_quality_bases} low quality bases in an UMI")
            logger.info(f"Discarding reads that after trimming are shorter than {self.min_length_trimming}")
            if self.umi_filter:
                logger.info(f"UMIs using filter: {self.umi_filter_template}")

    def run(self) -> None:
        """
        Runs the whole pipeline given the parameters present.
        It performs several sequential steps.
        It logs information and running time.
        It throws Exceptions if something went wrong.
        """
        # First adjust the intermediate files with the temp_folder path
        if self.temp_folder:
            for key, value in list(FILENAMES.items()):
                FILENAMES[key] = os.path.join(self.temp_folder, value)
            for key, value in list(FILENAMES_DISCARDED.items()):
                FILENAMES_DISCARDED[key] = os.path.join(self.temp_folder, value)

        # Get the starting time to compute total execution time
        globaltime = TimeStamper()

        # =================================================================
        # START PIPELINE
        # =================================================================
        start_exe_time = globaltime.get_timestamp()
        logger.info(f"Starting the pipeline: {start_exe_time}")

        # =================================================================
        # STEP: FILTERING
        # Applies different filters : sanity, quality, short, adaptors, UMI...
        # =================================================================
        # Get the barcode length
        barcode_length = len(list(read_barcode_file(self.ids).values())[0].sequence)
        if not self.disable_trimming:
            logger.info(f"Start filtering raw reads {globaltime.get_timestamp()}")
            try:
                stats = filter_input_data(
                    self.fastq_fw,
                    self.fastq_rv,
                    FILENAMES["quality_trimmed_R2"],
                    FILENAMES_DISCARDED["quality_trimmed_discarded"] if self.keep_discarded_files else None,
                    barcode_length,
                    self.barcode_start,
                    self.filter_AT_content,
                    self.filter_GC_content,
                    self.umi_start_position,
                    self.umi_end_position,
                    self.min_quality_trimming,
                    self.min_length_trimming,
                    self.remove_polyA_distance,
                    self.remove_polyT_distance,
                    self.remove_polyG_distance,
                    self.remove_polyC_distance,
                    self.remove_polyN_distance,
                    self.qual64,
                    self.umi_filter,
                    self.umi_filter_template,
                    self.umi_quality_bases,
                    self.adaptor_missmatches,
                    self.overhang,
                    self.disable_umi,
                    self.disable_barcode,
                )
                # update qa_stats
                self.qa_stats.input_reads_reverse = stats[0]
                self.qa_stats.reads_after_trimming_forward = stats[1]
                self.qa_stats.reads_after_trimming_reverse = stats[1]
            except Exception:
                raise

        # =================================================================
        # CONDITIONAL STEP: Filter out contaminated reads, e.g. rRNA(Optional)
        # =================================================================
        if self.contaminant_index:
            # To remove contaminants sequence we align the reads to the contaminant genome
            # and keep the un-mapped reads
            logger.info(f"Starting contaminant filter alignment {globaltime.get_timestamp()}")
            try:
                # Make the contaminant filter call
                alignReads(
                    FILENAMES["quality_trimmed_R2"],
                    self.contaminant_index,
                    FILENAMES_DISCARDED["contaminated_discarded"],
                    None,  # Do not pass the annotation file in contaminant filter
                    self.temp_folder,
                    self.trimming_rv,
                    self.inverse_trimming_rv,
                    self.threads,
                    1,  # Disable splice alignments in contaminant filter
                    1,  # Disable splice alignments in contaminant filter
                    False,  # Disable multimap in contaminant filter
                    False,  # Disable softclipping in contaminant filter
                    False,  # Disable 2-pass mode in contaminant filter
                    self.min_length_trimming,
                    True,  # Include un-aligned reads in the output
                    self.star_genome_loading,
                    self.star_sort_mem_limit,
                )
                # Extract the contaminant free reads (not aligned) from the output of STAR
                # NOTE: this will not be needed when STAR allows to chose the discarded
                # reads format (BAM)
                # We also need to set the NH tag to Null so to be able to run STAR again
                infile = pysam.AlignmentFile(FILENAMES_DISCARDED["contaminated_discarded"], "rb")
                out_unmap = pysam.AlignmentFile(FILENAMES["contaminated_clean"], "wb", template=infile)
                temp_name = os.path.join(self.temp_folder, next(tempfile._get_candidate_names()))
                out_map = pysam.AlignmentFile(temp_name, "wb", template=infile)
                for sam_record in infile.fetch(until_eof=True):
                    try:
                        sam_record.set_tag("NH", None)
                        sam_record.set_tag("HI", None)
                        sam_record.set_tag("AS", None)
                        sam_record.set_tag("nM", None)
                        # Do we need to remove more tags?
                        # Should we also reset the CIGAR string?
                    except KeyError:
                        # do nothing
                        pass
                    if sam_record.is_unmapped:
                        out_unmap.write(sam_record)
                    else:
                        out_map.write(sam_record)
                infile.close()
                out_map.close()
                out_unmap.close()
                shutil.move(temp_name, FILENAMES_DISCARDED["contaminated_discarded"])
            except Exception:
                raise

        # =================================================================
        # STEP: Maps against the genome using STAR
        # =================================================================
        if not self.disable_mapping:
            logger.info(f"Starting genome alignment {globaltime.get_timestamp()}")
            input_reads = FILENAMES["contaminated_clean"] if self.contaminant_index else FILENAMES["quality_trimmed_R2"]
            try:
                # Make the alignment call
                alignReads(
                    input_reads,
                    self.ref_map,
                    FILENAMES["mapped"],
                    None,  #  Do not annotate on the fly
                    self.temp_folder,  # type: ignore
                    self.trimming_rv,
                    self.inverse_trimming_rv,
                    self.threads,
                    self.min_intron_size,
                    self.max_intron_size,
                    self.disable_multimap,
                    self.disable_clipping,
                    self.two_pass_mode,
                    self.min_length_trimming,
                    self.keep_discarded_files,
                    self.star_genome_loading,
                    self.star_sort_mem_limit,
                )
                # Remove secondary alignments and un-mapped
                # NOTE: this will not be needed when STAR allows to chose the discarded
                # reads format (BAM)
                if self.keep_discarded_files:
                    temp_name = os.path.join(self.temp_folder, next(tempfile._get_candidate_names()))  # type: ignore
                    # Note use 260 to also discard multiple-alignments
                    command = "samtools view -b -h -F 4 -@ {} -o {} -U {} {}".format(
                        self.threads, temp_name, FILENAMES_DISCARDED["mapped_discarded"], FILENAMES["mapped"]
                    )
                    subprocess.check_call(command, shell=True)
                    os.rename(temp_name, FILENAMES["mapped"])
            except Exception:
                raise

        # =================================================================
        # STEP: DEMULTIPLEX READS Map against the barcodes (Optional)
        # =================================================================
        if not self.disable_barcode:
            logger.info(f"Starting barcode demultiplexing {globaltime.get_timestamp()}")
            try:
                stats = barcodeDemultiplexing(  # type: ignore
                    FILENAMES["mapped"],
                    self.ids,
                    self.allowed_missed,
                    self.allowed_kmer,
                    self.overhang,
                    self.taggd_metric,
                    self.taggd_multiple_hits_keep_one,
                    self.taggd_trim_sequences,
                    self.threads,
                    FILENAMES["demultiplexed_prefix"],  # Prefix for output files
                    self.keep_discarded_files,
                    self.taggd_chunk_size,
                )
                # pdate qa_stats
                self.qa_stats.reads_after_demultiplexing = stats  # type: ignore

                # TODO TaggD does not output the BAM file sorted
                command = "samtools sort -T {}/sort_bam -@ {} -o {} {}".format(
                    self.temp_folder,
                    self.threads,
                    FILENAMES["demultiplexed_matched"],
                    FILENAMES["demultiplexed_matched"],
                )
                subprocess.check_call(command, shell=True)
            except Exception:
                raise

        # =================================================================
        # STEP: annotate using htseq-count or the transcriptome (Optional)
        # =================================================================
        if not self.disable_annotation:
            input_file = (
                FILENAMES["demultiplexed_matched"] if FILENAMES["demultiplexed_matched"] else FILENAMES["mapped"]
            )
            if self.transcriptome:
                logger.info(f"Assigning gene names from transcriptome {globaltime.get_timestamp()}")
                # Iterate the BAM file to set the gene name as the transcriptome"s entry
                flag_read = "rb"
                flag_write = "wb"
                infile = pysam.AlignmentFile(input_file, flag_read)  # type: ignore
                outfile = pysam.AlignmentFile(FILENAMES["annotated"], flag_write, template=infile)  # type: ignore
                for rec in infile.fetch(until_eof=True):
                    # NOTE chrom may have to be trimmed to 250 characters max
                    chrom = infile.get_reference_name(rec.reference_id).split()[0]
                    rec.set_tag("XF", chrom, "Z")
                    outfile.write(rec)
                infile.close()
                outfile.close()
            else:
                logger.info(f"Starting annotation {globaltime.get_timestamp()}")
                try:
                    stats = annotateReads(  # type: ignore
                        input_file,
                        self.ref_annotation,  # type: ignore
                        FILENAMES["annotated"],
                        FILENAMES_DISCARDED["annotated_discarded"] if self.keep_discarded_files else None,
                        self.htseq_mode,
                        self.strandness,
                        self.htseq_no_ambiguous,
                        self.include_non_annotated,
                        self.htseq_features,
                    )
                    # update qa_stats
                    self.qa_stats.reads_after_annotation = stats  # type: ignore
                except Exception:
                    raise

        # =================================================================
        # STEP: compute saturation (Optional)
        # =================================================================
        # To compute saturation points we need the number of annotated reads
        # the fastest way is to get that information from the stats object
        if self.compute_saturation and os.path.isfile(FILENAMES["annotated"]):
            reads = (
                self.qa_stats.reads_after_annotation
                if not self.transcriptome
                else self.qa_stats.reads_after_demultiplexing
            )
            logger.info(f"Starting computing saturation points {globaltime.get_timestamp()}")
            try:
                compute_saturation(
                    reads,
                    FILENAMES["annotated"],
                    self.ref_annotation,  # type: ignore
                    self.umi_cluster_algorithm,
                    self.umi_allowed_mismatches,
                    self.umi_counting_offset,
                    self.disable_umi,
                    self.expName,
                    self.temp_folder,  # type: ignore
                    self.saturation_points,
                )
            except Exception:
                raise

        # =================================================================
        # STEP: Create dataset and remove duplicates
        # =================================================================
        if os.path.isfile(FILENAMES["annotated"]):
            logger.info(f"Starting creating dataset {globaltime.get_timestamp()}")
            try:
                stats = createDataset(  # type: ignore
                    FILENAMES["annotated"],
                    self.output_folder,  # type: ignore
                    self.ref_annotation,
                    self.umi_cluster_algorithm,
                    self.umi_allowed_mismatches,
                    self.umi_counting_offset,
                    self.disable_umi,
                    self.expName,
                    True,  # Verbose
                )
                # update qa_stats
                self.qa_stats.max_genes_feature = stats["max_genes_feature"]  # type: ignore
                self.qa_stats.min_genes_feature = stats["min_genes_feature"]  # type: ignore
                self.qa_stats.max_reads_feature = stats["max_reads_feature"]  # type: ignore
                self.qa_stats.min_reads_feature = stats["min_reads_feature"]  # type: ignore
                self.qa_stats.average_reads_feature = stats["average_reads_feature"]  # type: ignore
                self.qa_stats.average_genes_feature = stats["average_genes_feature"]  # type: ignore
                self.qa_stats.std_reads_feature = stats["std_reads_feature"]  # type: ignore
                self.qa_stats.std_genes_feature = stats["std_genes_feature"]  # type: ignore
                self.qa_stats.reads_after_duplicates_removal = stats["reads_after_duplicates_removal"]  # type: ignore
                self.qa_stats.barcodes_found = stats["barcodes_found"]  # type: ignore
                self.qa_stats.genes_found = stats["genes_found"]  # type: ignore
                self.qa_stats.duplicates_found = stats["duplicates_found"]  # type: ignore
            except Exception:
                raise

        # =================================================================
        # END PIPELINE
        # =================================================================
        # Write stats to JSON
        print(self.qa_stats)
        # self.qa_stats.write_json(os.path.join(self.output_folder, self.expName + "_qa_stats.json"))  # type: ignore

        finish_exe_time = globaltime.get_timestamp()
        total_exe_time = finish_exe_time - start_exe_time
        logger.info(f"Total Execution Time: {total_exe_time}")
