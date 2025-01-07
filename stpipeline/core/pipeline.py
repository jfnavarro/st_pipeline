# coding=utf-8
"""
This is the main object for the ST pipeline.
It contains the Pipeline object that has methods
to parse the input parameters, generate the input parameters,
do sanity check and ultimately run the pipeline.
"""

from stpipeline.common.utils import *
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing
from stpipeline.core.annotation import annotateReads
from stpipeline.common.stats import qa_stats
from stpipeline.common.dataset import createDataset
from stpipeline.common.saturation import computeSaturation
from stpipeline.version import version_number
from taggd.io.barcode_utils import read_barcode_file
from stpipeline.common.filter import filter_input_data
from stpipeline.common.stats import QAStats
import logging
import argparse
import sys
import shutil
import os
import tempfile
import subprocess
import pysam
import inspect
import re

FILENAMES = {"mapped": "mapped.bam",
             "annotated": "annotated.bam",
             "contaminated_clean": "contaminated_clean.bam",
             "demultiplexed_prefix": "demultiplexed",
             "demultiplexed_matched": "demultiplexed_matched.bam",
             "quality_trimmed_R2": "R2_quality_trimmed.bam"}

FILENAMES_DISCARDED = {"mapped_discarded": "mapping_discarded.bam",
                       "contaminated_discarded": "contaminated.bam",
                       "demultiplexed_ambiguos": "demultiplexed_ambiguos.bam",
                       "demultiplexed_unmatched": "demultiplexed_unmatched.bam",
                       "demultiplexed_results": "demultiplexed_results.tsv",
                       "quality_trimmed_discarded": "R2_quality_trimmed_discarded.fastq",
                       "annotated_discarded": "annotated_discarded.bam"}


class Pipeline():
    """
    This class contains all the ST pipeline
    attributes and a bunch of methods to parse
    the input parameters, do sanity check and
    run the pipeline steps.
    """
    LogName = "STPipeline"

    def __init__(self):
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
        self.ids = None
        self.ref_map = None
        self.ref_annotation = None
        self.expName = None
        self.htseq_mode = "intersection-nonempty"
        self.htseq_no_ambiguous = False
        self.htseq_features = ["exon"]
        self.qual64 = False
        self.contaminant_index = None
        self.fastq_fw = None
        self.fastq_rv = None
        self.logger = None
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
        self.qa_stats = QAStats()

    def clean_filenames(self):
        """ Just makes sure to remove
        all temp files
        """
        if self.clean:
            for file_name in list(FILENAMES.values()):
                safeRemove(file_name)
        if not self.keep_discarded_files:
            for file_name in list(FILENAMES_DISCARDED.values()):
                safeRemove(file_name)
        if self.temp_folder is not None and os.path.isdir(self.temp_folder):
            star_temp1 = os.path.join(self.temp_folder, "_STARgenome")
            star_temp2 = os.path.join(self.temp_folder, "_STARpass1")
            if os.path.isdir(star_temp1):
                shutil.rmtree(star_temp1)
            if os.path.isdir(star_temp2):
                shutil.rmtree(star_temp2)
            if self.clean and not self.keep_discarded_files \
                    and self.temp_folder != self.output_folder:
                shutil.rmtree(self.temp_folder)

    def sanityCheck(self):
        """ 
        Performs some basic sanity checks on the input parameters
        """

        if self.ref_annotation is not None and (not os.path.isfile(self.ref_annotation) or
                                                (not self.ref_annotation.endswith(".gtf")
                                                 and not self.ref_annotation.endswith(".gff3")
                                                 and not self.ref_annotation.endswith(".gff"))):
            error = "Error parsing parameters.\n" \
                    "Invalid annotation file {}".format(self.ref_annotation)
            self.logger.error(error)
            raise RuntimeError(error)

        if self.ref_annotation is None and not self.transcriptome:
            error = "Error, annotation file is missing but the transcriptome option is disabled"
            self.logger.error(error)
            raise RuntimeError(error)

        if self.ref_map is None and not self.disable_mapping:
            error = "Error, genome reference is missing but the disable_mapping option is disabled"
            self.logger.error(error)
            raise RuntimeError(error)

        if not os.path.isfile(self.fastq_fw) or not os.path.isfile(self.fastq_rv):
            error = "Error parsing parameters.\n" \
                    "Invalid input files {} {}".format(self.fastq_fw, self.fastq_rv)
            self.logger.error(error)
            raise RuntimeError(error)

        if (not self.fastq_fw.endswith(".fastq")
            and not self.fastq_fw.endswith(".fq")
            and not self.fastq_fw.endswith(".gz")
            and not self.fastq_rv.endswith(".bz2")) \
                or (not self.fastq_rv.endswith(".fastq")
                    and not self.fastq_rv.endswith(".fq")
                    and not self.fastq_rv.endswith(".gz")
                    and not self.fastq_rv.endswith(".bz2")):
            error = "Error parsing parameters.\n" \
                    "Incorrect format for input files {} {}".format(self.fastq_fw, self.fastq_rv)
            self.logger.error(error)
            raise RuntimeError(error)

        if not self.disable_barcode and not os.path.isfile(self.ids):
            error = "Error parsing parameters.\n" \
                    "Invalid IDs file {}".format(self.ids)
            self.logger.error(error)
            raise RuntimeError(error)

        if not self.disable_barcode and self.ids is None:
            error = "Error IDs file is missing but the option to disable the " \
                    "demultiplexing step is not activated\n"
            self.logger.error(error)
            raise RuntimeError(error)

        if self.saturation_points is not None and not self.compute_saturation:
            self.logger.warning("Saturation points are provided but the option"
                                "to compute saturation is disabled.")

        if not self.disable_umi and self.umi_filter:
            # Check template validity
            regex = "[^ACGTURYSWKMBDHVN]"
            if re.search(regex, self.umi_filter_template) is not None:
                error = "Error invalid UMI template given {}.".format(self.umi_filter_template)
                self.logger.error(error)
                raise RuntimeError(error)
            # Check template length
            if len(self.umi_filter_template) != (self.umi_end_position - self.umi_start_position):
                error = "Error the UMI template given does not have the same " \
                        "length as the UMIs {}.".format(self.umi_filter_template)
                self.logger.error(error)
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
            error = "Error starting the pipeline.\n" \
                    "Taggd allowed mismatches is bigger or equal than the Taggd k-mer size"
            self.logger.error(error)
            raise RuntimeError(error)

        if self.umi_start_position < self.barcode_start < self.umi_end_position \
                and not self.disable_barcode and not self.disable_umi:
            error = "Error starting the pipeline.\n" \
                    "The start position of the barcodes is between the UMIs start-end position"
            self.logger.error(error)
            raise RuntimeError(error)

        if self.umi_allowed_mismatches > (self.umi_end_position - self.umi_start_position) \
                and not self.disable_umi:
            error = "Error starting the pipeline.\n" \
                    "The allowed UMI mismatches is bigger than the UMI size"
            self.logger.error(error)
            raise RuntimeError(error)

        # Test the presence of the required tools
        required_scripts = set(["STAR"])
        unavailable_scripts = set()
        for script in required_scripts:
            if which_program(script) is None:
                unavailable_scripts.add(script)
        if len(unavailable_scripts) != 0:
            error = "Error starting the pipeline.\n" \
                    "Required software not found:\t".join(unavailable_scripts)
            self.logger.error(error)
            raise RuntimeError(error)

    def createParameters(self, parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
        """
        Adds the pipeline"s parameters to a given
        Argparse object and returns it.
        """
        class readable_dir(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                prospective_dir = values
                if not os.path.isdir(prospective_dir):
                    raise argparse.ArgumentTypeError("{0} is not a valid path".format(prospective_dir))
                if os.access(prospective_dir, os.R_OK):
                    setattr(namespace, self.dest, prospective_dir)
                else:
                    raise argparse.ArgumentTypeError("{0} is not a readable dir".format(prospective_dir))

        parser.add_argument("fastq_files", nargs=2)
        parser.add_argument("--ids",
                            metavar="[FILE]",
                            required=False,
                            help="Path to the file containing the map of barcodes to the array coordinates")
        parser.add_argument("--ref-map",
                            metavar="[FOLDER]",
                            action=readable_dir,
                            required=False,
                            help="Path to the folder with the STAR index "
                                 "for the genome that you want to use as reference")
        parser.add_argument("--ref-annotation",
                            metavar="[FILE]",
                            required=False,
                            help="Path to the reference annotation file "
                                 "(GTF or GFF format is required) to be used to annotated the mapped reads")
        parser.add_argument("--expName",
                            type=str,
                            metavar="[STRING]",
                            required=True,
                            help="Name of the dataset (The output files will prepend this name)")
        parser.add_argument("--contaminant-index",
                            metavar="[FOLDER]",
                            action=readable_dir,
                            default=None,
                            help="Path to the folder with a STAR index with a contaminant genome reference.\n"
                                 "Reads will be filtered using the specified genome and mapping reads will be discarded")
        parser.add_argument("--no-clean-up",
                            action="store_false",
                            default=True,
                            help="Do not remove temporary/intermediary files (useful for debugging)")
        parser.add_argument("--verbose",
                            action="store_true",
                            default=False,
                            help="Show extra information on the log file")
        parser.add_argument("--threads",
                            default=4,
                            metavar="[INT]",
                            type=int,
                            choices=range(1, 81),
                            help="Number of threads to use (default: %(default)s)")
        parser.add_argument("--bin-path",
                            metavar="[FOLDER]",
                            action=readable_dir,
                            default=None,
                            help="Path to folder where binary executables are present (system path by default)")
        parser.add_argument("--log-file",
                            metavar="[STR]",
                            default=None,
                            help="Name of the file that we want to use to store the logs (default output to screen)")
        parser.add_argument("--output-folder",
                            metavar="[FOLDER]",
                            action=readable_dir,
                            default=None,
                            help="Path of the output folder")
        parser.add_argument("--temp-folder",
                            metavar="[FOLDER]",
                            action=readable_dir,
                            default=None,
                            help="Path of the location for temporary files")
        parser.add_argument("--keep-discarded-files",
                            action="store_true",
                            default=False,
                            help="Keep files with discarded reads in every step")
        parser.add_argument("--qual-64",
                            action="store_true",
                            default=False,
                            help="Use phred-64 quality instead of phred-33(default) in the quality trimming step")
        parser.add_argument("--min-length-qual-trimming",
                            default=20,
                            metavar="[INT]",
                            type=int,
                            choices=range(5, 151),
                            help="Minimum length of the reads after trimming, "
                                 "shorter reads will be discarded (default: %(default)s)")
        parser.add_argument("--min-quality-trimming",
                            default=20,
                            metavar="[INT]",
                            type=int,
                            choices=range(1, 61),
                            help="Minimum phred quality a base must have in order to be kept "
                                 "in the quality trimming step (default: %(default)s)")
        parser.add_argument("--remove-polyA",
                            default=10,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 35),
                            help="Remove PolyA stretches of the given length from R2 "
                                 "(Use 0 to disable it) (default: %(default)s)")
        parser.add_argument("--remove-polyT",
                            default=10,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 35),
                            help="Remove PolyT stretches of the given length from R2 "
                                 "(Use 0 to disable it) (default: %(default)s)")
        parser.add_argument("--remove-polyG",
                            default=10,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 35),
                            help="Remove PolyG stretches of the given length from R2 "
                                 "(Use 0 to disable it) (default: %(default)s)")
        parser.add_argument("--remove-polyC",
                            default=10,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 35),
                            help="Remove PolyC stretches of the given length from R2 "
                                 "(Use 0 to disable it) (default: %(default)s)")
        parser.add_argument("--remove-polyN",
                            default=10,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 35),
                            help="Remove PolyN stretches of the given length from R2 "
                                 "(Use 0 to disable it) (default: %(default)s)")
        parser.add_argument("--homopolymer-mismatches",
                            default=0,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 9),
                            help="Number of mismatches allowed when removing "
                                 "homopolymers (A, T, G, C or N) (default: %(default)s)")
        parser.add_argument("--filter-AT-content",
                            default=90,
                            metavar="[INT%]",
                            type=int,
                            choices=range(0, 101),
                            help="Discards reads whose number of A and T bases in total are more "
                                 "or equal than the percentage given as input (0-100) (default: %(default)s)")
        parser.add_argument("--filter-GC-content",
                            default=90,
                            metavar="[INT%]",
                            type=int,
                            choices=range(0, 100),
                            help="Discards reads whose number of G and C bases in total are more "
                                 "or equal the percentage given as input (0-100) (default: %(default)s)")
        parser.add_argument("--mapping-rv-trimming",
                            default=0,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 51),
                            help="Number of bases to trim in the reverse reads (R2) for "
                                 "the mapping step (5" end) (default: %(default)s)")
        parser.add_argument("--inverse-mapping-rv-trimming",
                            default=0,
                            type=int,
                            metavar="[INT]",
                            choices=range(0, 51),
                            help="Number of bases to trim in the reverse reads (R2) for "
                                 "the mapping step (3" end) (default: %(default)s)")
        parser.add_argument("--disable-multimap",
                            action="store_true",
                            default=False,
                            help="If activated, multiple aligned reads obtained during mapping will be all discarded. "
                                 "Otherwise the highest scored one will be kept")
        parser.add_argument("--disable-clipping",
                            action="store_true",
                            default=False,
                            help="If activated, disable soft-clipping (local alignment) in the mapping step")
        parser.add_argument("--min-intron-size",
                            default=1,
                            metavar="[INT]",
                            type=int,
                            choices=range(1, 1000),
                            help="Minimum allowed intron size when searching for splice variants with STAR\n"
                                 "Splices alignments are disabled by default (=1) but to turn it on set this parameter\n"
                                 "to a bigger number, for example 10 or 20. (default: %(default)s)")
        parser.add_argument("--max-intron-size",
                            default=1,
                            metavar="[INT]",
                            type=int,
                            choices=range(1, 1000000),
                            help="Maximum allowed intron size when searching for splice variants with STAR\n"
                                 "Splices alignments are disabled by default (=1) but to turn it on set this parameter\n"
                                 "to a big number, for example 10000 or 100000. (default: %(default)s)")
        parser.add_argument("--star-two-pass-mode",
                            default=False,
                            action="store_true",
                            help="Activates the 2-pass mode in STAR to improve mapping accuracy")
        parser.add_argument("--star-genome-loading",
                            default="NoSharedMemory",
                            metavar="[STRING]",
                            type=str,
                            choices=["NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit"],
                            help="Similar to the STAR option --genomeLoad. It allows to load the genome index\n"
                                 "into memory so it can easily be shared by other jobs to save loading time.\n"
                                 "Read the STAR manual for more info on this. (default: %(default)s)")
        parser.add_argument("--star-sort-mem-limit",
                            default=0,
                            type=int,
                            help="The maximum available RAM for sorting BAM during mapping with STAR."
                                 "\nDefault is 0 which means that it will be set to the genome index size")
        parser.add_argument("--demultiplexing-mismatches",
                            default=2,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 31),
                            help="Number of allowed mismatches when demultiplexing the reads "
                                 "against the barcodes with TaggD (default: %(default)s)")
        parser.add_argument("--demultiplexing-kmer",
                            default=6,
                            metavar="[INT]",
                            type=int,
                            choices=range(1, 51),
                            help="KMer size to use when demultiplexing against the "
                                 "barcodes with TaggD (default: %(default)s)")
        parser.add_argument("--demultiplexing-overhang",
                            default=0,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 11),
                            help="Extra flanking bases added on each side of the barcode when demultiplexing against "
                                 "the barcodes with TaggD (default: %(default)s)")
        parser.add_argument("--demultiplexing-start",
                            default=0,
                            metavar="[INT]",
                            type=int,
                            help="Start position of the IDs (Barcodes) in R1 (counting from 0) (default: %(default)s)")
        parser.add_argument("--demultiplexing-metric",
                            default="Subglobal",
                            metavar="[STRING]",
                            type=str,
                            choices=["Subglobal", "Levenshtein", "Hamming"],
                            help="Distance metric to use for TaggD demultiplexing:\n"
                                 "Options:\n"
                                 "  Subglobal, Levenshtein or Hamming (default: Subglobal)")
        parser.add_argument("--demultiplexing-multiple-hits-keep-one",
                            default=False,
                            action="store_true",
                            help="When multiple ambiguous hits with same score are "
                                 "found in the demultiplexing step, keep only one (random).")
        parser.add_argument("--demultiplexing-trim-sequences",
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
                                 "Trimmng sequences can be given several times.")
        parser.add_argument("--htseq-mode",
                            default="intersection-nonempty",
                            type=str,
                            metavar="[STRING]",
                            choices=["union", "intersection-nonempty", "intersection-strict"],
                            help="Mode of annotation when using htseq-count. "
                                 "Modes = {union, intersection-nonempty(default), intersection-strict}")
        parser.add_argument("--htseq-no-ambiguous",
                            action="store_true",
                            default=False,
                            help="When using htseq-count discard reads annotating ambiguous genes (default False)")
        parser.add_argument("--htseq-features",
                            nargs="+",
                            default=["exon"],
                            type=str,
                            help="Which feature types to use from the GTF/GFF file in the annotation.\n "
                                 "Can be given more than one type (default exon)")
        parser.add_argument("--strandness",
                            default="yes",
                            type=str,
                            metavar="[STRING]",
                            choices=["no", "yes", "reverse"],
                            help="What strandness mode to use when annotating "
                                 "with htseq-count [no, yes(default), reverse]")
        parser.add_argument("--include-non-annotated",
                            action="store_true",
                            default=False,
                            help="Do not discard un-annotated reads (they will be labeled __no_feature)")
        parser.add_argument("--umi-cluster-algorithm",
                            default="AdjacentBi",
                            metavar="[STRING]",
                            type=str,
                            choices=["hierarchical", "Adjacent", "AdjacentBi"],
                            help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n"
                                 "Options = {hierarchical, Adjacent and AdjacentBi(default)}.")
        parser.add_argument("--umi-allowed-mismatches",
                            default=1,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 9),
                            help="Number of allowed mismatches (hamming distance) "
                                 "that UMIs of the same gene-spot must have in order to "
                                 "cluster together (default: %(default)s)")
        parser.add_argument("--umi-start-position",
                            default=18,
                            metavar="[INT]",
                            type=int,
                            help="Position in R1 (base wise) of the first base of the "
                                 "UMI (starting by 0) (default: %(default)s)")
        parser.add_argument("--umi-end-position",
                            default=27,
                            metavar="[INT]",
                            type=int,
                            help="Position in R1 (base wise) of the last base of the "
                                 "UMI (starting by 1) (default: %(default)s)")
        parser.add_argument("--umi-filter",
                            action="store_true",
                            default=False,
                            help="Enables the UMI quality filter based on the template given in --umi-filter-template")
        parser.add_argument("--umi-filter-template",
                            default="WSNNWSNNV",
                            type=str,
                            metavar="[STRING]",
                            help="UMI template (IUPAC nucleotide code) for the UMI filter, default = WSNNWSNNV")
        parser.add_argument("--umi-quality-bases",
                            default=6,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 13),
                            help="Maximum number of low quality bases allowed in an UMI (default: %(default)s)")
        parser.add_argument("--umi-counting-offset",
                            default=250,
                            metavar="[INT]",
                            type=int,
                            choices=range(0, 1001),
                            help="UMI count for each gene-spot combination is computed "
                                 "as the number of unique UMIs in each strand/start position. However "
                                 "some reads might have slightly different start positions due to "
                                 "amplification artifacts. This parameters allows to define an "
                                 "offset window from where to count unique UMIs. You can set it to a very "
                                 "high value +9999 to count unique UMIs for the whole gene (default: %(default)s)")
        parser.add_argument("--compute-saturation",
                            action="store_true",
                            default=False,
                            help="Performs a saturation curve computation by sub-sampling the annotated reads, computing "
                                 "unique UMIs and adding the stats to the log file (this can be used to plot saturation curves)")
        parser.add_argument("--saturation-points",
                            default=None,
                            nargs="+",
                            type=int,
                            help="Saturation points for the saturation curve computation can be "
                                 "provided instead of using default values.\n"
                                 "Provide a list of values like this for example: 10000 20000 50000 100000")
        parser.add_argument("--disable-trimming",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to skip the trimming step")
        parser.add_argument("--disable-mapping",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to skip the mapping step")
        parser.add_argument("--disable-annotation",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to skip the annotation")
        parser.add_argument("--disable-barcode",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to skip the barcode demultiplexing step")
        parser.add_argument("--disable-umi",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to skip the UMI filtering step")
        parser.add_argument("--transcriptome",
                            default=False,
                            action="store_true",
                            help="Use this flag if you want to use transcriptome instead of a genome, the gene tag will be "
                                 "obtained from the transcriptome file")
        parser.add_argument("--version", action="version", version="%(prog)s " + str(version_number))
        return parser

    def load_parameters(self, options: argparse.ArgumentParser):
        """
        Load the input parameters from the argparse object given as parameter
        :param options: a Argparse object
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
            self.output_folder = os.path.abspath(os.getcwd())
        if options.temp_folder is not None and os.path.isdir(options.temp_folder):
            self.temp_folder = os.path.abspath(options.temp_folder)
        else:
            self.temp_folder = tempfile.mkdtemp(prefix="st_pipeline_temp")
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
            self.saturation_points = [int(p) for p in options.saturation_points]
        # Assign class parameters to the QA stats object
        attributes = inspect.getmembers(self, lambda a: not (inspect.isroutine(a)))
        attributes_filtered = [a for a in attributes if not (a[0].startswith("__") and a[0].endswith("__"))]
        # Assign general parameters to the qa_stats object
        qa_stats.input_parameters = attributes_filtered
        qa_stats.annotation_tool = "htseq-count {}".format(getHTSeqCountVersion())
        qa_stats.demultiplex_tool = "Taggd {}".format(getTaggdCountVersion())
        qa_stats.pipeline_version = version_number
        qa_stats.mapper_tool = getSTARVersion()

    def createLogger(self):
        """
        Creates a logging object and logs some information from the input parameters
        """
        # Create a logger
        if self.logfile is not None:
            logging.basicConfig(filename=self.logfile, level=logging.DEBUG)
        else:
            logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

        self.logger = logging.getLogger(self.__class__.LogName)
        self.logger.info("ST Pipeline {}".format(version_number))

        # Some info
        self.logger.info("Output directory: {}".format(self.output_folder))
        self.logger.info("Temporary directory: {}".format(self.temp_folder))
        self.logger.info("Dataset name: {}".format(self.expName))
        self.logger.info("Forward(R1) input file: {}".format(self.fastq_fw))
        self.logger.info("Reverse(R2) input file: {}".format(self.fastq_rv))
        self.logger.info("Reference mapping STAR index folder: {}".format(self.ref_map))
        if self.ref_annotation is not None:
            self.logger.info("Reference annotation file: {}".format(self.ref_annotation))
        if self.contaminant_index is not None:
            self.logger.info("Using contamination filter STAR index: {}".format(self.contaminant_index))
        self.logger.info("CPU Nodes: {}".format(self.threads))
        if not self.disable_barcode:
            self.logger.info("Ids(barcodes) file: {}".format(self.ids))
            self.logger.info("TaggD allowed mismatches: {}".format(self.allowed_missed))
            self.logger.info("TaggD kmer size: {}".format(self.allowed_kmer))
            self.logger.info("TaggD overhang: {}".format(self.overhang))
            self.logger.info("TaggD metric: {}".format(self.taggd_metric))
            if self.taggd_multiple_hits_keep_one:
                self.logger.info("TaggD multiple hits keep one (random) is enabled")
            if self.taggd_trim_sequences is not None:
                self.logger.info(
                    "TaggD trimming from the barcodes " + "-".join(str(x) for x in self.taggd_trim_sequences))
        if not self.disable_mapping:
            self.logger.info("Mapping reverse trimming: {}".format(self.trimming_rv))
            self.logger.info("Mapping inverse reverse trimming: {}".format(self.inverse_trimming_rv))
            self.logger.info("Mapping tool: STAR")
            self.logger.info(
                "Mapping minimum intron size allowed (splice alignments) with STAR: {}".format(self.min_intron_size))
            self.logger.info(
                "Mapping maximum intron size allowed (splice alignments) with STAR: {}".format(self.max_intron_size))
            self.logger.info("STAR genome loading strategy {}".format(self.star_genome_loading))
            if self.disable_clipping:
                self.logger.info("Not allowing soft clipping when mapping with STAR")
            if self.disable_multimap:
                self.logger.info("Not allowing multiple alignments when mapping with STAR")
            if self.two_pass_mode:
                self.logger.info("Using the STAR 2-pass mode for the mapping step")
        if not self.disable_annotation:
            self.logger.info("Annotation tool: HTSeq")
            self.logger.info("Annotation mode: {}".format(self.htseq_mode))
            self.logger.info("Annotation strandness {}".format(self.strandness))
            self.logger.info("Annotation feature types {}".format(",".join(self.htseq_features)))
            if self.include_non_annotated:
                self.logger.info("Including non annotated reads in the output")
        if self.compute_saturation:
            self.logger.info("Computing saturation curve with several sub-samples...")
            if self.saturation_points is not None:
                self.logger.info(
                    "Using the following points {}".format(" ".join(str(p) for p in self.saturation_points)))
        if not self.disable_umi:
            self.logger.info("UMIs start position: {}".format(self.umi_start_position))
            self.logger.info("UMIs end position: {}".format(self.umi_end_position))
            self.logger.info("UMIs allowed mismatches: {}".format(self.umi_allowed_mismatches))
            self.logger.info("UMIs clustering algorithm: {}".format(self.umi_cluster_algorithm))
            self.logger.info("Allowing an offset of {} when clustering UMIs "
                             "by strand-start in a gene-spot".format(self.umi_counting_offset))
            self.logger.info("Allowing {} low quality bases in an UMI".format(self.umi_quality_bases))
            self.logger.info(
                "Discarding reads that after trimming are shorter than {}".format(self.min_length_trimming))
            if self.umi_filter:
                self.logger.info("UMIs using filter: {}".format(self.umi_filter_template))
        if not self.disable_trimming:
            if self.remove_polyA_distance > 0:
                self.logger.info(
                    "Removing polyA sequences of a length of at least: {}".format(self.remove_polyA_distance))
            if self.remove_polyT_distance > 0:
                self.logger.info(
                    "Removing polyT sequences of a length of at least: {}".format(self.remove_polyT_distance))
            if self.remove_polyG_distance > 0:
                self.logger.info(
                    "Removing polyG sequences of a length of at least: {}".format(self.remove_polyG_distance))
            if self.remove_polyC_distance > 0:
                self.logger.info(
                    "Removing polyC sequences of a length of at least: {}".format(self.remove_polyC_distance))
            if self.remove_polyN_distance > 0:
                self.logger.info(
                    "Removing polyN sequences of a length of at least: {}".format(self.remove_polyN_distance))
            self.logger.info("Allowing {} mismatches when removing homopolymers".format(self.adaptor_missmatches))
            self.logger.info("Remove reads whose AT content is {}%".format(self.filter_AT_content))
            self.logger.info("Remove reads whose GC content is {}%".format(self.filter_GC_content))

    def run(self):
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
        start_exe_time = globaltime.getTimestamp()
        self.logger.info("Starting the pipeline: {}".format(start_exe_time))

        # =================================================================
        # STEP: FILTERING
        # Applies different filters : sanity, quality, short, adaptors, UMI...
        # =================================================================
        # Get the barcode length
        barcode_length = len(list(read_barcode_file(self.ids).values())[0].sequence)
        if not self.disable_trimming:
            self.logger.info("Start filtering raw reads {}".format(globaltime.getTimestamp()))
            try:
                stats = filter_input_data(self.fastq_fw,
                                 self.fastq_rv,
                                 FILENAMES["quality_trimmed_R2"],
                                 FILENAMES_DISCARDED[
                                     "quality_trimmed_discarded"] if self.keep_discarded_files else None,
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
                                 self.disable_barcode)
                # TODO update qa_stats
            except Exception:
                raise

        # =================================================================
        # CONDITIONAL STEP: Filter out contaminated reads, e.g. typically bacterial rRNA
        # =================================================================
        if self.contaminant_index:
            # To remove contaminants sequence we align the reads to the contaminant genome
            # and keep the un-mapped reads
            self.logger.info("Starting contaminant filter alignment {}".format(globaltime.getTimestamp()))
            try:
                # Make the contaminant filter call
                alignReads(FILENAMES["quality_trimmed_R2"],
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
                           self.star_sort_mem_limit)
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
            self.logger.info("Starting genome alignment {}".format(globaltime.getTimestamp()))
            input_reads = FILENAMES["contaminated_clean"] if self.contaminant_index else FILENAMES["quality_trimmed_R2"]
            try:
                # Make the alignment call
                alignReads(input_reads,
                           self.ref_map,
                           FILENAMES["mapped"],
                           None,  #  Do not annotate on the fly
                           self.temp_folder,
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
                           self.star_sort_mem_limit)
                # Remove secondary alignments and un-mapped
                # NOTE: this will not be needed when STAR allows to chose the discarded
                # reads format (BAM)
                if self.keep_discarded_files:
                    temp_name = os.path.join(self.temp_folder, next(tempfile._get_candidate_names()))
                    # Note use 260 to also discard multiple-alignments
                    command = "samtools view -b -h -F 4 -@ {} -o {} -U {} {}".format(self.threads,
                                                                                     temp_name,
                                                                                     FILENAMES_DISCARDED[
                                                                                         "mapped_discarded"],
                                                                                     FILENAMES["mapped"])
                    subprocess.check_call(command, shell=True)
                    os.rename(temp_name, FILENAMES["mapped"])
            except Exception:
                raise

        # =================================================================
        # STEP: DEMULTIPLEX READS Map against the barcodes
        # =================================================================
        if not self.disable_barcode:
            self.logger.info("Starting barcode demultiplexing {}".format(globaltime.getTimestamp()))
            try:
                stats = barcodeDemultiplexing(FILENAMES["mapped"],
                                      self.ids,
                                      self.allowed_missed,
                                      self.allowed_kmer,
                                      self.overhang,
                                      self.taggd_metric,
                                      self.taggd_multiple_hits_keep_one,
                                      self.taggd_trim_sequences,
                                      self.threads,
                                      FILENAMES["demultiplexed_prefix"],  # Prefix for output files
                                      self.keep_discarded_files)
                # TODO update qa_stats.reads_after_demultiplexing
                # TODO TaggD does not output the BAM file sorted
                command = "samtools sort -T {}/sort_bam -@ {} -o {} {}".format(self.temp_folder,
                                                                               self.threads,
                                                                               FILENAMES["demultiplexed_matched"],
                                                                               FILENAMES["demultiplexed_matched"])
                subprocess.check_call(command, shell=True)
            except Exception:
                raise

        # =================================================================
        # STEP: annotate using htseq-count or the transcriptome
        # =================================================================
        if not self.disable_annotation:
            input_file = FILENAMES["demultiplexed_matched"] if FILENAMES["demultiplexed_matched"] else FILENAMES[
                "mapped"]
            if self.transcriptome:
                self.logger.info("Assigning gene names from transcriptome {}".format(globaltime.getTimestamp()))
                # Iterate the BAM file to set the gene name as the transcriptome"s entry
                flag_read = "rb"
                flag_write = "wb"
                infile = pysam.AlignmentFile(input_file, flag_read)
                outfile = pysam.AlignmentFile(FILENAMES["annotated"], flag_write, template=infile)
                for rec in infile.fetch(until_eof=True):
                    # NOTE chrom may have to be trimmed to 250 characters max
                    chrom = infile.getrname(rec.reference_id).split()[0]
                    rec.set_tag("XF", chrom, "Z")
                    outfile.write(rec)
                infile.close()
                outfile.close()
            else:
                self.logger.info("Starting annotation {}".format(globaltime.getTimestamp()))
                try:
                    stats = annotateReads(input_file,
                                  self.ref_annotation,
                                  FILENAMES["annotated"],
                                  FILENAMES_DISCARDED["annotated_discarded"] if self.keep_discarded_files else None,
                                  self.htseq_mode,
                                  self.strandness,
                                  self.htseq_no_ambiguous,
                                  self.include_non_annotated,
                                  self.htseq_features)
                    # TODO update qa_stats.reads_after_annotation
                except Exception:
                    raise

        # =================================================================
        # STEP: compute saturation (Optional)
        # =================================================================
        # To compute saturation points we need the number of annotated reads
        # the fastest way is to get that information from the stats object
        if self.compute_saturation and os.path.isfile(FILENAMES["annotated"]):
            reads = qa_stats.reads_after_annotation if not self.transcriptome else qa_stats.reads_after_demultiplexing
            self.logger.info("Starting computing saturation points {}".format(globaltime.getTimestamp()))
            try:
                computeSaturation(reads,
                                  FILENAMES["annotated"],
                                  self.ref_annotation,
                                  self.umi_cluster_algorithm,
                                  self.umi_allowed_mismatches,
                                  self.umi_counting_offset,
                                  self.disable_umi,
                                  self.expName,
                                  self.temp_folder,
                                  self.saturation_points)
            except Exception:
                raise

        # =================================================================
        # STEP: Create dataset and remove duplicates
        # =================================================================
        if os.path.isfile(FILENAMES["annotated"]):
            self.logger.info("Starting creating dataset {}".format(globaltime.getTimestamp()))
            try:
                stats = createDataset(FILENAMES["annotated"],
                              self.ref_annotation,
                              self.umi_cluster_algorithm,
                              self.umi_allowed_mismatches,
                              self.umi_counting_offset,
                              self.disable_umi,
                              self.output_folder,
                              self.expName,
                              True)  # Verbose
                # TODO update qa_stats
                
            except Exception:
                raise

        # =================================================================
        # END PIPELINE
        # =================================================================
        # Write stats to JSON
        print(qa_stats)
        qa_stats.writeJSON(os.path.join(self.output_folder, self.expName + "_qa_stats.json"))

        finish_exe_time = globaltime.getTimestamp()
        total_exe_time = finish_exe_time - start_exe_time
        self.logger.info("Total Execution Time: {}".format(total_exe_time))
