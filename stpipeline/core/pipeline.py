""" 
This is the main object for the ST pipeline.
It contains the Pipeline object that has methods
to parse the input parameters, generate the input parameters,
do sanity check and ultimately run the pipeline.
"""

from stpipeline.common.utils import *
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing
from stpipeline.core.annotation import annotateReads
from stpipeline.common.fastq_utils import filterInputReads, hashDemultiplexedReads
from stpipeline.common.sam_utils import filterMappedReads
from stpipeline.common.stats import qa_stats
from stpipeline.common.dataset import createDataset
from stpipeline.common.saturation import computeSaturation
from stpipeline.version import version_number
import logging
import argparse
import sys
import shutil
import os
import gzip
import bz2
import tempfile
import shutil

FILENAMES = {"mapped" : "mapped.bam",
             "annotated" : "annotated.bam",
             "contaminated_clean" : "contaminated_clean.fastq",
             "demultiplexed_prefix" : "demultiplexed",
             "demultiplexed_matched" : "demultiplexed_matched.fastq",
             "mapped_filtered" : "mapped_filtered.bam",
             "quality_trimmed_R1" : "R1_quality_trimmed.fastq",
             "quality_trimmed_R2" : "R2_quality_trimmed.fastq",
             "two_pass_splices" : "SJ.out.tab"}

FILENAMES_DISCARDED = {"mapped_discarded" : "mapping_discarded.fastq",
                       "contaminated_discarded" : "contaminated.bam",
                       "demultiplexed_ambiguos" : "demultiplexed_ambiguos.fastq",
                       "demultiplexed_unmatched" : "demultiplexed_unmatched.fastq",
                       "demultiplexed_results" : "demultiplexed_results.tsv",
                       "mapped_filtered_discarded" : "mapped_filtered_discarded.bam",
                       "quality_trimmed_discarded" : "R2_quality_trimmed_discarded.fastq"}

class Pipeline():
    """ This class contains all the ST pipeline
    attributes and a bunch of methods to parse
    the input parameters, do sanity check and
    run the pipeline steps.
    """
    LogName = "STPipeline"
    
    def __init__(self):
        self.allowed_missed = 2
        self.allowed_kmer = 6
        self.overhang = 2
        self.min_length_trimming = 25
        self.trimming_rv = 0
        self.min_quality_trimming = 20 
        self.clean = True
        self.barcode_start = 0
        self.barcode_length = 18
        self.threads = 8
        self.verbose = False
        self.ids = None
        self.ref_map = None
        self.ref_annotation = None
        self.expName = None
        self.htseq_mode = "intersection-nonempty"
        self.htseq_no_ambiguous = False
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
        self.umi_min_cluster_size = 2
        self.keep_discarded_files = False
        self.remove_polyA_distance = 0
        self.remove_polyT_distance = 0
        self.remove_polyG_distance = 0
        self.remove_polyC_distance = 0
        self.filter_AT_content = 90
        self.filter_GC_content = 90
        self.disable_clipping = False
        self.disable_multimap = False
        self.umi_cluster_algorithm = "naive"
        self.min_intron_size = 20
        self.max_intron_size = 1000000
        self.umi_filter = False
        self.umi_filter_template = "WSNNWSNNV"
        self.compute_saturation = False
        self.include_non_annotated = False
        self.inverse_trimming_rv = 0
        self.low_memory = False
        self.two_pass_mode = False
        self.strandness = "yes"
        self.umi_quality_bases = 4
        self.umi_counting_offset = 50
        
    def clean_filenames(self):
        """ Just makes sure to remove
        all temp files
        """
        if self.clean:
            for file_name in FILENAMES.itervalues():
                safeRemove(file_name)
        if not self.keep_discarded_files:
            for file_name in FILENAMES_DISCARDED.itervalues():
                safeRemove(file_name)
        if self.temp_folder is not None and os.path.isdir(self.temp_folder):
            safeRemove(os.path.join(self.temp_folder,"unzipped_fastq_fw.fastq"))
            safeRemove(os.path.join(self.temp_folder,"unzipped_fastq_rv.fastq"))
            if os.path.isdir(os.path.join(self.temp_folder, "_STARgenome")):
                shutil.rmtree(os.path.join(self.temp_folder, "_STARgenome"))
            if os.path.isdir(os.path.join(self.temp_folder, "_STARpass1")):
                shutil.rmtree(os.path.join(self.temp_folder, "_STARpass1"))
            if self.clean and not self.keep_discarded_files and os.path.isdir(self.temp_folder):
                shutil.rmtree(self.temp_folder)
          
    def sanityCheck(self):
        """ 
        Performs some basic sanity checks on the input parameters
        """

        if not os.path.isfile(self.ref_annotation) or \
        (not self.ref_annotation.endswith(".gtf") \
        and not self.ref_annotation.endswith(".gff3") \
        and not self.ref_annotation.endswith(".gff")):
            error = "Error parsing parameters.\n" \
            "Invalid annotation file {}".format(self.ref_annotation)
            self.logger.error(error)
            raise RuntimeError(error)
          
        if not os.path.isfile(self.fastq_fw) or not os.path.isfile(self.fastq_rv):
            error = "Error parsing parameters.\n" \
            "Invalid input files {} {}".format(self.fastq_fw, self.fastq_rv)
            self.logger.error(error)
            raise RuntimeError(error)
                     
        if (not self.fastq_fw.endswith(".fastq") \
        and not self.fastq_fw.endswith(".fq") \
        and not self.fastq_fw.endswith(".gz") \
        and not self.fastq_rv.endswith(".bz2")) \
        or (not self.fastq_rv.endswith(".fastq") \
        and not self.fastq_rv.endswith(".fq") \
        and not self.fastq_rv.endswith(".gz") \
        and not self.fastq_rv.endswith(".bz2")):
            error = "Error parsing parameters.\n" \
            "Incorrect format for input files {} {}".format(self.fastq_fw, self.fastq_rv)
            self.logger.error(error)
            raise RuntimeError(error)
             
        if not os.path.isfile(self.ids):
            error = "Error parsing parameters.\n" \
            "Invalid IDs file {}".format(self.ids)
            self.logger.error(error)
            raise RuntimeError(error)      
           
        if self.umi_filter:
            # Check template validity
            import re
            regex = "[^ACGTURYSWKMBDHVN]"
            if re.search(regex, self.umi_filter_template) is not None:
                error = "Error invalid UMI template given {}.\n".format(self.umi_filter_template)
                self.logger.error(error)
                raise RuntimeError(error)                
            # Check template length
            if len(self.umi_filter_template) != (self.umi_end_position - self.umi_start_position):
                error = "Error the UMI template given does not have the same " \
                "length as the UMIs {}.\n".format(self.umi_filter_template)
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
                         
        # TODO add checks for trimming parameters, demultiplex parameters and UMI parameters
                
        # Test the presence of the scripts 
        required_scripts = set(['STAR'])
        unavailable_scripts = set()
        for script in required_scripts:
            if which_program(script) is None: 
                unavailable_scripts.add(script)
        if len(unavailable_scripts) != 0:
            error = "Error starting the pipeline.\n" \
            "Required software not found:\t".join(unavailable_scripts)
            self.logger.error(error)
            raise RuntimeError(error) 
       
    def createParameters(self, parser):
        """
        Adds the pipeline's parameters to a given
        Argparse object
        :param parser: the Argparse object
        :return: a new Argparse object with the parameters
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

        parser.add_argument('fastq_files', nargs=2)
        parser.add_argument('--ids', metavar="[FILE]", required=True,
                            help='Path to the file containing the map of barcodes to the array coordinates')
        parser.add_argument('--ref-map', metavar="[FOLDER]", action=readable_dir, required=True,
                            help="Path to the folder with the STAR index " \
                            "for the genome that you want to use to align the reads")
        parser.add_argument('--ref-annotation', metavar="[FILE]", required=True,
                            help="Path to the reference annotation file " \
                            "(GTF or GFF format is required) to be used to annotated the reads")
        parser.add_argument('--expName', type=str, metavar="[STRING]", required=True,
                            help="Name of the experiment/dataset (The output files will prepend this name)")
        parser.add_argument('--allowed-missed', default=2, metavar="[INT]", type=int, choices=range(0, 7),
                            help="Number of allowed mismatches when demultiplexing " \
                            "against the barcodes with TaggD (default: %(default)s)")
        parser.add_argument('--allowed-kmer', default=6, metavar="[INT]", type=int, choices=range(1, 10),
                            help="KMer length when demultiplexing against the barcodes with TaggD (default: %(default)s)")
        parser.add_argument('--overhang', default=2, metavar="[INT]", type=int, choices=range(0, 7),
                            help="Extra flanking bases added when demultiplexing against the barcodes")
        parser.add_argument('--min-length-qual-trimming', default=25, metavar="[INT]", type=int, choices=range(10, 101),
                            help="Minimum length of the reads after trimming, " \
                            "shorter reads will be discarded (default: %(default)s)")
        parser.add_argument('--mapping-rv-trimming', default=0, metavar="[INT]", type=int, choices=range(0, 101),
                            help="Number of bases to trim in the reverse reads for the mapping step (5' end) (default: %(default)s)")
        parser.add_argument('--length-id', default=18, type=int, metavar="[INT]", choices=[18, 21, 24, 27],
                            help="Length of IDs (the length of the barcodes) (default: %(default)s)")
        parser.add_argument('--contaminant-index', metavar="[FOLDER]", action=readable_dir, default=None,
                            help="Path to the folder with a STAR index with a contaminant genome. Reads will be filtered "
                            "against the specified genome and mapping reads will be descarded")
        parser.add_argument('--qual-64', action="store_true", default=False,
                            help="Use phred-64 quality instead of phred-33(default)")
        parser.add_argument('--htseq-mode', default="intersection-nonempty", type=str, metavar="[STRING]",
                            choices=["union", "intersection-nonempty", "intersection-strict"],
                            help="Mode of Annotation when using HTSeq. "
                            "Modes = {union,intersection-nonempty(default),intersection-strict}")
        parser.add_argument('--htseq-no-ambiguous', action="store_true", default=False,
                            help="When using htseq discard reads annotating ambiguous genes (default False)")
        parser.add_argument('--start-id', default=0, metavar="[INT]", type=int, choices=range(0, 27),
                            help="Start position of the IDs (Barcodes) in the R1 (counting from 0) (default: %(default)s)")
        parser.add_argument('--no-clean-up', action="store_false", default=True,
                            help="Do not remove temporary/intermediary files (useful for debugging)")
        parser.add_argument('--verbose', action="store_true", default=False,
                            help="Show extra information on the log file")
        parser.add_argument('--mapping-threads', default=4, metavar="[INT]", type=int, choices=range(1, 17),
                            help="Number of threads to use in the mapping step (default: %(default)s)")
        parser.add_argument('--min-quality-trimming', default=20, metavar="[INT]", type=int, choices=range(10, 61),
                            help="Minimum phred quality a base must have in the trimming step (default: %(default)s)")
        parser.add_argument('--bin-path', metavar="[FOLDER]", action=readable_dir, default=None,
                            help="Path to folder where binary executables are present (system path by default)")
        parser.add_argument('--log-file', metavar="[STR]", default=None,
                            help="Name of the file that we want to use to store the logs (default output to screen)")
        parser.add_argument('--output-folder', metavar="[FOLDER]", action=readable_dir, default=None,
                            help='Path of the output folder')
        parser.add_argument('--temp-folder', metavar="[FOLDER]", action=readable_dir, default=None,
                            help='Path of the location for temporary files')
        parser.add_argument('--umi-allowed-mismatches', default=1, metavar="[INT]", type=int, choices=range(0, 5),
                            help="Number of allowed mismatches (hamming distance) " \
                            "that UMIs of the same gene-spot must have in order to cluster together (default: %(default)s)")
        parser.add_argument('--umi-start-position', default=18, metavar="[INT]", type=int, choices=range(0, 43),
                            help="Position in R1 (base wise) of the first base of the " \
                            "UMI (starting by 0) (default: %(default)s)")
        parser.add_argument('--umi-end-position', default=27, metavar="[INT]", type=int, choices=range(8, 51),
                            help="Position in R1 (base wise) of the last base of the "\
                            "UMI (starting by 1) (default: %(default)s)")
        parser.add_argument('--umi-min-cluster-size', default=2, metavar="[INT]", type=int, choices=range(1, 11),
                            help="Min number of equal UMIs to count" \
                            " as a cluster (duplicate) given the allowed mismatches (default: %(default)s)")
        parser.add_argument('--keep-discarded-files', action="store_true", default=False,
                            help='Writes down unaligned, un-annotated and un-demultiplexed reads to files')
        parser.add_argument('--remove-polyA', default=0, metavar="[INT]", type=int, choices=range(0, 50),
                            help="Remove PolyA stretches of the given length from R2 (default: %(default)s)")
        parser.add_argument('--remove-polyT', default=0, metavar="[INT]", type=int, choices=range(0, 50),
                            help="Remove PolyT stretches of the given length from R2 (default: %(default)s)")
        parser.add_argument('--remove-polyG', default=0, metavar="[INT]", type=int, choices=range(0, 50),
                            help="Remove PolyG stretches of the given length from R2 (default: %(default)s)")
        parser.add_argument('--remove-polyC', default=0, metavar="[INT]", type=int, choices=range(0, 50),
                            help="Remove PolyC stretches of the given length from R2 (default: %(default)s)")
        parser.add_argument('--filter-AT-content', default=90, metavar="[INT%]", type=int, choices=range(0, 100),
                            help="Discards reads whose number of A and T bases in total are more " \
                            "or equal than the number given in percentage (default: %(default)s)")
        parser.add_argument('--filter-GC-content', default=90, metavar="[INT%]", type=int, choices=range(0, 100),
                            help="Discards reads whose number of G and C bases in total are more " \
                            "or equal than the number given in percentage (default: %(default)s)")
        parser.add_argument('--disable-multimap', action="store_true", default=False,
                            help="If activated, multiple aligned reads obtained during mapping will be all discarded. " \
                            "Otherwise the highest scored one will be kept")
        parser.add_argument('--disable-clipping', action="store_true", default=False,
                            help="If activated, disable soft-clipping (local alignment) in the mapping step")
        parser.add_argument('--umi-cluster-algorithm', default="naive", metavar="[STRING]", 
                            type=str, choices=["naive", "hierarchical"],
                            help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                            "Modes = {naive(default), hierarchical}")
        parser.add_argument('--min-intron-size', default=20, metavar="[INT]", type=int, choices=range(0, 1000),
                            help="Minimum allowed intron size when searching for splice " \
                            "variants in the mapping step (default: %(default)s)")
        parser.add_argument('--max-intron-size', default=100000, metavar="[INT]", type=int, choices=range(0, 1000000),
                            help="Maximum allowed intron size when searching for splice " \
                            "variants in the mapping step (default: %(default)s)")
        parser.add_argument('--umi-filter', action="store_true", default=False,
                            help="Enables the UMI quality filter based on the template given in --umi-filter-template")
        parser.add_argument('--umi-filter-template', default="WSNNWSNNV", type=str, metavar="[STRING]",
                            help="UMI template (IUPAC nucleotide code) for the UMI filter, default = WSNNWSNNV")
        parser.add_argument('--compute-saturation', action="store_true", default=False,
                            help="Performs a saturation curve computation by sub-sampling the annotated reads, computing " \
                            "unique molecules and then a saturation curve (included in the log file)")
        parser.add_argument('--include-non-annotated', action="store_true", default=False,
                            help="Do not discard un-annotated reads (they will be labeled __no_feature)")
        parser.add_argument('--inverse-mapping-rv-trimming', default=0, type=int, metavar="[INT]", choices=range(0, 100),
                            help="Number of bases to trim in the reverse reads for the mapping step on the 3' end")
        parser.add_argument('--low-memory', default=False, action="store_true",
                            help="Writes temporary records into disk in order to save memory but gaining a speed penalty")
        parser.add_argument('--two-pass-mode', default=False, action="store_true",
                            help="Activates the 2 pass mode in STAR to also map against splice variants")
        parser.add_argument('--strandness', default="yes", type=str, metavar="[STRING]", choices=["no", "yes", "reverse"],
                            help="What strandness mode to use when annotating with htseq-count [no, yes(default), reverse]")
        parser.add_argument('--umi-quality-bases', default=4, metavar="[INT]", type=int, choices=range(0, 10),
                            help="Maximum number of low quality bases allowed in an UMI (default: %(default)s)")
        parser.add_argument('--umi-counting-offset', default=50, metavar="[INT]", type=int, choices=range(0, 300),
                            help="Expression count for each gene-spot combination is expressed " \
                            "as the number of unique UMIs in each strand/start position. However " \
                            "some reads might have slightly different start positions due to " \
                            "amplification artifacts. This parameters allows to define an " \
                            "offset from where to count unique UMIs (default: %(default)s)") 
        parser.add_argument('--version', action='version', version='%(prog)s ' + str(version_number))
        return parser
         
    def load_parameters(self, options):
        """
        Load the input parameters from the argparse object given as parameter
        :param options: a Argparse object
        """
        self.allowed_missed = options.allowed_missed
        self.allowed_kmer = options.allowed_kmer
        self.overhang = options.overhang
        self.min_length_trimming = options.min_length_qual_trimming
        self.trimming_rv = options.mapping_rv_trimming
        self.min_quality_trimming = options.min_quality_trimming
        self.clean = options.no_clean_up
        self.barcode_start = options.start_id
        self.barcode_length = options.length_id
        self.threads = options.mapping_threads
        self.verbose = options.verbose
        self.ids = os.path.abspath(options.ids)
        self.ref_map = os.path.abspath(options.ref_map)
        self.ref_annotation = os.path.abspath(options.ref_annotation)
        self.expName = options.expName
        self.htseq_mode = options.htseq_mode
        self.htseq_no_ambiguous = options.htseq_no_ambiguous
        self.qual64 = options.qual_64
        self.contaminant_index = options.contaminant_index
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
        self.umi_min_cluster_size = options.umi_min_cluster_size
        self.keep_discarded_files = options.keep_discarded_files
        self.remove_polyA_distance = options.remove_polyA
        self.remove_polyT_distance = options.remove_polyT
        self.remove_polyG_distance = options.remove_polyG
        self.remove_polyC_distance = options.remove_polyC
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
        self.low_memory = options.low_memory
        self.two_pass_mode = options.two_pass_mode
        self.strandness = options.strandness
        self.umi_quality_bases = options.umi_quality_bases
        self.umi_counting_offset = options.umi_counting_offset
        
        # Assign class parameters to the QA stats object
        import inspect
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        attributes_filtered = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
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
        self.logger.info("Temp directory: {}".format(self.temp_folder))
        self.logger.info("Experiment name: {}".format(self.expName))
        self.logger.info("Forward input file: {}".format(self.fastq_fw))
        self.logger.info("Reverse input file: {}".format(self.fastq_rv))
        self.logger.info("Reference mapping index folder: {}".format(self.ref_map))
        self.logger.info("Reference annotation file: {}".format(self.ref_annotation))
        if self.contaminant_index is not None:
            self.logger.info("Using contamination filter: {}".format(self.contaminant_index))
        self.logger.info("CPU Nodes: {}".format(self.threads))
        self.logger.info("Ids file: {}".format(self.ids))
        self.logger.info("TaggD allowed mismatches: {}".format(self.allowed_missed))
        self.logger.info("TaggD barcode length: {}".format(self.barcode_length))
        self.logger.info("TaggD kmer size: {}".format(self.allowed_kmer))
        self.logger.info("TaggD overhang: {}".format(self.overhang))
        self.logger.info("Mapping reverse trimming: {}".format(self.trimming_rv))
        self.logger.info("Mapping inverse reverse trimming: {}".format(self.inverse_trimming_rv))
        self.logger.info("Mapping tool: STAR")
        self.logger.info("Annotation tool: HTSeq")
        self.logger.info("Annotation mode: {}".format(self.htseq_mode))
        self.logger.info("Annotation strandness {}".format(self.strandness))
        self.logger.info("Filter of AT content in reads: {}".format(self.filter_AT_content))
        self.logger.info("Filter of GC content in reads: {}".format(self.filter_GC_content))
        if self.disable_clipping:
            self.logger.info("Not allowing soft clipping when mapping")
        if self.disable_multimap:
            self.logger.info("Not allowing multiple alignments when mapping")
        self.logger.info("Mapping minimum intron size: {}".format(self.min_intron_size))
        self.logger.info("Mapping maximum intron size: {}".format(self.max_intron_size))
        if self.compute_saturation:
            self.logger.info("Computing saturation curve")
        if self.include_non_annotated:
            self.logger.info("Including non annotated reads")
        self.logger.info("Using UMIs to remove duplicates")
        self.logger.info("Using UMIs start position: {}".format(self.umi_start_position))
        self.logger.info("Using UMIs end position: {}".format(self.umi_end_position))
        self.logger.info("Using UMIs min cluster size: {}".format(self.umi_min_cluster_size))
        self.logger.info("Using UMIs allowed mismatches: {}".format(self.umi_allowed_mismatches))
        self.logger.info("Using UMIs clustering algorithm: {}".format(self.umi_cluster_algorithm))
        self.logger.info("Allowing an offset of {} when clustering UMIs " \
                         "by strand-start in a gene-spot".format(self.umi_counting_offset))
        self.logger.info("Allowing {} low quality bases in an UMI".format(self.umi_quality_bases))
        if self.umi_filter:
            self.logger.info("UMIs using filter: {}".format(self.umi_filter_template))    
        if self.remove_polyA_distance > 0:
            self.logger.info("Removing polyA adaptors of a length of at least: {}".format(self.remove_polyA_distance))                        
        if self.remove_polyT_distance > 0:
            self.logger.info("Removing polyT adaptors of a length of at least: {}".format(self.remove_polyT_distance))
        if self.remove_polyG_distance > 0:
            self.logger.info("Removing polyG adaptors of a length of at least: {}".format(self.remove_polyG_distance))
        if self.remove_polyC_distance > 0:
            self.logger.info("Removing polyC adaptors of a length of at least: {}".format(self.remove_polyC_distance))
        if self.low_memory:
            self.logger.info("Using a SQL based container to save memory")
        if self.two_pass_mode :
            self.logger.info("Using the STAR 2-pass mode for the mapping step")
        
    def run(self):
        """ 
        Runs the whole pipeline given the parameters present.
        It performs several sequencial steps.
        It logs information and running time.
        It throws Exceptions if something went wrong.
        """
        # First adjust the intermediate files with the temp_folder path
        if self.temp_folder:
            for key, value in FILENAMES.iteritems():
                FILENAMES[key] = os.path.join(self.temp_folder, value)
            for key, value in FILENAMES_DISCARDED.iteritems():
                FILENAMES_DISCARDED[key] = os.path.join(self.temp_folder, value)
                      
        # Get the starting time to compute total execution time
        globaltime = TimeStamper()

        #=================================================================
        # START PIPELINE
        #=================================================================
        start_exe_time = globaltime.getTimestamp()
        self.logger.info("Starting the pipeline: {}".format(start_exe_time))

        # Check if input fastq files are compressed
        # TODO it is faster to make a system call with gunzip/bzip
        # TODO reliable way to test if files are compressed (something more robust than just file name endings)
        try:
            if self.fastq_fw.endswith(".gz"):
                temp_fastq_fw = os.path.join(self.temp_folder, "unzipped_fastq_fw.fastq")
                with gzip.open(self.fastq_fw, "rb") as filehandler_read:
                    with open(temp_fastq_fw, "w") as filehandler_write:
                        for line in filehandler_read:
                            filehandler_write.write(line)
                self.fastq_fw = temp_fastq_fw
            elif self.fastq_fw.endswith(".bz2"):
                temp_fastq_fw = os.path.join(self.temp_folder, "unzipped_fastq_fw.fastq")
                with bz2.BZ2File(self.fastq_fw, "rb") as filehandler_read:
                    with open(temp_fastq_fw, "w") as filehandler_write:
                        for line in filehandler_read:
                            filehandler_write.write(line)
                self.fastq_fw = temp_fastq_fw
                
            if self.fastq_rv.endswith(".gz"):
                temp_fastq_rv = os.path.join(self.temp_folder, "unzipped_fastq_rv.fastq")
                with gzip.open(self.fastq_rv, "rb") as filehandler_read:
                    with open(temp_fastq_rv, "w") as filehandler_write:
                        for line in filehandler_read:
                            filehandler_write.write(line)
                self.fastq_rv = temp_fastq_rv
            elif self.fastq_rv.endswith(".bz2"):
                    temp_fastq_rv = os.path.join(self.temp_folder, "unzipped_fastq_rv.fastq")
                    with bz2.BZ2File(self.fastq_rv, "rb") as filehandler_read:
                        with open(temp_fastq_rv, "w") as filehandler_write:
                            for line in filehandler_read:
                                filehandler_write.write(line)
                    self.fastq_rv = temp_fastq_rv
                    
        except Exception as e:
            self.logger.error("Error decompressing GZIP/BZIP2 input files {0} {1}".format(self.fastq_fw, self.fastq_rv))
            raise e

        #=================================================================
        # STEP: FILTERING 
        # Applies different filters : sanity, quality, short, adaptors, UMI...
        #=================================================================
        self.logger.info("Start filtering raw reads {}".format(globaltime.getTimestamp()))
        try: 
            filterInputReads(self.fastq_fw,
                             self.fastq_rv,
                             FILENAMES["quality_trimmed_R1"],
                             FILENAMES["quality_trimmed_R2"],
                             FILENAMES_DISCARDED["quality_trimmed_discarded"] if self.keep_discarded_files else None,
                             self.barcode_start,
                             self.barcode_length,
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
                             self.qual64,
                             self.umi_filter,
                             self.umi_filter_template,
                             self.umi_quality_bases)
        except Exception:
            raise
          
        #=================================================================
        # CONDITIONAL STEP: Filter out contaminated reads, e.g. typically bacterial rRNA
        #=================================================================
        if self.contaminant_index:
            # To remove contaminants sequence we align the reads to the contaminant genome
            # and keep the un-mapped reads
            self.logger.info("Starting contaminant filter alignment {}".format(globaltime.getTimestamp()))
            try:
                alignReads(FILENAMES["quality_trimmed_R2"],
                           self.contaminant_index,
                           FILENAMES_DISCARDED["contaminated_discarded"],
                           None, # Do not pass the annotation file in contaminant filter
                           FILENAMES["contaminated_clean"],
                           self.temp_folder,
                           self.trimming_rv,
                           self.inverse_trimming_rv,
                           self.threads,
                           self.min_intron_size,
                           self.max_intron_size,
                           False, # Enable multimap in contaminant filter
                           True, # Disable softclipping in contaminant filter
                           False) # Disable 2-pass mode in contaminant filter
            except Exception:
                raise
            
        #=================================================================
        # STEP: Maps against the genome using STAR
        #=================================================================
        self.logger.info("Starting genome alignment {}".format(globaltime.getTimestamp()))
        input_reads = FILENAMES["contaminated_clean"] if self.contaminant_index else FILENAMES["quality_trimmed_R2"]
        try:
            alignReads(input_reads,
                       self.ref_map,
                       FILENAMES["mapped"],
                       self.ref_annotation,
                       FILENAMES_DISCARDED["mapped_discarded"],
                       self.temp_folder,
                       self.trimming_rv,
                       self.inverse_trimming_rv,
                       self.threads,
                       self.min_intron_size,
                       self.max_intron_size,
                       self.disable_multimap,
                       self.disable_clipping,
                       self.two_pass_mode)
        except Exception:
            raise
            
        #=================================================================
        # STEP: DEMULTIPLEX READS Map against the barcodes
        #=================================================================
        self.logger.info("Starting barcode demultiplexing {}".format(globaltime.getTimestamp()))
        try:
            barcodeDemultiplexing(FILENAMES["quality_trimmed_R1"],
                                  self.ids,
                                  self.allowed_missed,
                                  self.allowed_kmer,
                                  self.barcode_start,
                                  self.overhang,
                                  self.threads,
                                  FILENAMES["demultiplexed_prefix"], # Prefix for output files
                                  self.keep_discarded_files)
        except Exception:
            raise
        
        #=================================================================
        # STEP: OBTAIN HASH OF DEMULTIPLEXED READS
        # Hash demultiplexed reads to obtain a hash of read_name => (barcode,x,y,umi) 
        #=================================================================
        self.logger.info("Parsing demultiplexed reads {}".format(globaltime.getTimestamp()))
        hash_reads = hashDemultiplexedReads(FILENAMES["demultiplexed_matched"], 
                                            self.umi_start_position,
                                            self.umi_end_position,
                                            self.low_memory)
        
        #================================================================
        # STEP: filters mapped reads and add the (Barcode,x,y,umi) as SAM tags
        #================================================================
        self.logger.info("Starting processing aligned reads {}".format(globaltime.getTimestamp()))
        try:
            filterMappedReads(FILENAMES["mapped"],
                              hash_reads,
                              FILENAMES["mapped_filtered"],
                              FILENAMES_DISCARDED["mapped_filtered_discarded"] if self.keep_discarded_files else None,
                              self.min_length_trimming)
        except Exception:
            raise
        finally:
            if self.low_memory: hash_reads.close() 
                
        #=================================================================
        # STEP: annotate using htseq-count
        #=================================================================
        self.logger.info("Starting annotation {}".format(globaltime.getTimestamp()))
        try:
            annotateReads(FILENAMES["mapped_filtered"],
                          self.ref_annotation,
                          FILENAMES["annotated"],
                          self.htseq_mode,
                          self.strandness,
                          self.htseq_no_ambiguous,
                          self.include_non_annotated)
        except Exception:
            raise

        #=================================================================
        # STEP: compute saturation (Optional)
        #=================================================================
        # To compute saturation points we need the number of annotated reads
        # the fastest way is to get that information from the stats object
        if self.compute_saturation:
            annotated_reads = qa_stats.reads_after_annotation
            self.logger.info("Starting computing saturation points {}".format(globaltime.getTimestamp()))
            try:
                computeSaturation(annotated_reads,
                                  FILENAMES["annotated"],
                                  self.umi_cluster_algorithm,
                                  self.umi_allowed_mismatches,
                                  self.umi_min_cluster_size,
                                  self.umi_counting_offset,
                                  self.expName,
                                  self.temp_folder)
            except Exception:
                raise
                
        #=================================================================
        # STEP: Create dataset and remove duplicates
        #=================================================================
        self.logger.info("Starting creating dataset {}".format(globaltime.getTimestamp()))
        try:
            createDataset(FILENAMES["annotated"],
                          qa_stats, # Passed as reference
                          self.umi_cluster_algorithm,
                          self.umi_allowed_mismatches,
                          self.umi_min_cluster_size,
                          self.umi_counting_offset,
                          self.output_folder,
                          self.expName,
                          True) # Verbose
        except Exception:
            raise

        #=================================================================
        # END PIPELINE
        #=================================================================
        # Write stats to JSON
        qa_stats.writeJSON(os.path.join(self.output_folder, self.expName + "_qa_stats.json"))
        
        finish_exe_time = globaltime.getTimestamp()
        total_exe_time = finish_exe_time - start_exe_time
        self.logger.info("Total Execution Time: {}".format(total_exe_time))
