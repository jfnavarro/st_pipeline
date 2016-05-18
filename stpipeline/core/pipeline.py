#!/usr/bin/env python
""" 
This is the main API for the ST pipeline, it needs a set of files and parameters in order
to run the jobs, input files are fastq, output files are json. It logs status into a file.
"""

from stpipeline.common.utils import *
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing, createIndex
from stpipeline.core.annotation import annotateReads
from stpipeline.common.fastq_utils import reformatRawReads, hashDemultiplexedReads
from stpipeline.common.sam_utils import filterMappedReads
from stpipeline.common.stats import Stats
from stpipeline.common.dataset import createDataset
from stpipeline.common.saturation import computeSaturation
from stpipeline.version import version_number
import logging
import gzip
import argparse

class Pipeline():
    
    LogName = "STPipeline"
    DefaultLogLevel = 'DEBUG'
    
    def __init__(self):
        self.allowed_missed = 2
        self.allowed_kmer = 6
        self.overhang = 2
        self.min_length_trimming = 28
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
        self.path = None
        self.logger = None
        self.logfile = None
        self.output_folder = None
        self.temp_folder = None
        self.molecular_barcodes = False
        self.mc_allowed_mismatches = 1
        self.mc_start_position = 18
        self.mc_end_position = 27
        self.min_cluster_size = 2
        self.keep_discarded_files = False
        self.remove_polyA_distance = 0
        self.remove_polyT_distance = 0
        self.remove_polyG_distance = 0
        self.remove_polyC_distance = 0
        self.filter_AT_content = 90
        self.sam_type = "BAM"
        self.disable_clipping = False
        self.disable_multimap = False
        self.mc_cluster = "naive"
        self.min_intron_size = 20
        self.max_intron_size = 1000000
        self.max_gap_size = 1000000
        # Object to store stats during the run
        self.qa_stats = Stats()
        self.umi_filter = False
        self.umi_filter_template = "WSNNWSNNV"
        self.compute_saturation = False
        self.include_non_annotated = False
        self.inverse_trimming_rv = 0
        self.low_memory = False
        self.two_pass_mode = False
        self.two_pass_mode_genome = None
        self.strandness = "yes"
        
    def sanityCheck(self):
        """ 
        Performs some basic sanity checks in the input parameters
        """

        conds = {"forward_file": fileOk(self.fastq_fw), 
                 "reverse_file": fileOk(self.fastq_rv), 
                 "ref_annotation": fileOk(self.ref_annotation), 
                 "ids_file": fileOk(self.ids),
                 "ref_mapping": self.ref_map is not None and os.path.isdir(self.ref_map), 
                 "Exp Name":  self.expName is not None}
        conds["annotation_extension"] = self.ref_annotation.endswith("gtf") \
                                        or self.ref_annotation.endswith("gff3") \
                                        or self.ref_annotation.endswith("gff")
        conds["annotation_mode"] = self.htseq_mode in ["union","intersection-nonempty","intersection-strict"]
        conds["mc_cluster"] = self.mc_cluster in ["naive", "hierarchical"]
        conds["forward_extension"] = self.fastq_fw.endswith("fastq") \
                                     or self.fastq_fw.endswith("fq") \
                                     or self.fastq_fw.endswith("gz")
        conds["reverse_extension"] = self.fastq_rv.endswith("fastq") \
                                     or self.fastq_rv.endswith("fq") \
                                     or self.fastq_rv.endswith("gz")
        conds["polyA"] = self.remove_polyA_distance >= 0 
        conds["polyT"] = self.remove_polyT_distance >= 0 
        conds["polyG"] = self.remove_polyG_distance >= 0
        conds["polyC"] = self.remove_polyC_distance >= 0
        conds["barcode_start"] = self.barcode_start >= 0
        conds["barcode_length"] = self.barcode_length > 0
        conds["allowed_missed"] = self.allowed_missed >= 0 and self.allowed_missed < self.barcode_length
        conds["allowed_kmer"] = self.allowed_kmer >= 0 and self.allowed_kmer <= self.barcode_length
        conds["overhang"] = self.overhang >= 0 and self.overhang <= self.barcode_length
        conds["min_intron_size"] = self.min_intron_size >= 0 and self.min_intron_size < self.max_intron_size
        conds["max_intron_size"] = self.max_intron_size > 0
        conds["max_gap_size"] = self.max_gap_size > 0
        conds["filter_AT_content"] = self.filter_AT_content >= 0 and self.filter_AT_content <= 100
        conds["sam_type"] = self.sam_type in ["BAM", "SAM"]
        conds["strandness"] = self.strandness in ["yes", "no", "reverse"]
        
        if self.two_pass_mode and self.two_pass_mode_genome and not os.path.isfile(self.two_pass_mode_genome):
            conds["two_pass_mode"] = False
            
        # TODO add checks for trimming parameters
        
        if self.molecular_barcodes:
            conds["molecular_barcodes"] = self.mc_end_position > self.mc_start_position
           
        if self.umi_filter:
            # Check template validity
            import re
            regex = "[^ACGTURYSWKMBDHVN]"
            conds["umi_filter"] = re.search(regex, self.umi_filter_template) is None and self.molecular_barcodes
            
        if not all(conds.values()):
            error = "Error starting the pipeline. Input file or parameters are incorrect %s" % (str(conds))
            self.logger.error(error)
            raise RuntimeError(error)

        # Test the presence of the scripts 
        required_scripts = set(['STAR'])

        unavailable_scripts = set()
        for script in required_scripts:
            if which_program(script) is None: 
                unavailable_scripts.add(script)
         
        if len(unavailable_scripts) == 0:
            self.logger.info("All tools present...Starting the pipeline")
        else:
            error = "Error starting the pipeline. Required software not found:\t".join(unavailable_scripts)
            self.logger.error(error)
            raise RuntimeError(error)
        
        # Add scripts versions to QA Stats
        self.qa_stats.pipeline_version = version_number
        self.qa_stats.mapper_tool = getSTARVersion()
        if self.qa_stats.mapper_tool.find("2.5") == -1:
            error = "Error starting the pipeline. You need STAR 2.5.x or bigger\n"
            self.logger.error(error)
            raise RuntimeError(error)            
        self.qa_stats.annotation_tool = "HTSeqCount " + getHTSeqCountVersion()
        self.qa_stats.demultiplex_tool = "Taggd " + getTaggdCountVersion()
       
    def createParameters(self, parser):
            """
            Adds the pipeline's parameters to a given
            Argparse object 
            """
            class readable_dir(argparse.Action):
                def __call__(self,parser, namespace, values, option_string=None):
                    prospective_dir=values
                    if not os.path.isdir(prospective_dir):
                        raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
                    if os.access(prospective_dir, os.R_OK):
                        setattr(namespace,self.dest,prospective_dir)
                    else:
                        raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

            parser.add_argument('fastq_files', nargs=2)
            parser.add_argument('--ids', metavar="[FILE]", type=argparse.FileType('r'), required=True,
                                help='Path to the file containing the barcodes and the array coordinates.')
            parser.add_argument('--ref-map', metavar="[FOLDER]", action=readable_dir, required=True,
                                help="Path to the folder with the genome STAR index " \
                                "for the genome that you want to use to align the reads")
            parser.add_argument('--ref-annotation', metavar="[FILE]", type=argparse.FileType('r'), required=True,
                                help="Path to the reference annotation file " \
                                "(GTF or GFF format is required)")
            parser.add_argument('--expName', type=str, metavar="[STRING]", required=True,
                                help="Name of the experiment/dataset (The output files will prepend this name)")
            parser.add_argument('--allowed-missed', default=2, metavar="[INT]", type=int, choices=range(0,6),
                                help="Number of allowed mismatches when demultiplexing against the barcodes")
            parser.add_argument('--allowed-kmer', default=6, metavar="[INT]", type=int, choices=range(1,9),
                                help="KMer length when demultiplexing against the barcodes")
            parser.add_argument('--overhang', default=2, metavar="[INT]", type=int, choices=range(0,6),
                                help="Extra flanking bases added when demultiplexing against the barcodes")
            parser.add_argument('--min-length-qual-trimming', default=28, metavar="[INT]", type=int, choices=range(20,100),
                                help="Minimum length of the reads after trimming, " \
                                "shorter reads will be discarded")
            parser.add_argument('--mapping-rv-trimming', default=0, metavar="[INT]", type=int, choices=range(0,100),
                                help="Number of bases to trim in the reverse reads for the mapping step")
            parser.add_argument('--length-id', default=18, type=int, metavar="[INT]", choices=[18,21,24,27],
                                help="Length of IDs (the length of the barcodes)")
            parser.add_argument('--contaminant-index', metavar="[FOLDER]", action=readable_dir,
                                help="Path to the folder with a STAR index with a contaminant genome. Reads will be filtered "
                                "against the specified genome and mapping reads will be descarded")
            parser.add_argument('--qual-64', action="store_true", default=False,
                                help="Use phred-64 quality instead of phred-33(default)")
            parser.add_argument('--htseq-mode', default="intersection-nonempty", type=str, metavar="[STRING]", 
                                choices=["union", "intersection-nonempty", "intersection-strict"],
                                help="Mode of Annotation when using HTSeq. "
                                "Modes = {union,intersection-nonempty(default),intersection-strict}")
            parser.add_argument('--htseq-no-ambiguous', action="store_true", default=False,
                                help="When using htseq discard reads annotating ambiguous genes (default false)")
            parser.add_argument('--start-id', default=0, metavar="[INT]", type=int, choices=range(0,25),
                                help="Start position of the IDs (Barcodes) in the read 1 (counting from 0)")
            parser.add_argument('--no-clean-up', action="store_false", default=True,
                                help="Do not remove temporary/intermediary files (useful for debugging)")
            parser.add_argument('--verbose', action="store_true", default=False,
                                help="Show extra information on the log file")
            parser.add_argument('--mapping-threads', default=4, metavar="[INT]", type=int, choices=range(1,16),
                                help="Number of threads to use in the mapping step")
            parser.add_argument('--min-quality-trimming', default=20, metavar="[INT]", type=int, choices=range(10,60),
                                help="Minimum phred quality for trimming bases in the trimming step")
            parser.add_argument('--bin-path', metavar="[FOLDER]", action=readable_dir,
                                help="Path to folder where binary executables are present (system path by default)")
            parser.add_argument('--log-file', metavar="[FILE]", type=argparse.FileType('r'),
                                help="Name of the file that we want to use to store the logs (default output to screen)")
            parser.add_argument('--output-folder', metavar="[FOLDER]", action=readable_dir,
                                help='Path of the output folder')
            parser.add_argument('--temp-folder', metavar="[FOLDER]", action=readable_dir,
                                help='Path of the location for temporary files')
            parser.add_argument('--molecular-barcodes',
                                action="store_true", 
                                help="Activates the molecular barcodes (UMI) duplicates filter")
            parser.add_argument('--mc-allowed-mismatches', default=1, metavar="[INT]", type=int, choices=range(0,4),
                                help='Number of allowed mismatches when applying the molecular barcodes (UMI) duplicates filter')
            parser.add_argument('--mc-start-position', default=18, metavar="[INT]", type=int, choices=range(0,42),
                                help='Position (base wise) of the first base of the molecular barcodes (starting by 0)')
            parser.add_argument('--mc-end-position', default=27, metavar="[INT]", type=int, choices=range(8,50),
                                help='Position (base wise) of the last base of the molecular barcodes (starting by 1)')
            parser.add_argument('--min-cluster-size', default=2, metavar="[INT]", type=int, choices=range(1,10),
                                help='Min number of equal molecular barcodes to count as a cluster (duplicate)')
            parser.add_argument('--keep-discarded-files', action="store_true", default=False,
                                help='Writes down discarded reads and barcodes into files')
            parser.add_argument('--remove-polyA', default=0, metavar="[INT]", type=int, choices=range(0,25),
                                help="Remove PolyAs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyT', default=0, metavar="[INT]", type=int, choices=range(0,25),
                                help="Remove PolyTs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyG', default=0, metavar="[INT]", type=int, choices=range(0,25),
                                help="Remove PolyGs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyC', default=0, metavar="[INT]", type=int, choices=range(0,25),
                                help="Remove PolyCs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--filter-AT-content', default=90, metavar="[INT%]", type=int, choices=range(1,99),
                                help="Discards reads whose number of A and T bases in total are more " \
                                "or equal than the number given in percentage")
            parser.add_argument('--sam-type', default="BAM", metavar="[STRING]", type=str, choices=["SAM", "BAM"],
                                help="Type of SAM format for intermediate files (SAM or BAM)")
            parser.add_argument('--disable-multimap', action="store_true", default=False,
                                help="If activated, multiple aligned reads obtained during mapping will be all discarded. " \
                                "Otherwise the highest scored one will be kept")
            parser.add_argument('--disable-clipping', action="store_true", default=False,
                                help="If activated, disable soft-clipping (local alignment) in the mapping")
            parser.add_argument('--mc-cluster', default="naive", metavar="[STRING]", type=str, choices=["naive", "hierarchical"],
                                help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                                "Modes = {naive(default), hierarchical}")
            parser.add_argument('--min-intron-size', default=20, metavar="[INT]", type=int, choices=range(0,1000),
                                help="Minimum allowed intron size when searching for splice variants in the mapping step")
            parser.add_argument('--max-intron-size', default=100000, metavar="[INT]", type=int, choices=range(0,1000000),
                                help="Maximum allowed intron size when searching for splice variants in the mapping step")
            parser.add_argument('--max-gap-size', default=1000000, metavar="[INT]", type=int, choices=range(0,1000000),
                                help="Maximum allowed distance between pairs in the mapping step")
            parser.add_argument('--umi-filter', action="store_true", default=False,
                                help="Enables the UMI quality filter based on the template given in --umi-filter-template")
            parser.add_argument('--umi-filter-template', default="WSNNWSNNV", type=str, metavar="[STRING]",
                                help="UMI template for the UMI filter, default = WSNNWSNNV")
            parser.add_argument('--compute-saturation', action="store_true", default=False,
                                help="Performs a saturation curve computation by sub-sampling the annotated reads, computing " \
                                "unique molecules and then a saturation curve (included in the log file)")
            parser.add_argument('--include-non-annotated', action="store_true", default=False,
                                help="Do not discard un-annotated reads (they will be labeled no_feature)")
            parser.add_argument('--inverse-mapping-rv-trimming', default=0, type=int, metavar="[INT]", choices=range(0,100),
                                help="Number of bases to trim in the reverse reads for the mapping step on the 3' end")
            parser.add_argument('--low-memory', default=False, action="store_true",
                                help="Writes temporary records into disk in order to save memory but gaining a speed penalty")
            parser.add_argument('--two-pass-mode', default=False, action="store_true",
                                help="Activates the 2 pass mode in STAR to also map against splice variants")
            parser.add_argument('--two-pass-mode-genome', metavar="[FILE]", type=argparse.FileType('r'),
                                help="Path to a genome file in fasta format. \n" \
                                "When using the two pass mode the path of the fasta file with the genome is needed")
            parser.add_argument('--strandness', default="yes", type=str, metavar="[STRING]", choices=["no", "yes", "reverse"],
                                help="What strandness mode to use when annotating with htseq-count [no, yes(default), reverse]")
            parser.add_argument('--version', action='version',  version='%(prog)s ' + str(version_number))
            return parser
         
    def load_parameters(self, options):
        """
        Load the input parameters from the argparse object
        """
        self.allowed_missed = int(options.allowed_missed)
        self.allowed_kmer = int(options.allowed_kmer)
        self.overhang = int(options.overhang)
        self.min_length_trimming = int(options.min_length_qual_trimming)
        self.trimming_rv = int(options.mapping_rv_trimming)
        self.min_quality_trimming = int(options.min_quality_trimming) 
        self.clean = options.no_clean_up
        self.barcode_start = int(options.start_id)
        self.barcode_length = int(options.length_id)
        self.threads = int(options.mapping_threads)
        self.verbose = options.verbose
        self.ids = os.path.abspath(options.ids)
        self.ref_map = os.path.abspath(options.ref_map)
        self.ref_annotation = os.path.abspath(options.ref_annotation)
        self.expName = options.expName
        self.htseq_mode = options.htseq_mode
        self.htseq_no_ambiguous = options.htseq_no_ambiguous
        self.qual64 = options.qual_64
        self.contaminant_index = options.contaminant_index
        self.path = options.bin_path
        # Load the given path into the system PATH
        if self.path is not None and os.path.isdir(self.path): 
            os.environ["PATH"] += os.pathsep + self.path
        if options.log_file is not None:
            self.logfile = os.path.abspath(options.log_file)  
        self.fastq_fw = options.fastq_files[0]
        self.fastq_rv = options.fastq_files[1]
        if options.output_folder is not None and os.path.isdir(options.output_folder):
            self.output_folder = os.path.abspath(options.output_folder)
        if options.temp_folder is not None and os.path.isdir(options.temp_folder): 
            self.temp_folder = os.path.abspath(options.temp_folder)
        # Set default output and temp folders if erroneous given
        if self.output_folder is None or not os.path.isdir(self.output_folder):
            self.logger.info("Invalid path for output directory -- using current directory instead")
            self.output_folder = os.path.abspath(os.getcwd())
        if self.temp_folder is None or not os.path.isdir(self.temp_folder):
            self.logger.info("Invalid path for temp directory -- using current directory instead")
            self.temp_folder = os.path.abspath(os.getcwd())
        self.molecular_barcodes = options.molecular_barcodes
        self.mc_allowed_mismatches = int(options.mc_allowed_mismatches)
        self.mc_start_position = int(options.mc_start_position)
        self.mc_end_position = int(options.mc_end_position)
        self.min_cluster_size = int(options.min_cluster_size)
        self.keep_discarded_files = options.keep_discarded_files
        self.remove_polyA_distance = int(options.remove_polyA)
        self.remove_polyT_distance = int(options.remove_polyT)
        self.remove_polyG_distance = int(options.remove_polyG)
        self.remove_polyC_distance = int(options.remove_polyC)
        self.filter_AT_content = int(options.filter_AT_content)
        self.sam_type = str(options.sam_type)
        self.disable_multimap = options.disable_multimap
        self.disable_clipping = options.disable_clipping
        self.mc_cluster = str(options.mc_cluster)
        self.min_intron_size = int(options.min_intron_size)
        self.max_intron_size = int(options.max_intron_size)
        self.max_gap_size = int(options.max_gap_size)
        self.umi_filter = options.umi_filter
        self.umi_filter_template = str(options.umi_filter_template).upper()
        self.compute_saturation = options.compute_saturation
        self.include_non_annotated = options.include_non_annotated
        self.inverse_trimming_rv = int(options.inverse_mapping_rv_trimming)
        self.low_memory = options.low_memory
        self.two_pass_mode = options.two_pass_mode
        self.two_pass_mode_genome = options.two_pass_mode_genome
        self.strandness = str(options.strandness)
        # Assign class parameters to the QA stats object
        import inspect
        attributes = inspect.getmembers(self, lambda a:not(inspect.isroutine(a)))
        attributes_filtered = [a for a in attributes if not(a[0].startswith('__') and a[0].endswith('__'))]
        self.qa_stats.input_parameters = attributes_filtered
        
    def createLogger(self):
        """
        Creates a logging object and logs some information from the input parameters
        """    
        # Create a logger
        if self.logfile is not None:
            logging.basicConfig(filename=self.logfile ,level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.DEBUG)
        self.logger = logging.getLogger(self.__class__.LogName)

        self.logger.info("ST Pipeline " + str(version_number))
        
        # Some info
        self.logger.info("Output directory : " + self.output_folder)
        self.logger.info("Temp directory : " + self.temp_folder)
        self.logger.info("Experiment : " + str(self.expName))
        self.logger.info("Forward reads file : " + str(self.fastq_fw))
        self.logger.info("Reverse reads file : " + str(self.fastq_rv))
        self.logger.info("Reference mapping file : " + str(self.ref_map))
        self.logger.info("Reference annotation file : " + str(self.ref_annotation))
        if self.contaminant_index is not None:
            self.logger.info("Using contamination filter with " + str(self.contaminant_index))
        self.logger.info("Nodes : " + str(self.threads))
        self.logger.info("Ids file : " + str(self.ids))
        self.logger.info("TaggD allowed mismatches " + str(self.allowed_missed))
        self.logger.info("TaggD barcode length " + str(self.barcode_length))
        self.logger.info("TaggD kmer size " + str(self.allowed_kmer))
        self.logger.info("TaggD overhang " + str(self.overhang))
        self.logger.info("Mapping reverse trimming " + str(self.trimming_rv))
        self.logger.info("Mapping inverse reverse trimming " + str(self.inverse_trimming_rv))
        self.logger.info("Mapper : STAR")
        self.logger.info("Annotation Tool : HTSeq")
        self.logger.info("Annotation mode " + str(self.htseq_mode))
        self.logger.info("Annotation strandness " + str(self.strandness))
        self.logger.info("Filter of AT content in reads : " + str(self.filter_AT_content))
        self.logger.info("Sam type : " + str(self.sam_type))
        if self.disable_clipping:
            self.logger.info("Not allowing soft clipping when mapping")
        if self.disable_multimap:
            self.logger.info("Not allowing multiple alignments")
        self.logger.info("Mapping minimum intron size " + str(self.min_intron_size))
        self.logger.info("Mapping maximum intron size " + str(self.max_intron_size))
        self.logger.info("Mapping maximum gap size " + str(self.max_gap_size))
        if self.compute_saturation:
            self.logger.info("Computing saturation curve")
        if self.include_non_annotated:
            self.logger.info("Including non annotated reads")
        if self.molecular_barcodes:
            self.logger.info("Using Molecular Barcodes")
            self.logger.info("Molecular Barcode start position " + str(self.mc_start_position))
            self.logger.info("Molecular Barcode end position " + str(self.mc_end_position))
            self.logger.info("Molecular Barcode min cluster size " + str(self.min_cluster_size))
            self.logger.info("Molecular Barcode allowed mismatches " + str(self.mc_allowed_mismatches))
            self.logger.info("Molecular Barcode clustering algorithm " + str(self.mc_cluster))
            if self.umi_filter:
                self.logger.info("Molecular Barcodes using filter " + str(self.umi_filter_template))    
        if self.remove_polyA_distance > 0:
            self.logger.info("Removing polyA adaptors of a length at least " + str(self.remove_polyA_distance))                        
        if self.remove_polyT_distance > 0:
            self.logger.info("Removing polyT adaptors of a length at least " + str(self.remove_polyT_distance))
        if self.remove_polyG_distance > 0:
            self.logger.info("Removing polyG adaptors of a length at least " + str(self.remove_polyG_distance))
        if self.remove_polyC_distance > 0:
            self.logger.info("Removing polyC adaptors of a length at least " + str(self.remove_polyC_distance))
        if self.low_memory:
            self.logger.info("Using the low memory option")
        if self.two_pass_mode :
            self.logger.info("Using the STAR 2 pass mode with genome " + str(self.two_pass_mode_genome))
        
    def run(self):
        """ 
        Runs the whole pipeline given the parameters present
        """
        globaltime = TimeStamper()

        #=================================================================
        # START PIPELINE
        #=================================================================
        start_exe_time = globaltime.getTimestamp()
        self.logger.info("Starting the pipeline : " + str(start_exe_time))

        # Check if input fastq files are gzipped
        # TODO it is faster to make a system call to gunzip
        # TODO add support to bzip 
        if self.fastq_fw.endswith("gz"):
            temp_fastq_fw = os.path.join(self.temp_folder, "unzipped_fastq_fw.fastq")
            with gzip.open(self.fastq_fw, "rb") as filehandler_read:
                with open(temp_fastq_fw, "w") as filehandler_write:
                    for line in filehandler_read:
                        filehandler_write.write(line)
            self.fastq_fw = temp_fastq_fw
        if self.fastq_rv.endswith("gz"):
            temp_fastq_rv = os.path.join(self.temp_folder, "unzipped_fastq_rv.fastq")
            with gzip.open(self.fastq_rv, "rb") as filehandler_read:
                with open(temp_fastq_rv, "w") as filehandler_write:
                    for line in filehandler_read:
                        filehandler_write.write(line)
            self.fastq_rv = temp_fastq_rv
              
        #=================================================================
        # STEP: FILTERING 
        # Applies different filters : sanity, quality, short, adaptors, UMI...
        #=================================================================
        # NOTE after the trimming :
        #    - discarded reads will be replaced by Ns
        #    - The trimming is only performed in the reverse reads
        self.logger.info("Start Quality Filtering raw reads " + str(globaltime.getTimestamp()))
        fastq_rv_trimmed = reformatRawReads(self.fastq_fw,
                                            self.fastq_rv,
                                            self.qa_stats,
                                            self.barcode_start,
                                            self.barcode_length,
                                            self.filter_AT_content,
                                            self.molecular_barcodes,
                                            self.mc_start_position,
                                            self.mc_end_position,
                                            self.trimming_rv,
                                            self.min_quality_trimming,
                                            self.min_length_trimming,
                                            self.remove_polyA_distance,
                                            self.remove_polyT_distance,
                                            self.remove_polyG_distance,
                                            self.remove_polyC_distance,
                                            self.qual64,
                                            self.temp_folder,
                                            self.keep_discarded_files,
                                            self.umi_filter,
                                            self.umi_filter_template)
          
        #=================================================================
        # CONDITIONAL STEP: Filter out contaminated reads, e.g. typically bacterial rRNA
        #=================================================================
        if self.contaminant_index:
            self.logger.info("rRNA filter")
            # To remove contaminants sequence we align the reads to the contaminant genome
            # and keep the un-mapped reads
            old_rv_trimmed = fastq_rv_trimmed
            self.logger.info("Starting rRNA filter alignment " + str(globaltime.getTimestamp()))
            contaminated_sam, fastq_rv_trimmed = alignReads(fastq_rv_trimmed,
                                                            self.contaminant_index,
                                                            self.trimming_rv,
                                                            self.threads,
                                                            "rRNA_filter",
                                                            self.min_intron_size,
                                                            self.max_intron_size,
                                                            self.max_gap_size,
                                                            False, # disable splice variants alignments
                                                            self.sam_type,
                                                            False, # enable multimap in rRNA filter
                                                            True, # disable softclipping in rRNA filter
                                                            0, # do not use the inverse filter for now
                                                            self.temp_folder)
            if self.clean: safeRemove(old_rv_trimmed)
            if not self.keep_discarded_files: safeRemove(contaminated_sam)
            
        #=================================================================
        # STEP: Maps against the genome using STAR
        #=================================================================
        self.logger.info("Genome mapping")
        self.logger.info("Starting genome alignment " + str(globaltime.getTimestamp()))
        sam_mapped, unmapped_reverse = alignReads(fastq_rv_trimmed,
                                                  self.ref_map,
                                                  self.trimming_rv,
                                                  self.threads,
                                                  "alignment",
                                                  self.min_intron_size,
                                                  self.max_intron_size,
                                                  self.max_gap_size,
                                                  True, # enable splice variants alignments
                                                  self.sam_type,
                                                  self.disable_multimap,
                                                  self.disable_clipping,
                                                  self.inverse_trimming_rv,
                                                  self.temp_folder)
        
        if self.two_pass_mode:
            self.logger.info("2 Pass mode creating index " + str(globaltime.getTimestamp()))
            tmp_index = createIndex(self.two_pass_mode_genome,
                                    self.threads,
                                    self.temp_folder)
            safeRemove(sam_mapped)
            safeRemove(unmapped_reverse)
            self.logger.info("2 Pass mode re-alignment " + str(globaltime.getTimestamp()))
            sam_mapped, unmapped_reverse = alignReads(fastq_rv_trimmed,
                                                      tmp_index,
                                                      self.trimming_rv,
                                                      self.threads,
                                                      "alignment",
                                                      self.min_intron_size,
                                                      self.max_intron_size,
                                                      self.max_gap_size,
                                                      True, # enable splice variants alignments
                                                      self.sam_type,
                                                      self.disable_multimap,
                                                      self.disable_clipping,
                                                      self.inverse_trimming_rv,
                                                      self.temp_folder)
            import shutil
            if os.path.exists(tmp_index): shutil.rmtree(tmp_index)
            
        if self.clean: safeRemove(fastq_rv_trimmed)
        if not self.keep_discarded_files: safeRemove(unmapped_reverse)
            
        #=================================================================
        # STEP: DEMULTIPLEX READS Map against the barcodes
        #=================================================================
        self.logger.info("Starting barcode demultiplexing " + str(globaltime.getTimestamp()))
        fastq_fw_demultiplexed = barcodeDemultiplexing(self.fastq_fw,
                                                       self.ids,
                                                       self.qa_stats,
                                                       self.allowed_missed,
                                                       self.allowed_kmer,
                                                       self.barcode_start,
                                                       self.overhang,
                                                       self.threads,
                                                       self.temp_folder,
                                                       self.keep_discarded_files)
        
        #=================================================================
        # STEP: OBTAIN HASH OF DEMULTIPLEXED READS
        # Hash demultiplexed reads to obtain a hash of read_name => (barcode,x,y,umi) 
        #=================================================================
        self.logger.info("Parsing demultiplexed reads " + str(globaltime.getTimestamp()))
        hash_reads = hashDemultiplexedReads(fastq_fw_demultiplexed, 
                                            self.molecular_barcodes, 
                                            self.mc_start_position,
                                            self.mc_end_position,
                                            self.low_memory)
        if self.clean: safeRemove(fastq_fw_demultiplexed)
        
        #================================================================
        # STEP: filters mapped reads and add the (Barcode,x,y,umi) as SAM tags
        #================================================================
        self.logger.info("Starting genome alignment processing unmapped reads " + str(globaltime.getTimestamp()))
        sam_mapped_clean = filterMappedReads(sam_mapped,
                                             self.qa_stats,
                                             hash_reads,
                                             self.min_length_trimming,
                                             self.temp_folder, 
                                             self.keep_discarded_files)
        
        if self.clean: safeRemove(sam_mapped)
        if self.low_memory: 
            hash_reads.close() 
        else: 
            del hash_reads
        
        #=================================================================
        # STEP: SORT sam file with mapped reads by read name or position
        #=================================================================
        # HTSeq-count needs a sorted file but the file comes from the aligner sorted already so
        # NO need to sort for now
        # sam_sorted = sortSamFile(sam_mapped_clean, self.temp_folder)
        # if self.clean: safeRemove(sam_mapped_clean)
        
        #=================================================================
        # STEP: annotate using htseq count
        #=================================================================
        self.logger.info("Starting annotation " + str(globaltime.getTimestamp()))
        annotatedFilteredFile = annotateReads(sam_mapped_clean,
                                              self.ref_annotation,
                                              self.qa_stats,
                                              self.htseq_mode,
                                              self.strandness,
                                              self.htseq_no_ambiguous,
                                              self.include_non_annotated,
                                              self.temp_folder)
        if self.clean: safeRemove(sam_mapped_clean) 

        #=================================================================
        # STEP: filter out un-annotated reads
        # We use a pathed version of htseq-count which does the filtering already
        #=================================================================

        if self.compute_saturation:
            self.logger.info("Starting computing saturation points " + str(globaltime.getTimestamp()))
            computeSaturation(self.qa_stats.reads_after_annotation,
                              annotatedFilteredFile,
                              self.molecular_barcodes,
                              self.mc_cluster,
                              self.mc_allowed_mismatches,
                              self.min_cluster_size,
                              self.expName,
                              self.low_memory,
                              self.temp_folder)
                
        #=================================================================
        # STEP: Create dataset and remove duplicates
        #=================================================================
        self.logger.info("Starting creating dataset " + str(globaltime.getTimestamp()))
        createDataset(annotatedFilteredFile,
                      self.qa_stats,
                      self.molecular_barcodes,
                      self.mc_cluster,
                      self.mc_allowed_mismatches,
                      self.min_cluster_size,
                      self.output_folder,
                      self.expName,
                      True, # Verbose
                      self.low_memory)
        if self.clean: safeRemove(annotatedFilteredFile)

        #=================================================================
        # END PIPELINE
        #=================================================================
        # Write stats to JSON
        self.qa_stats.writeJSON(os.path.join(self.output_folder, self.expName + "_qa_stats.json"))
        
        finish_exe_time = globaltime.getTimestamp()
        total_exe_time = finish_exe_time - start_exe_time
        self.logger.info("Total Execution Time : " + str(total_exe_time))
