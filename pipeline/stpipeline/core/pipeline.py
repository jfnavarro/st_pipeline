#!/usr/bin/env python
""" 
This is the main API for the ST pipeline, it needs a bunch of files and parameters in order
to run the jobs, input files are fastq, output files are json. It logs status into a file.
"""

from stpipeline.common.utils import *
from stpipeline.core.mapping import alignReads, barcodeDemultiplexing
from stpipeline.core.annotation import annotateReads
from stpipeline.common.fastq_utils import reformatRawReads, filter_rRNA_reads
from stpipeline.common.sam_utils import filterAnnotatedReads, filterMappedReads
from stpipeline.common.stats import Stats
from stpipeline.version import version_number
import logging
import subprocess
import gzip
from collections import defaultdict

class Pipeline():
    
    LogName = "STPipeline"
    DefaultLogLevel = 'DEBUG'
    
    def __init__(self):
        self.allowed_missed = 3
        self.allowed_kmer = 6
        self.overhang = 2
        self.min_length_trimming = 28
        self.trimming_fw = 42
        self.trimming_rv = 5
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
        self.pair_mode_keep = "reverse"
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
        #TODO for now htseq-count does not support BAM files
        #http://sourceforge.net/p/htseq/bugs/12/
        self.sam_type = "SAM"
        self.disable_clipping = False
        self.disable_multimap = False
        self.mc_cluster = "naive"
        self.min_intron_size = 20
        self.max_intron_size = 1000000
        self.max_gap_size = 1000000
        # To store stats from different parts
        self.qa_stats = Stats()
        self.umi_filter = False
        self.umi_filter_template = "WSNNWSNNV"
        self.compute_saturation = False
        
    def sanityCheck(self):
        """ 
        Performs some basic sanity checks in the input parameters
        """

        conds = {"forward_file": fileOk(self.fastq_fw), 
                 "reverse_file": fileOk(self.fastq_rv), 
                 "ref_annotation": fileOk(self.ref_annotation), 
                 "ref_mapping": self.ref_map is not None, 
                 "Exp Name":  self.expName is not None}
        
        conds["annotation_extension"] = self.ref_annotation.endswith("gtf") \
                                        or self.ref_annotation.endswith("gff3") \
                                        or self.ref_annotation.endswith("gff")
        conds["annotation_mode"] = self.htseq_mode in ["union","intersection-nonempty","intersection-strict"]
        conds["mc_cluster"] = self.mc_cluster in ["naive", "counttrie", "hierarchical"]
        conds["forward_extension"] = self.fastq_fw.endswith("fastq") \
                                     or self.fastq_fw.endswith("fq") \
                                     or self.fastq_fw.endswith("gz")
        conds["forward_extension"] = self.fastq_rv.endswith("fastq") \
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
        conds["pair_keep_mode"] = self.pair_mode_keep in ["forward", "reverse", "both"]
        conds["filter_AT_content"] = self.filter_AT_content >= 0 and self.filter_AT_content <= 100
        conds["sam_type"] = self.sam_type in ["BAM", "SAM"]
        
        # TODO add checks for trimming parameters
        
        if self.molecular_barcodes:
            conds["molecular_barcodes"] = self.mc_end_position > self.mc_start_position
           
        if self.umi_filter:
            #Check template validity
            import re
            regex = "[^ACGTURYSWKMBDHVN]"
            conds["umi_filter"] = re.search(regex, self.umi_filter_template) is None
            
        if not all(conds.values()):
            error = "Error: required file/s and or parameters not " \
            " found or incorrect parameters :" + str(conds)
            self.logger.error(error)
            raise RuntimeError(error)

        # Test the presence of the scripts 
        required_scripts = set(['STAR', 'taggd_demultiplex.py'])

        unavailable_scripts = set()
        for script in required_scripts:
            if which(script) is None: 
                unavailable_scripts.add(script)
         
        if len(unavailable_scripts) == 0:
            self.logger.info("All tools present..starting the analysis")
        else:
            error = "Error, these programs not found:\t".join(unavailable_scripts)
            self.logger.error(error)
            raise RuntimeError(error)
        
        # Add scripts versions to QA Stats
        self.qa_stats.pipeline_version = version_number
        self.qa_stats.mapper_tool = getSTARVersion()
        self.qa_stats.annotation_tool = "HTSeqCount " + getHTSeqCountVersion()
        self.qa_stats.demultiplex_tool = "Taggd " + getTaggdCountVersion()
       
    def createParameters(self, parser):
            """
            Adds the pipeline's parameters to a given
            Argparse object 
            """
            parser.add_argument('fastq_files', nargs=2)
            parser.add_argument('--ids',
                                help='The name of the file containing the barcodes and the coordinates. Keep empty for normal RNA analysis')
            parser.add_argument('--ref-map',
                                help="<path_to_genome_indexes> = Reference genome index " \
                                "for the genome that you want to use to align the reads")
            parser.add_argument('--ref-annotation',
                                help="Path to the reference annotation file " \
                                "(htseq requires a GTF file annotation file that you want to use to annotate")
            parser.add_argument('--expName', help="Name of the experiment (output file name)")
            parser.add_argument('--allowed-missed', default=6, 
                                help="Number of allowed mismatches when mapping against the barcodes")
            parser.add_argument('--allowed-kmer', default=7, 
                                help="KMer length when mapping against the barcodes")
            parser.add_argument('--overhang', default=2,
                                help="Extra flanking bases added when mapping against the barcodes")
            parser.add_argument('--min-length-qual-trimming', default=28,
                                help="Minimum length of the sequence for mapping after trimming, " \
                                "shorter reads will be discarded")
            parser.add_argument('--mapping-fw-trimming', default=42,
                                help="Number of bases to trim in the forward reads for the Mapping [24 + ID_LENGTH]")
            parser.add_argument('--mapping-rv-trimming', default=5,
                                help="Number of bases to trim in the reverse reads for the Mapping")
            parser.add_argument('--length-id', default=18, help="Length of IDs (the length of the barcodes)")
            parser.add_argument('--contaminant-index',
                                help="<path_to_genome_indexes> = When provided, reads will be filtered "
                                "against the specified genome index, non-mapping reads will be saved and demultiplexed")
            parser.add_argument('--qual-64', action="store_true", default=False,
                                help="Use phred-64 quality instead of phred-33(default)")
            parser.add_argument('--htseq-mode', default="intersection-nonempty",
                                help="Mode of Annotation when using HTSeq. "
                                "Modes = {union,intersection-nonempty(default),intersection-strict}")
            parser.add_argument('--htseq-no-ambiguous', action="store_true",
                                help="When using htseq discard reads annotating ambiguous genes")
            parser.add_argument('--start-id', default=0, help="Start position of the IDs (Barcodes) in the reads from 0")
            parser.add_argument('--no-clean-up', action="store_false", default=True,
                                help="Do not remove temporary files at the end (useful for debugging)")
            parser.add_argument('--verbose', action="store_true", default=False,
                                help="Show extra information on the log")
            parser.add_argument('--mapping-threads', default=8, help="Number of threads to use in the mapping step")
            parser.add_argument('--min-quality-trimming', default=20, help="Minimum quality for trimming")
            parser.add_argument('--bin-path', 
                                help="Path to folder where binary executables are present (system path by default)")
            parser.add_argument('--log-file', 
                                help="Name of the file that we want to use to store the logs (default output to screen)")
            parser.add_argument('--output-folder', help='Path of the output folder')
            parser.add_argument('--temp-folder', help='Path of the location for temporary files')
            parser.add_argument('--molecular-barcodes',
                                action="store_true", 
                                help="Activates the molecular barcodes PCR duplicates filter")
            parser.add_argument('--mc-allowed-mismatches', default=1,
                                help='Number of allowed mismatches when applying the molecular barcodes PCR filter')
            parser.add_argument('--mc-start-position', type=int, default=18,
                                help='Position (base wise) of the first base of the molecular barcodes (starting by 0)')
            parser.add_argument('--mc-end-position', default=27,
                                help='Position (base wise) of the last base of the molecular barcodes (starting by 1)')
            parser.add_argument('--min-cluster-size', default=2,
                                help='Min number of equal molecular barcodes to count as a cluster')
            parser.add_argument('--keep-discarded-files', action="store_true", default=False,
                                help='Writes down discarded reads and barcodes into files')
            parser.add_argument('--remove-polyA', default=0, 
                                help="Remove PolyAs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyT', default=0, 
                                help="Remove PolyTs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyG', default=0, 
                                help="Remove PolyGs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--remove-polyC', default=0,
                                help="Remove PolyCs and everything after it in the reads of a length at least as given number")
            parser.add_argument('--pair-mode-keep', default="reverse", 
                                help="When filtering out un-mapped reads if both strands map, what action to perform to keep (forward, reverse or both)")
            parser.add_argument('--filter-AT-content', default=90,
                                help="Discards reads whose number of A and T bases in total are more or equal than the number given in percentage")
            #TODO for now htseq-count does not support BAM files
            #http://sourceforge.net/p/htseq/bugs/12/
            #parser.add_argument('--sam-type', default="BAM",
            #                    help="Type of SAM format for intermediate files (SAM or BAM)")
            parser.add_argument('--disable-multimap', action="store_true", default=False,
                                help="If activated, multiple aligned reads obtained during mapping will be all discarded. Otherwise the highest scored one will be kept")
            parser.add_argument('--disable-clipping', action="store_true", default=False,
                                help="If activated, disable soft-clipping (local alignment) in the mapping")
            parser.add_argument('--mc-cluster', default="naive",
                                help="Type of clustering algorithm to use when performing UMIs duplicates removal.\n" \
                                "Modes = {naive(default), counttrie, hierarchical}")
            parser.add_argument('--min-intron-size', default=20,
                                help="Minimum allowed intron size when searching for splice variants in the mapping step")
            parser.add_argument('--max-intron-size', default=100000,
                                help="Maximum allowed intron size when searching for splice variants in the mapping step")
            parser.add_argument('--max-gap-size', default=1000000,
                                help="Maximum allowed distance between pairs in the mapping step")
            parser.add_argument('--umi-filter', action="store_true", default=False,
                                help="Enables the UMI quality filter based on the template given in --umi-filter-template")
            parser.add_argument('--umi-filter-template', default="WSNNWSNNV",
                                help="UMI template for the UMI filter, default = WSNNWSNNV")
            parser.add_argument('--compute-saturation', action="store_true", default=False,
                                help="Performs a saturation curve computation by sub-sampling the annotated reads, computing" \
                                "unique molecules and then a saturation curve")
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
        self.trimming_fw = int(options.mapping_fw_trimming)
        self.trimming_rv = int(options.mapping_rv_trimming)
        self.min_quality_trimming = int(options.min_quality_trimming) 
        self.clean = options.no_clean_up
        self.barcode_start = int(options.start_id)
        self.barcode_length = int(options.length_id)
        self.threads = int(options.mapping_threads)
        self.verbose = options.verbose
        if options.ids is not None:
            self.ids = os.path.abspath(options.ids)
        self.ref_map = os.path.abspath(options.ref_map)
        self.ref_annotation = os.path.abspath(options.ref_annotation)
        self.expName = options.expName
        self.htseq_mode = options.htseq_mode
        self.htseq_no_ambiguous = options.htseq_no_ambiguous
        self.qual64 = options.qual_64
        self.contaminant_index = options.contaminant_index
        self.path = options.bin_path
        if options.log_file is not None:
            self.logfile = os.path.abspath(options.log_file)  
        self.fastq_fw = options.fastq_files[0]
        self.fastq_rv = options.fastq_files[1]
        if options.output_folder is not None and os.path.isdir(options.output_folder):
            self.output_folder = os.path.abspath(options.output_folder)
        if options.temp_folder is not None and os.path.isdir(options.temp_folder): 
            self.temp_folder = os.path.abspath(options.temp_folder)
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
        self.pair_mode_keep = str(options.pair_mode_keep)
        self.filter_AT_content = int(options.filter_AT_content)
        #TODO for now htseq-count does not support BAM files
        #http://sourceforge.net/p/htseq/bugs/12/
        #self.sam_type = str(options.sam_type)
        self.disable_multimap = options.disable_multimap
        self.disable_clipping = options.disable_clipping
        self.mc_cluster = str(options.mc_cluster)
        self.min_intron_size = int(options.min_intron_size)
        self.max_intron_size = int(options.max_intron_size)
        self.max_gap_size = int(options.max_gap_size)
        self.umi_filter = options.umi_filter
        self.umi_filter_template = str(options.umi_filter_template).upper()
        self.compute_saturation = options.compute_saturation
        
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
        
        # Load the given path into the system PATH
        if self.path is not None and os.path.isdir(self.path): 
            os.environ["PATH"] += os.pathsep + self.path

        self.logger.info("ST Pipeline " + str(version_number))
        
        # Set default output and temp folders if erroneous given
        if self.output_folder is None or not os.path.isdir(self.output_folder):
            self.logger.info("Invalid path for output directory -- using current directory instead")
            self.output_folder = os.path.abspath(os.getcwd())
        if self.temp_folder is None or not os.path.isdir(self.temp_folder):
            self.logger.info("Invalid path for temp directory -- using current directory instead")
            self.temp_folder = os.path.abspath(os.getcwd())
        
        # Some info
        self.logger.info("Output directory : " + self.output_folder)
        self.logger.info("Temp directory : " + self.temp_folder)
        self.logger.info("Experiment : " + str(self.expName))
        self.logger.info("Forward reads file : " + str(self.fastq_fw))
        self.logger.info("Reverse reads file : " + str(self.fastq_rv))
        if self.ids is not None:
            self.logger.info("Ids file : " + str(self.ids))
        self.logger.info("Reference mapping file : " + str(self.ref_map))
        self.logger.info("Reference annotation file : " + str(self.ref_annotation))
        if self.contaminant_index is not None:
            self.logger.info("Using contamination filter with " + str(self.contaminant_index))
        self.logger.info("Nodes : " + str(self.threads))
        
        if self.ids is not None:
            self.logger.info("TaggD allowed mismatches " + str(self.allowed_missed))
            self.logger.info("TaggD barcode legnth " + str(self.barcode_length))
            self.logger.info("TaggD kmer size " + str(self.allowed_kmer))
            self.logger.info("TaggD overhang " + str(self.overhang))
        
        self.logger.info("Mapping forward trimming " + str(self.trimming_fw))
        self.logger.info("Mapping reverse trimming " + str(self.trimming_rv))
        self.logger.info("Mapper : STAR")
        self.logger.info("Annotation Tool : HTSeq")
        self.logger.info("When both strands are mapped, keep = " + str(self.pair_mode_keep))
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
            
        if self.molecular_barcodes and not self.ids is None:
            self.logger.info("Using Molecular Barcodes")
            self.logger.info("Molecular Barcode start position " + str(self.mc_start_position))
            self.logger.info("Molecular Barcode end position " + str(self.mc_end_position))
            self.logger.info("Molecular Barcode min cluster size " + str(self.min_cluster_size))
            self.logger.info("Molecular Barcode allowed mismatches " + str(self.mc_allowed_mismatches))
            self.logger.info("Molecular Barcode clustering algorithm " + str(self.mc_cluster))
            if self.umi_filter:
                self.logger.info("Molecular Barcodes using filter " + str(self.umi_filter_template))
                
        elif self.molecular_barcodes:
            self.logger.warning("Molecular barcodes cannot be used in normal RNA-Seq")
              
        if self.remove_polyA_distance > 0:
            self.logger.info("Removing polyA adaptors of a length at least " + str(self.remove_polyA_distance))                        
        if self.remove_polyT_distance > 0:
            self.logger.info("Removing polyT adaptors of a length at least " + str(self.remove_polyT_distance))
        if self.remove_polyG_distance > 0:
            self.logger.info("Removing polyG adaptors of a length at least " + str(self.remove_polyG_distance))
        if self.remove_polyC_distance > 0:
            self.logger.info("Removing polyC adaptors of a length at least " + str(self.remove_polyC_distance))

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
        if self.fastq_fw.endswith("gz"):
            temp_fastq_fw = os.path.join(self.temp_folder, "unzipped_fastq_fw.fastq")
            with gzip.open(self.fastq_fw, "rb") as filehandler_read:
                with open(temp_fastq_fw, "w") as filehandler_write:
                    filehandler_write.write(filehandler_read.read())
            self.fastq_fw = temp_fastq_fw
        if self.fastq_rv.endswith("gz"):
            temp_fastq_rv = os.path.join(self.temp_folder, "unzipped_fastq_rv.fastq")
            with gzip.open(self.fastq_rv, "rb") as filehandler_read:
                with open(temp_fastq_rv, "w") as filehandler_write:
                    filehandler_write.write(filehandler_read.read())
            self.fastq_rv = temp_fastq_rv
               
        #=================================================================
        # STEP: add BC and PolyT from FW reads to the RW reads and apply quality filter
        # also applies quality trimming and adaptor removal
        #=================================================================
        #NOTE after the trimming :
        #    - discarded reads will be replaced by Ns
        #    - reverse reads will contain the the 3 end the barcode and molecular barcode (if any)
        #      from the forward reads
        #    - if apply adaptors will be removed before the trimming
        self.logger.info("Start Reformatting and Filtering raw reads " + str(globaltime.getTimestamp()))
        fastq_fw_trimmed, fastq_rv_trimmed = reformatRawReads(self.fastq_fw,
                                                              self.fastq_rv,
                                                              self.qa_stats,
                                                              self.barcode_start,
                                                              self.barcode_length,
                                                              self.filter_AT_content,
                                                              self.molecular_barcodes,
                                                              self.mc_start_position,
                                                              self.mc_end_position,
                                                              self.trimming_fw,
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
          
        # Reverse reads contain now the barcodes and molecular barcodes (if any) so
        # the trimming needs to be adjusted for the mapping
        adjusted_forward_trimming = self.trimming_rv + self.barcode_length
        if self.molecular_barcodes:
            adjusted_forward_trimming += (self.mc_end_position - self.mc_start_position)  
           
        #=================================================================
        # CONDITIONAL STEP: Filter out contaminated reads, e.g. typically bacterial rRNA
        #=================================================================
        if self.contaminant_index:
            self.logger.info("rRNA filter")
            # To remove contaminants sequence with align the reads to the contaminant genome
            # and keep the un-mapped reads
            old_fw_trimmed = fastq_fw_trimmed
            old_rv_trimmed = fastq_rv_trimmed
            self.logger.info("Starting rRNA filter alignment " + str(globaltime.getTimestamp()))
            contaminated_sam, fastq_fw_trimmed, fastq_rv_trimmed = alignReads(fastq_fw_trimmed,
                                                                              fastq_rv_trimmed,
                                                                              self.contaminant_index,
                                                                              self.trimming_fw,
                                                                              adjusted_forward_trimming,
                                                                              self.threads,
                                                                              "rRNA_filter",
                                                                              self.min_intron_size,
                                                                              self.max_intron_size,
                                                                              self.max_gap_size,
                                                                              False, # disable splice variants alignments
                                                                              self.sam_type,
                                                                              False, # enable multimap in rRNA filter
                                                                              True, # disable softclipping in rRNA filter
                                                                              self.temp_folder)
            if self.clean:
                safeRemove(old_fw_trimmed)
                safeRemove(old_rv_trimmed)
            if not self.keep_discarded_files: safeRemove(contaminated_sam)
            
            # Unfortunately STAR will include all the reads in the output, including un-aligned reads
            # Therefore, we need to extract the un-aligned reads
            old_fw_trimmed = fastq_fw_trimmed
            old_rv_trimmed = fastq_rv_trimmed
            self.logger.info("Starting rRNA filter discarding unmapped reads " + str(globaltime.getTimestamp()))
            fastq_fw_trimmed, fastq_rv_trimmed = filter_rRNA_reads(fastq_fw_trimmed, 
                                                                   fastq_rv_trimmed,
                                                                   self.qa_stats, 
                                                                   self.temp_folder)
            safeRemove(old_fw_trimmed)
            safeRemove(old_rv_trimmed)
            
        #=================================================================
        # STEP: maps against the genome using STAR
        #=================================================================
        self.logger.info("Genome mapping")
        self.logger.info("Starting genome alignment " + str(globaltime.getTimestamp()))
        sam_mapped, unmapped_forward, unmapped_reverse = alignReads(fastq_fw_trimmed,
                                                                    fastq_rv_trimmed,
                                                                    self.ref_map,
                                                                    self.trimming_fw,
                                                                    adjusted_forward_trimming,
                                                                    self.threads,
                                                                    "alignment",
                                                                    self.min_intron_size,
                                                                    self.max_intron_size,
                                                                    self.max_gap_size,
                                                                    True, #enable splice variants alignments
                                                                    self.sam_type,
                                                                    self.disable_multimap,
                                                                    self.disable_clipping,
                                                                    self.temp_folder)
        if self.clean: 
            safeRemove(fastq_fw_trimmed)
            safeRemove(fastq_rv_trimmed)
        if not self.keep_discarded_files: 
            safeRemove(unmapped_forward)
            safeRemove(unmapped_reverse)
         
        #================================================================
        # STEP: filters mapped reads
        #================================================================
        self.logger.info("Starting genome alignment processing unmapped reads " + str(globaltime.getTimestamp()))
        sam_mapped_clean = filterMappedReads(sam_mapped,
                                             self.qa_stats,
                                             self.min_length_trimming,
                                             self.pair_mode_keep, 
                                             self.temp_folder, 
                                             self.keep_discarded_files)
        
        if self.clean: safeRemove(sam_mapped)
        
        #=================================================================
        # STEP: SORT sam file with mapped reads by read name
        #=================================================================
        #NO need to sort for now
        #sam_sorted = sortSamFile(sam_mapped_clean, self.temp_folder)
        #if self.clean: safeRemove(sam_mapped_clean)
        
        #=================================================================
        # STEP: annotate using htseq count
        #=================================================================
        self.logger.info("Starting annotation " + str(globaltime.getTimestamp()))
        annotatedFile = annotateReads(sam_mapped_clean, 
                                      self.ref_annotation,
                                      self.htseq_mode,
                                      self.temp_folder)
        if self.clean: safeRemove(sam_mapped_clean) 

        #=================================================================
        # STEP: filter out un-annotated reads
        #=================================================================
        self.logger.info("Starting filtering annotated reads " + str(globaltime.getTimestamp()))
        annotatedFilteredFile = filterAnnotatedReads(annotatedFile,
                                                     self.qa_stats, 
                                                     self.htseq_no_ambiguous,
                                                     self.temp_folder, 
                                                     self.keep_discarded_files)

        if self.clean: safeRemove(annotatedFile)
        
            
        if self.ids is not None:
            #=================================================================
            # STEP: Map against the barcodes
            #=================================================================
            self.logger.info("Starting barcode demultiplexing " + str(globaltime.getTimestamp()))
            mapFile = barcodeDemultiplexing(annotatedFilteredFile,
                                            self.ids,
                                            self.qa_stats,
                                            self.allowed_missed, 
                                            self.allowed_kmer, 
                                            self.barcode_start,  
                                            self.overhang,
                                            self.threads,
                                            self.temp_folder,
                                            self.keep_discarded_files)
            if self.clean: safeRemove(annotatedFilteredFile)
        
            if self.compute_saturation:
                self.logger.info("Starting computing saturation points " + str(globaltime.getTimestamp()))
                # #TODO move this to external functions
            
                import pysam
                import random
                import math
                #=================================================================
                # STEP: compute saturation curve
                #=================================================================
                
                # If we use count() to get number of reads from the PySAM object
                # it will require a seek() operation afterwards
                nreads = self.qa_stats.reads_after_annotation
                #saturation_points = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                saturation_points = list()
                for x in xrange(0,15):
                    spoint = int(math.floor(10e6 + (math.exp(x) * 10e6)))
                    if spoint >= nreads:
                        break
                    saturation_points.append(10e6 + (math.exp(x) * 10e6))

                files = dict()
                file_names = dict()
                subsampling = dict()
                file_ext = getExtension(mapFile).lower()
                flag_read = "rb"
                flag_write = "wbh"
                if file_ext == "sam":
                    flag_read = "r"
                    flag_write = "wh"
                
                annotated_sam = pysam.AlignmentFile(mapFile, flag_read)
            
                # Generate subsamples and files
                for spoint in saturation_points:
                    file_name = "subsample_" + str(spoint) + "." + file_ext
                    if self.temp_folder is not None and os.path.isdir(self.temp_folder):
                        file_name = os.path.join(self.temp_folder, file_name)
                    output_sam = pysam.AlignmentFile(file_name, flag_write, template=annotated_sam)
                    file_names[spoint] = file_name
                    files[spoint] = output_sam
                
                    indices = list(xrange(nreads))
                    random.shuffle(indices)
                    #amount = int(nreads * spoint)
                    amount = spoint
                    subbed = indices[0:amount]
                    subbed.sort()
                    subsampling[spoint] = subbed
            
                # Write subsamples
                index = 0
                sub_indexes = defaultdict(int)
                for read in annotated_sam:
                    for spoint in saturation_points:
                        if sub_indexes[spoint] == subsampling[spoint]:
                            files[spoint].write(read)
                            sub_indexes[spoint] += 1
                    index += 1  
            
                # Close the files
                annotated_sam.close()
                for file_sam in files.itervalues():
                    file_sam.close()
                
                # Compute saturation points
                saturation_points_values_unique_events = list()
                saturation_points_values_reads = list()
                for spoint in saturation_points:
                    stats = Stats()
                    input_file = file_names[spoint]
                    self.createDataset(input_file,
                                       stats,
                                       self.molecular_barcodes,
                                       self.mc_cluster,
                                       self.mc_allowed_mismatches,
                                       self.mc_start_position,
                                       self.mc_end_position,
                                       self.min_cluster_size,
                                       self.expName,
                                       False)
                    saturation_points_values_unique_events.append(stats.unique_events)
                    saturation_points_values_reads.append(stats.reads_after_duplicates_removal)
                
                if self.clean:
                    # Remove the files
                    for file_sam in file_names.itervalues():
                        safeRemove(file_sam)
                        
                # Generate plot
                self.logger.info("Saturation points:")
                self.logger.info(', '.join(str(a) for a in saturation_points))
                self.logger.info("Unique events per saturation point")
                self.logger.info(', '.join(str(a) for a in saturation_points_values_unique_events))
                self.logger.info("Reads per saturation point")
                self.logger.info(', '.join(str(a) for a in saturation_points_values_reads))
            
            #=================================================================
            # STEP: create json files with the results
            #=================================================================
            self.logger.info("Starting creating dataset " + str(globaltime.getTimestamp()))
            self.createDataset(mapFile,
                               self.qa_stats,
                               self.molecular_barcodes,
                               self.mc_cluster,
                               self.mc_allowed_mismatches,
                               self.mc_start_position,
                               self.mc_end_position,
                               self.min_cluster_size,
                               self.expName)
            if self.clean: safeRemove(mapFile)

        #=================================================================
        # END PIPELINE
        #=================================================================
        # Write stats to JSON
        self.qa_stats.writeJSON(os.path.join(self.output_folder, self.expName + "_qa_stats.json"))
        
        finish_exe_time = globaltime.getTimestamp()
        total_exe_time = finish_exe_time - start_exe_time
        self.logger.info("Total Execution Time : " + str(total_exe_time))

    def createDataset(self,
                      input_name,
                      qa_stats,
                      molecular_barcodes = False,
                      mc_cluster = "naive",
                      allowed_mismatches = 1,
                      start_position = 18, 
                      end_position = 27, 
                      min_cluster_size = 2,
                      output_template = None,
                      verbose = True):
        """ 
        parse annotated and mapped reads with the reads that contain barcodes to
        create json files with the barcodes and coordinates and json file with the raw reads
        and some useful stats and plots
        It also allows to remove PCR Duplicates using molecular barcodes
        We passes the number of forward bases trimmed for mapping to get a clean read
        in the output
        """
        args = ['createDataset.py', '--input', str(input_name)]
        
        if molecular_barcodes:
            args += ['--molecular-barcodes',
                     '--mc-allowed-mismatches', allowed_mismatches,
                     '--mc-start-position', start_position, 
                     '--mc-end-position', end_position, 
                     '--min-cluster-size', min_cluster_size,
                     '--mc-cluster', mc_cluster]
            
        if self.output_folder is not None:
            args += ['--output-folder', self.output_folder]
     
        if output_template is not None:
            args += ['--output-file-template', output_template]
            
        try:
            proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (stdout, errmsg) = proc.communicate()
        except Exception as e:
            error = "Error creating dataset: createDataset execution failed"
            self.logger.error(error)
            self.logger.error(e)
            raise
        
        if len(errmsg) > 0:
            error = "Error, There was an error creating the dataset: \n" + stdout + " \n" + errmsg
            self.logger.error(error)
            raise RuntimeError(error)    
              
        procOut = stdout.split("\n")
        for line in procOut:
            # Write QA stats
            # TODO a more efficient way perhaps to use the line numbers 
            if line.find("Number of Transcripts with Barcode present:") != -1:
                qa_stats.reads_after_duplicates_removal = int(line.split()[-1])
            if line.find("Number of unique events present:") != -1:
                qa_stats.unique_events = int(line.split()[-1])
            if line.find("Number of unique Barcodes present:") != -1:
                qa_stats.genes_found = int(line.split()[-1])
            if line.find("Number of unique Genes present:") != -1:
                qa_stats.barcodes_found = int(line.split()[-1])
            if line.find("Number of discarded reads (possible PCR duplicates):") != -1:
                qa_stats.duplicates_found = int(line.split()[-1])
            if line.find("Max number of genes over all features:") != -1:
                qa_stats.max_genes_feature = int(line.split()[-1])
            if line.find("Min number of genes over all features:") != -1:
                qa_stats.min_genes_feature = int(line.split()[-1])
            if line.find("Max number of reads over all features:") != -1:
                qa_stats.max_reads_feature = int(line.split()[-1])
            if line.find("Min number of reads over all features:") != -1:
                qa_stats.min_reads_feature = int(line.split()[-1])
            if line.find("Max number of reads over all unique events:") != -1:
                qa_stats.max_reads_unique_event = int(line.split()[-1])
            if line.find("Min number of reads over all unique events:") != -1:
                qa_stats.min_reads_unique_event = int(line.split()[-1])
            if verbose:   
                self.logger.info(str(line))
