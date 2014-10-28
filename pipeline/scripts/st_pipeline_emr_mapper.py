#!/usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""

''' This is the Amazon EMR wrapper for the mapper function to run the Map Reduce version of the ST pipeline'''

import sys
sys.path.append("./")
from stpipeline.core.pipeline import *
from stpipeline.common.json_utils import json_iterator
import argparse
                  
def main(argv):
    
    #TODO find a way to retrieve this from st_pipeline_run
    # code is duplicated
 
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fastq_files', nargs=2)
    parser.add_argument('--chunks', default=10000, help='Number of reads per chunk to split (Hadoop Map Reduced)')
    parser.add_argument('--output-file', default="merged.json", help='Name of the output file (Hadoop Map Reduced)')
    parser.add_argument('--ids', help='The name of the file containing the barcodes and the coordinates')
    parser.add_argument('--ref-map', help= "<path_to_bowtie2_indexes> = Reference genome name for the genome that you want to use to align the reads")
    parser.add_argument('--ref-annotation', help="Path to the reference annotation file (htseq requires a GTF file annotation file that you want to use to annotate")
    parser.add_argument('--expName', help="Name of the experiment (output file name)")
    parser.add_argument('--allowed-missed', default=6, help="Number of allowed mismatches when mapping against the barcodes")
    parser.add_argument('--allowed-kimer', default=7, help="KMer length when mapping against the barcodes")
    parser.add_argument('--min-length-qual-trimming', default=28, help="Minimum length of the sequence for mapping after trimming, shorter reads will be discarded")
    parser.add_argument('--mapping-fw-trimming', default=42, help="Number of bases to trim in the forward reads for the Mapping [24 + ID_LENGTH]")
    parser.add_argument('--mapping-rv-trimming', default=5, help="Number of bases to trim in the reverse reads for the Mapping")
    parser.add_argument('--length-id', default=18, help="Length of ID, a.k.a. the length of the barcodes")
    parser.add_argument('--contaminant-bowtie2-index', help="<path_to_bowtie2_indexes> = When provided, reads will be filtered against the specified bowtie2 index, non-mapping reads will be saved and demultiplexed")
    parser.add_argument('--qual-64', action="store_true", default=False, help="Use phred-64 quality instead of phred-33(default)")
    parser.add_argument('--htseq-mode', default="intersection-nonempty", help="Mode of Annotation when using HTSeq. Modes = {union,intersection-nonempty(default),intersection-strict}")
    parser.add_argument('--htseq-no-ambiguous', action="store_true", help="When using htseq discard reads annotating ambiguous genes")
    parser.add_argument('--start-id', default=0, help="Start position of BARCODE ID [0]")
    parser.add_argument('--error-id', default=0, help="Id positional error [0]")
    parser.add_argument('--no-clean-up', action="store_false", default=True, help="Do not remove temporary files at the end (useful for debugging)")
    parser.add_argument('--verbose', action="store_true", default=False, help="Show extra information on the log")
    parser.add_argument('--bowtie-threads', default=8, help="Number of threads to use in the mapping step")
    parser.add_argument('--min-quality-trimming', default=20, help="Minimum quality for trimming")
    parser.add_argument('--discard-fw', action="store_true", default=False, help="Discard forwards reads that maps uniquely")
    parser.add_argument('--discard-rv', action="store_true", default=False, help="Discard reverse reads that maps uniquely")
    parser.add_argument('--bowtie2-discordant', action="store_true", default=False, help="Discard non-discordant alignments when mapping")
    parser.add_argument('--bin-path', help="Path to folder where binary executables are present (system path by default)")
    parser.add_argument('--log-file', help="Name of the file that we want to use to store the logs (default output to screen)")
    parser.add_argument('--output-folder', help='Path of the output folder')
    parser.add_argument('--temp-folder', help='Path of the location for temporary files')
    parser.add_argument('--molecular-barcodes', 
                        action="store_true", help="Activates the molecular barcodes PCR duplicates filter")
    parser.add_argument('--mc-allowed-missmatches', default=3,
                        help='Number of allowed missmatches when applying the molecular barcodes PCR filter')
    parser.add_argument('--mc-start-position', type=int, default=19,
                        help='Position (base wise) of the first base of the molecular barcodes')
    parser.add_argument('--mc-end-position', default=28,
                        help='Position (base wise) of the last base of the molecular barcodes')
    parser.add_argument('--min-cluster-size', default=10,
                        help='Min number of equal molecular barcodes to count as a cluster')
    
    #parse arguments
    options = parser.parse_args()
    
    #local variables
    batch = []
    chunks = [] 
    
    pipeline = Pipeline()  
    #init pipeline arguments
    pipeline.allowed_missed = int(options.allowed_missed)
    pipeline.allowed_kimera = int(options.allowed_kimer)
    pipeline.min_length_trimming = int(options.min_length_qual_trimming)
    pipeline.trimming_fw_bowtie = int(options.mapping_fw_trimming)
    pipeline.trimming_rw_bowtie = int(options.mapping_rv_trimming)
    pipeline.min_quality_trimming = int(options.min_quality_trimming) 
    pipeline.clean = options.no_clean_up
    pipeline.s = int(options.start_id)
    pipeline.l = int(options.length_id)
    pipeline.e = int(options.error_id)
    pipeline.threads = int(options.bowtie_threads)
    pipeline.verbose = options.verbose
    pipeline.ids = os.path.abspath(options.ids)
    pipeline.ref_map = os.path.abspath(options.ref_map)
    pipeline.ref_annotation = os.path.abspath(options.ref_annotation)
    pipeline.expName = options.expName
    pipeline.htseq_mode = options.htseq_mode
    pipeline.htseq_no_ambiguous = options.htseq_no_ambiguous
    pipeline.qual64 = options.qual_64
    pipeline.discard_fw = options.discard_fw
    pipeline.discard_rv = options.discard_rv
    pipeline.discordant = options.bowtie2_discordant
    pipeline.contaminant_bt2_index = options.contaminant_bowtie2_index
    pipeline.path = options.bin_path
    if options.log_file is not None:
        pipeline.logfile = os.path.abspath(options.log_file)         
    if options.output_folder is not None and os.path.isdir(options.output_folder):
        pipeline.output_folder = os.path.abspath(options.output_folder)
    if options.temp_folder is not None and os.path.isdir(options.temp_folder): 
        pipeline.temp_folder = os.path.abspath(options.temp_folder)
    pipeline.molecular_barcodes = options.molecular_barcodes
    pipeline.mc_allowed_missmatches = int(options.mc_allowed_missmatches)
    pipeline.mc_start_position = int(options.mc_start_position)
    pipeline.mc_end_position = int(options.mc_end_position)
    pipeline.min_cluster_size = int(options.min_cluster_size)
        
    #test and load parameters
    pipeline.load_parameters()
    
    #parse input streaming to create chunks
    line = sys.stdin.readline()
    try:
        while line:
            batch.append(line + '\n')
            if len(batch) == int(options.chunks): 
                chunks.append(''.join(batch))
                batch = []
            line = sys.stdin.readline()
    except "end of file":
        return None
    
    #some lines left behind?
    if len(batch) > 0:
        chunks.append(''.join(batch))
        
    #call pipeline and output streams
    for key,value in pipeline.run_pipeline(chunks):
        print '%s%s%d' % (key, "\t", value)
    
if __name__ == "__main__":
    main(sys.argv)
