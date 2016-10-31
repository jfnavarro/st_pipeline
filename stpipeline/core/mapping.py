""" 
This module contains functions related to sequence alignment and barcode
demultiplexing in the ST pipeline
"""
import logging 
import subprocess
from subprocess import CalledProcessError
from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import uuid
import os
import shutil

def createIndex(genome,
                log_sj,
                threads,
                tmpFolder):
    """
    When running STAR in two pass mode a new index needs to be created
    using the discovered splice variants. This functions generates
    the index and returns the location to the new index.
    :param genome: the path to the fasta file with the original genome
    :param threads: number of threads to be used
    :param tmpFolder: where to place the temporary files and to find the file with variants
    :param log_sj: the path to the file containing the splices
    :type genome: str
    :type threads: int
    :type tmpFolder: str
    :type log_sj: str
    :return: the path to the new index
    :raises: RuntimeError,ValueError,OSError,CalledProcessError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(log_sj):
        error = "Error creating index with STAR.\n" \
        "input file not present {}\n".format(log_sj)
        logger.error(error)
        raise RuntimeError(error)
    
    # Unique identifier to the index folder
    genome_dir = str(uuid.uuid4())
    log = "Log.out"
    log_final = "Log.final.out"
    if tmpFolder is not None and os.path.isdir(tmpFolder):
        genome_dir = os.path.join(tmpFolder, genome_dir)
        log_final = os.path.join(tmpFolder, log_final)
        log = os.path.join(tmpFolder, log)
    
    os.mkdir(genome_dir)
    if not os.path.isdir(genome_dir):
        error = "Error creating index with STAR.\n "
        "There was a problem creating a temp folder to place the index\n"
        logger.error(error)
        raise RuntimeError(error)   
    
    args = ['STAR']
    args += ["--runMode", "genomeGenerate",
             "--genomeDir", genome_dir,
             "--genomeFastaFiles", genome,
             "--sjdbOverhang", 100,
             "--runThreadN", threads,
             "--sjdbFileChrStartEnd", log_sj]
    try:
        proc = subprocess.Popen([str(i) for i in args],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                close_fds=True, shell=False)
        (stdout, errmsg) = proc.communicate()
    except ValueError as e:
        logger.error("Error creating index with STAR\n Incorrect arguments.")
        raise e
    except OSError as e:
        logger.error("Error creating index with STAR\n Executable not found.")
        raise e
    except CalledProcessError as e:
        logger.error("Error creating index with STAR\n Program returned error.")
        raise e
    
    if os.path.isfile(log_sj): os.remove(log_sj)
    if os.path.isfile(log): os.remove(log)
    if os.path.isfile(log_final): os.remove(log_final)
    
    if len(errmsg) > 0:
        logger.warning("STAR outputted an error message while " \
                       "creating the index.\n{}\n".format(errmsg))
    
    if not os.path.isdir(genome_dir):
        error = "Error creating index with STAR.\nThe index output folder " 
        "is not present\n{}\n".format(errmsg)
        logger.error(error)
        raise RuntimeError(error)
    
    return genome_dir

def alignReads(reverse_reads, 
               ref_map,
               outputFile,
               outputFileDiscarded=None,
               outputFolder=None,
               trimReverse=0, 
               cores=4,
               min_intron_size=20,
               max_intron_size=1000000,
               max_gap_size=1000000,
               use_splice_juntions=True,
               disable_multimap=False,
               diable_softclipping=False,
               invTrimReverse=0,
               sortedBAMOutput=True):
    """
    This function will perform a sequence alignment using STAR.
    Mapped and unmapped reads are written to the paths given as
    parameters. It needs the path of the STAR genome index. 
    :param reverse_reads: file containing reverse reads in fastq format (Illumina pair end)
    :param ref_map: a path to the genome/transcriptome STAR index
    :param outputFile: the name of the SAM/BAM output file to write the alignments to
    :param outputFileDiscarded: the name of the SAM/BAM output file to write discarded alignments
    :param outputFolder: the path of the output folder
    :param trimReverse: the number of bases to trim in the reverse reads (to not map)
    :param cores: the number of cores to use to speed up the alignment
    :param file_name_patter: indicates how the output files will be named
    :param min_intron_size: min allowed intron size when spanning splice junctions
    :param max_intron size: max allowed intron size when spanning splice junctions
    :param max_gap_size: max allowed gap between pairs
    :param use_splice_junctions: whether to use splice aware alignment or not
    :param disable_multimap: if True no multiple alignments will be allowed
    :param diable_softclipping: it True no local alignment allowed
    :param invTrimReverse: number of bases to trim in the 5' of the read2
    :param sortedBAMOutput: True if the BAM output must be sorted
    :type reverse_reads: str
    :type ref_map: str
    :type outputFile: str
    :type outputFileDiscarded: str
    :type outputFolder: str
    :type trimReverse: int
    :type cores: int
    :type file_name_patter: str
    :type min_intron_size: int
    :type max_intron: int
    :type max_gap_size: int
    :type use_splice_junctions: bool
    :type disable_multimap: bool
    :type diable_softclipping: bool
    :type invTrimReverse: int
    :type sortedBAMOutput: bool
    :raises: RuntimeError,ValueError,OSError,CalledProcessError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(reverse_reads):
        error = "Error mapping with STAR, input file not present {}\n".format(reverse_reads)
        logger.error(error)
        raise RuntimeError(error)
    
    # STAR has predefined output names for the files
    tmpOutputFile = "Aligned.sortedByCoord.out.bam" if sortedBAMOutput else "Aligned.out.bam"
    tmpOutputFileDiscarded = "Unmapped.out.mate1"
    log_std = "Log.std.out"
    log = "Log.out"
    log_sj = "SJ.out.tab"
    log_final = "Log.final.out"
    log_progress = "Log.progress.out"
    
    if outputFolder is not None and os.path.isdir(outputFolder):
        tmpOutputFile = os.path.join(outputFolder, tmpOutputFile)
        tmpOutputFileDiscarded = os.path.join(outputFolder, tmpOutputFileDiscarded)
        log_std = os.path.join(outputFolder, log_std)
        log = os.path.join(outputFolder, log)
        log_sj = os.path.join(outputFolder, log_sj)
        log_final = os.path.join(outputFolder, log_final)
        log_progress = os.path.join(outputFolder, log_progress)
    
    # Options
    # outFilterType(BySJout) this will keep only reads 
    #     that contains junctions present in SJ.out.tab
    # outSamOrder(Paired) one mate after the other 
    # outSAMprimaryFlag(OneBestScore) only one alignment with the best score is primary
    # outFilterMultimapNmax 
    #     read alignments will be output only if the read maps fewer than this value
    # outFilterMismatchNmax = alignment will be output only if 
    #     it has fewer mismatches than this value
    # outFilterMismatchNoverLmax = alignment will be output only if 
    #     its ratio of mismatches to *mapped* length is less than this value
    # alignIntronMin minimum intron size: genomic gap is considered intron 
    #     if its length>=alignIntronMin, otherwise it is considered Deletion
    # alignIntronMax maximum intron size, if 0, max intron size will be 
    #     determined by (2 to the power of winBinNbits)*winAnchorDistNbins
    # alignMatesGapMax maximum gap between two mates, if 0, max intron gap will 
    #     be determined by (2 to the power of winBinNbits)*winAnchorDistNbins
    # alignEndsType Local standard local alignment with soft-clipping allowed EndToEnd: 
    #     force end-to-end read alignment, do not soft-clip
    # chimSegmentMin if >0 To switch on detection of chimeric (fusion) alignments
    # --outMultimapperOrder Random multimap are written in Random order
    # --outSAMmultNmax Number of multimap that we want to output 
    # put to 1 to not include multiple mappings (default 10)
    
    multi_map_number = 1 if disable_multimap else 10
    alignment_mode = "EndToEnd" if diable_softclipping else "Local"
    sjdb_overhang = 100 if use_splice_juntions else 0
    bam_sorting = "SortedByCoordinate" if sortedBAMOutput else "Unsorted"
    
    core_flags = ["--runThreadN", str(max(cores, 1))]
    trim_flags = ["--clip5pNbases", trimReverse, 
                  "--clip3pNbases", invTrimReverse]
    io_flags   = ["--outFilterType", "Normal", 
                  "--outSAMtype", "BAM", bam_sorting, 
                  "--alignEndsType", alignment_mode, 
                  "--outSAMunmapped", "None", # unmapped reads not included in main output
                  "--outSAMorder", "Paired",    
                  "--outSAMprimaryFlag", "OneBestScore", 
                  "--outFilterMultimapNmax", multi_map_number, 
                  "--alignSJoverhangMin", 5, # default is 5
                  "--alignSJDBoverhangMin", 3, # default is 3
                  "--sjdbOverhang", sjdb_overhang, # 0 to not use splice junction database
                  "--outFilterMismatchNmax", 10, # large number switches it off (default 10)
                  "--outFilterMismatchNoverLmax", 0.3, # default is 0.3
                  "--alignIntronMin", min_intron_size,
                  "--alignIntronMax", max_intron_size, 
                  "--alignMatesGapMax", max_gap_size,
                  "--winBinNbits", 16,
                  "--winAnchorDistNbins", 9,
                  "--chimSegmentMin", 0,
                  "--readMatesLengthsIn", "NotEqual",
                  "--genomeLoad", "NoSharedMemory"] # Options to use share remove can be given here 

    args = ['STAR']
    args += trim_flags
    args += core_flags
    args += io_flags
    args += ["--genomeDir", ref_map,
             "--readFilesIn", reverse_reads,
             "--outFileNamePrefix", outputFolder + os.sep]  # MUST ENSURE AT LEAST ONE SLASH
    args += ["--outReadsUnmapped", "Fastx"]
    
    try:
        proc = subprocess.Popen([str(i) for i in args],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                close_fds=True, shell=False)
        (stdout, errmsg) = proc.communicate()
    except ValueError as e:
        logger.error("Error mapping with STAR\n Incorrect arguments.")
        raise e
    except OSError as e:
        logger.error("Error mapping with STAR\n Executable not found.")
        raise e
    except CalledProcessError as e:
        logger.error("Error mapping with STAR\n Program returned error.")
        raise e
        
    if not fileOk(tmpOutputFile):
        error = "Error mapping with STAR.\n" \
        "Output file not present {}\n{}\n".format(tmpOutputFile, errmsg)
        logger.error(error)
        raise RuntimeError(error)

    if len(errmsg) > 0:
        logger.warning("STAR has generated error messages during mapping.\n{}\n".format(errmsg))
        
    # Rename output files.
    shutil.move(tmpOutputFile, outputFile)
    shutil.move(tmpOutputFileDiscarded, outputFileDiscarded)
        
    # Remove temp files from STAR
    if os.path.isfile(log_std): os.remove(log_std)
    if os.path.isfile(log): os.remove(log)
    if os.path.isfile(log_progress): os.remove(log_progress)
    # Do not remove to use it for computing a new index in 2pass mode
    # if os.path.isfile(log_sj): os.remove(log_sj)
    
    if not os.path.isfile(log_final):
        logger.warning("Log output file from STAR is not present")
    else:
        logger.info("Mapping stats: ")
        logger.info("Mapping stats are computed from all the pair reads present in the raw files")
        uniquely_mapped = 0
        multiple_mapped = 0
        # Parse log file from STAR to get stats
        # TODO find a cleaner way to do this
        with open(log_final, "r") as star_log:
            for line in star_log.readlines():
                if line.find("Uniquely mapped reads %") != -1 \
                or line.find("Uniquely mapped reads number") != -1 \
                or line.find("Number of reads mapped to multiple loci") != -1 \
                or line.find("% of reads mapped to multiple loci") != -1 \
                or line.find("% of reads unmapped: too short") != -1:
                    logger.info(str(line).rstrip())
                # Some duplicated code here; TODO refactor
                if line.find("Uniquely mapped reads number") != -1:
                    uniquely_mapped = int(str(line).rstrip().split()[-1])
                if line.find("Number of reads mapped to multiple loci") != -1:
                    multiple_mapped = int(str(line).rstrip().split()[-1])
            logger.info("Total mapped reads: {}".format(uniquely_mapped + multiple_mapped))   
             
    # Remove log file       
    if os.path.isfile(log_final): os.remove(log_final)

def barcodeDemultiplexing(reads, 
                          idFile,
                          mismatches,
                          kmer, 
                          start_positon,
                          over_hang,
                          cores,
                          outputFilePrefix,
                          keep_discarded_files=False):
    """ 
    This functions performs a demultiplexing using Taggd. Input reads will be filtered
    out looking at their barcodes. Only the ones that contain a barcode
    that is matched in the barcodes files will be kept.
    Information about the barcode and the array coordinates will be added
    to the output file. 
    :param reads: a file in SAM/BAM/FASTQ format containing reads with barcodes
    :param idFile: a tab delimited file (BARCODE - X - Y) containing all the barcodes
    :param mismatches: the number of allowed mismatches
    :param kmer: the kmer length
    :param start_positon: the start position of the barcode
    :param over_hang: the number of bases to allow for overhang
    :param outputFilePrefix: location and prefix for the output files
    :param keep_discarded_files: if True files with the non demultiplexed reads will be generated
    :type reads: str
    :type idFile: str
    :type mismatches: int
    :type kmer: int
    :type start_positon: int
    :type over_hang: int
    :type outputFilePrefix: str
    :type keep_discarded_files: bool
    :raises: RuntimeError,ValueError,OSError,CalledProcessError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(reads):
        error = "Error, input file not present {}\n".format(reads)
        logger.error(error)
        raise RuntimeError(error)
    
    # Taggd options
    #--metric (subglobal (default) , Levenshtein or Hamming)
    #--slider-increment (space between kmer searches, 0 is default = kmer length)
    #--seed
    #--no-multiprocessing
    #--overhang additional flanking bases around read barcode to allow
    #--estimate-min-edit-distance is set estimate the min edit distance among true barcodes
    #--no-offset-speedup turns off speed up, 
    #  it might yield more hits (exactly as findIndexes)
    #--homopolymer-filter if set excludes erads where barcode 
    #  contains a homolopymer of the given length (0 no filter), default 8
    args = ['taggd_demultiplex.py',
            "--max-edit-distance", mismatches,
            "--k", kmer,
            "--start-position", start_positon,
            "--homopolymer-filter", 8,
            "--subprocesses", cores,
            "--overhang", over_hang]
    
    if not keep_discarded_files:
        args.append("--no-unmatched-output")
        args.append("--no-ambiguous-output")
        args.append("--no-results-output")
        
    args += [idFile, reads, outputFilePrefix]

    try:
        proc = subprocess.Popen([str(i) for i in args], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                close_fds=True, shell=False)
        (stdout, errmsg) = proc.communicate()
    except ValueError as e:
        logger.error("Error demultiplexing with TAGGD\n Incorrect arguments.")
        raise e
    except OSError as e:
        logger.error("Error demultiplexing with TAGGD\n Executable not found.")
        raise e
    except CalledProcessError as e:
        logger.error("Error demultiplexing with TAGGD\n Program returned error.")
        raise e
    
    # We know the output file from the prefix and suffix
    outputFile = "{}_matched{}".format(outputFilePrefix, os.path.splitext(reads)[1].lower())
    if not fileOk(outputFile):
        error = "Error demultiplexing with TAGGD.\n" \
        "Output file is not present {}\n{}\n".format(outputFile, errmsg)
        logger.error(error)
        raise RuntimeError(error)
 
    if len(errmsg) > 0:
        logger.warning("Taggd has generated error messages during " \
                       "demultiplexing.\n{}\n".format(errmsg))
           
    # TODO must be a cleaner way to get the stats from the output file
    procOut = stdout.split("\n")
    logger.info("Barcode Mapping stats:")
    for line in procOut: 
        if line.find("Total reads:") != -1:
            logger.info(str(line))
        if line.find("Total reads written:") != -1:
            # Update the QA stats
            # TODO find a cleaner way to to this
            qa_stats.reads_after_demultiplexing = int(line.split()[-1])
            logger.info(str(line))
        if line.find("Perfect Matches:") != -1:
            logger.info(str(line))
        if line.find("Imperfect Matches") != -1:
            logger.info(str(line))
        if line.find("Ambiguous matches:") != -1:
            logger.info(str(line))
        if line.find("Non-unique ambiguous matches:") != -1:
            logger.info(str(line))
        if line.find("Unmatched:") != -1:
            logger.info(str(line))