""" 
This module contains functions related to sequence alignment and barcode
demultiplexing in the ST pipeline
"""
import logging 
import subprocess
from subprocess import CalledProcessError
from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import os
import shutil

def alignReads(reverse_reads, 
               ref_map,
               outputFile,
               annotation,
               outputFileDiscarded,
               outputFolder,
               trimReverse,
               invTrimReverse,
               cores,
               min_intron_size,
               max_intron_size,
               disable_multimap,
               diable_softclipping,
               twopassMode,
               min_length):
    """
    This function will perform a sequence alignment using STAR.
    Mapped and unmapped reads are written to the paths given as
    parameters. It needs the path of the STAR genome index.
    It allows to perform the 2-Pass mode.
    It needs the annotation file to use the on-the-fly mode.
    :param reverse_reads: file containing reverse reads in fastq format (Illumina pair end)
    :param ref_map: a path to the genome/transcriptome STAR index
    :param outputFile: the name of the SAM/BAM output file to write the alignments to
    :param annotation: the annotation file in GTF
    :param outputFileDiscarded: the name of the SAM/BAM output file to write discarded alignments
    :param outputFolder: the path of the output folder
    :param trimReverse: the number of bases to trim in the reverse reads (from 5')
    :param invTrimReverse: number of bases to trim from the 3'
    :param cores: the number of cores to use to speed up the alignment
    :param min_intron_size: min allowed intron size when spanning splice junctions
    :param max_intron size: max allowed intron size when spanning splice junctions 
    :param disable_multimap: if True no multiple alignments will be allowed
    :param diable_softclipping: it True no local alignment allowed
    :param twopassMode: True to use the 2-pass mode
    :param min_length: the min allowed read length (mapped bases)
    :type reverse_reads: str
    :type ref_map: str
    :type outputFile: str
    :type annotation: str
    :type outputFileDiscarded: str
    :type outputFolder: str
    :type trimReverse: int
    :type invTrimReverse: int
    :type cores: int
    :type min_intron_size: int
    :type max_intron: int
    :type disable_multimap: bool
    :type diable_softclipping: bool
    :type twopassMode: bool
    :type min_length: str
    :raises: RuntimeError,ValueError,OSError,CalledProcessError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(reverse_reads):
        error = "Error mapping with STAR, input file not present {}\n".format(reverse_reads)
        logger.error(error)
        raise RuntimeError(error)
    
    # STAR has predefined output names for the files
    tmpOutputFile = "Aligned.sortedByCoord.out.bam"
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
    
    multi_map_number = 1 if disable_multimap else 20 # 10 is the STAR default
    alignment_mode = "EndToEnd" if diable_softclipping else "Local"
    
    flags = ["--clip3pNbases", invTrimReverse,
             "--clip5pNbases", trimReverse,
             "--runThreadN", str(max(cores, 1)),
             "--outFilterType", "Normal", 
             "--outSAMtype", "BAM", "SortedByCoordinate",
             "--alignEndsType", alignment_mode, 
             "--outSAMunmapped", "None", # unmapped reads not included in main output
             "--outSAMorder", "Paired",    
             "--outSAMprimaryFlag", "OneBestScore", 
             "--outFilterMultimapNmax", multi_map_number,
             "--alignIntronMin", min_intron_size,
             "--alignIntronMax", max_intron_size,
             "--outFilterMatchNmin", min_length,
             "--outSAMmultNmax", 1,
             "--readMatesLengthsIn", "NotEqual",
             "--outFilterMismatchNoverLmax", 0.1, ## (0.3 default)
             "--genomeLoad", "NoSharedMemory"] 
    
    if twopassMode:
        flags += ["--twopassMode", "Basic"]

    if annotation is not None:
        flags += ["--sjdbGTFfile", annotation]
        
    args = ["STAR",
            "--genomeDir", ref_map,
            "--readFilesIn", reverse_reads,
            "--outFileNamePrefix", outputFolder + os.sep, # MUST ENSURE AT LEAST ONE SLASH
            "--outReadsUnmapped", "Fastx"]  
    args += flags
    
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
    if os.path.isfile(log_sj): os.remove(log_sj)
    
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
                          taggd_metric,
                          taggd_multiple_hits_keep_one,
                          taggd_trim_sequences,
                          cores,
                          outputFilePrefix,
                          keep_discarded_files=False):
    """ 
    This functions performs a demultiplexing using Taggd. Input reads will be filtered
    out looking at their barcodes. Only the ones that contain a barcode
    that is matched in the barcodes files will be kept.
    Information about the barcode and the array coordinates will be added
    to the output file. 
    :param reads: a file in FASTQ/BAM format containing reads with barcodes
    :param idFile: a tab delimited file (BARCODE - X - Y) containing all the barcodes
    :param mismatches: the number of allowed mismatches
    :param kmer: the kmer length
    :param start_positon: the start position of the barcode
    :param over_hang: the number of bases to allow for overhang
    :param taggd_metric: the distance metric algorithm (Subglobal, Levensthein or Hamming)
    :param taggd_multiple_hits_keep_one: when True keep one random hit when multiple candidates
    :param taggd_trim_sequences: coordinates to trim in the barcode
    :param outputFilePrefix: location and prefix for the output files
    :param keep_discarded_files: if True files with the non demultiplexed reads will be generated
    :type reads: str
    :type idFile: str
    :type mismatches: int
    :type kmer: int
    :type start_positon: int
    :type over_hang: int
    :type taggd_metric: str
    :type taggd_multiple_hits_keep_one: bool
    :type taggd_trim_sequences: list
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
    #--overhang additional flanking bases around read barcode to allow
    #--estimate-min-edit-distance is set estimate the min edit distance among true barcodes
    #--no-offset-speedup turns off speed up, 
    #  it might yield more hits (exactly as findIndexes)
    #--homopolymer-filter if set excludes reads where barcode 
    #  contains a homolopymer of the given length (0 no filter), default 8
    
    if taggd_metric == "Hamming": over_hang = 0 
    args = ['taggd_demultiplex.py']
    
    if taggd_trim_sequences is not None:
        args.append("--trim-sequences") 
        for pos in taggd_trim_sequences:
            args.append(pos) 
            
    args += ["--max-edit-distance", mismatches,
            "--k", kmer,
            #"--barcode-tag", "B0", # if input if BAM we tell taggd what tag contains the barcode
            "--start-position", start_positon,
            "--homopolymer-filter", 0,
            "--subprocesses", cores,
            "--metric", taggd_metric,
            "--overhang", over_hang]
            
    if taggd_multiple_hits_keep_one:
        args.append("--multiple-hits-keep-one")  
            
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
    logger.info("Demultiplexing Mapping stats:")
    for line in procOut: 
        if line.find("Total reads:") != -1:
            logger.info(str(line))
        if line.find("Total reads written:") != -1:
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