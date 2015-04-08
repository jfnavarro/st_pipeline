#!/usr/bin/env python
""" 
This module contains functions related to sequence alignment and barcode
demultiplexing in the ST pipeline
"""

import logging 
import subprocess
import os
from stpipeline.common.utils import *

def alignReads(forward_reads, 
               reverse_reads, 
               ref_map, 
               trimForward, 
               trimReverse, 
               cores, 
               outputFolder=None):
    """
    :param forward_reads file containing forward reads in fastq format for pair end sequences
    :param reverse_reads file containing reverse reads in fastq format for pair end sequences
    :param ref_map a path to the genome/transcriptome indexes
    :param trimForward the number of bases to trim in the forward reads (to not map)
    :param trimReverse the number of bases to trim in the reverse reasd (to not map)
    :param cores the number of cores to use to speed up the alignment
    :param outputFolder if set all output files will be placed there
    This function will perform a sequence alignement using STAR
    mapped and unmapped reads are returned and the set of parameters
    used are described
    """

    logger = logging.getLogger("STPipeline")
    logger.info("Start STAR Mapping")
    
    # STAR has predefined output names
    tmpOutputFile = "Aligned.out.sam"
    tmpOutputFileDiscarded1 = "Unmapped.out.mate1"
    tmpOutputFileDiscarded2 = "Unmapped.out.mate2"
    outputFile = replaceExtension(getCleanFileName(forward_reads),"_mapped.sam")
    outputFileDiscarded1 = replaceExtension(getCleanFileName(forward_reads),"_mapped_discarded.fastq")
    outputFileDiscarded2 = replaceExtension(getCleanFileName(reverse_reads),"_mapped_discarded.fastq")
    log_std = "Log.std.out"
    log = "Log.out"
    log_sj = "SJ.out.tab"
    log_final = "Log.final.out"
    log_progress = "Log.progress.out"
    
    if outputFolder is not None and os.path.isdir(outputFolder):
        tmpOutputFile = os.path.join(outputFolder, tmpOutputFile)
        tmpOutputFileDiscarded1 = os.path.join(outputFolder, tmpOutputFileDiscarded1)
        tmpOutputFileDiscarded2 = os.path.join(outputFolder, tmpOutputFileDiscarded2)
        outputFile = os.path.join(outputFolder, outputFile)
        outputFileDiscarded1 = os.path.join( outputFileDiscarded1)
        outputFileDiscarded2 = os.path.join(outputFolder, outputFileDiscarded2)
        log_std = os.path.join(outputFolder, log_std)
        log = os.path.join(outputFolder, log)
        log_sj = os.path.join(outputFolder, log_sj)
        log_final = os.path.join(outputFolder, log_final)
        log_progress = os.path.join(outputFolder, log_progress)
    
    # Options
    #outFilterType(BySJout) this will keep only reads that contains junctions present in SJ.out.tab
    #outSamOrder(Paired) one mate after the other 
    #outSAMprimaryFlag(OneBestScore) only one alignment with the best score is primary
    #outFilterMultimapNmax = read alignments will be output only if the read maps fewer than this value
    #outFilterMismatchNmax = alignment will be output only if it has fewer mismatches than this value
    #outFilterMismatchNoverLmax = alignment will be output only if its ratio of mismatches to *mapped* length is less than this value
    #alignIntronMin minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
    #alignIntronMax maximum intron size, if 0, max intron size will be determined by (2 to the power of winBinNbits)*winAnchorDistNbins
    #alignMatesGapMax maximum gap between two mates, if 0, max intron gap will be determined by (2 to the power of winBinNbits)*winAnchorDistNbins
    #alignEndsType Local standard local alignment with soft-clipping allowed EndToEnd: force end-to-end read alignment, do not soft-clip
    #chimSegmentMin if >0 To switch on detection of chimeric (fusion) alignments
    core_flags = ["--runThreadN", str(max(cores, 1))]
    trim_flags = ["--clip5pNbases", trimForward, trimReverse] 
    io_flags   = ["--outFilterType", "Normal", 
                  "--outSAMtype", "SAM",
                  "--outSAMorder", "Paired",    
                  "--outSAMprimaryFlag", "OneBestScore", 
                  "--outFilterMultimapNmax", 20, 
                  "--alignSJoverhangMin", 8,
                  "--alignSJDBoverhangMin", 1,
                  "--outFilterMismatchNmax", 10, 
                  "--outFilterMismatchNoverLmax", 0.3,
                  "--alignIntronMin", 20,
                  "--alignIntronMax", 0, 
                  "--alignMatesGapMax", 0,
                  "--alignEndsType", "Local", 
                  "--winBinNbits", 16,
                  "--winAnchorDistNbins", 9,
                  "--chimSegmentMin", 0]

    # Main parameters
    #--outSAMtype SAM SortedByCoordinate will output a SAM file sorted by coordinate
    #--quantMode TranscriptomeSAM to output in transcriptome coordinates to use with RSEM
    args = ['STAR']
    args += trim_flags
    args += core_flags
    args += io_flags
    args += ["--genomeDir", ref_map,
             "--readFilesIn", forward_reads, reverse_reads,
             "--outFileNamePrefix", outputFolder + "/"]  # MUST ENSURE AT LEAST ONE SLASH
    args += ["--outReadsUnmapped", "Fastx"]
    
    try:
        proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        error = "Error mapping: STAR execution failed"
        logger.info(error)
        logger.info(e)
        raise
    
    #TODO STAR will output error message if something went wrong
    #should check that here too
    if not fileOk(tmpOutputFile):
        error = "Error mapping: output file is not present : " + tmpOutputFile
        logger.error(error)
        logger.error(stdout)
        logger.error(errmsg)
        raise RuntimeError(error + "\n")
    else:
        # Rename files.
        os.rename(tmpOutputFile, outputFile)
        os.rename(tmpOutputFileDiscarded1, outputFileDiscarded1)
        os.rename(tmpOutputFileDiscarded2, outputFileDiscarded2)
        
        #remove temp files from STAR
        if os.path.isfile(log_std):
            os.remove(log_std)
        if os.path.isfile(log):
            os.remove(log)
        if os.path.isfile(log_sj):
            os.remove(log_sj)
        if os.path.isfile(log_progress):
            os.remove(log_progress)
            
        if not os.path.isfile(log_final):
            logger.info("Warning, log output from STAR is not present")
        else:
            logger.info("Mapping stats: ")
            with open(log_final, "r") as star_log:
                #TODO should only print the useful stats
                for line in star_log.readlines():
                    if str(line) != "":
                        logger.info(str(line))
            os.remove(log_final)
            
    logger.info("Finish STAR Mapping")
    return outputFile, outputFileDiscarded1, outputFileDiscarded2

def barcodeDemultiplexing(readsContainingTr, 
                          idFile, 
                          miss_matches, 
                          kmer, 
                          start_positon,
                          over_hang,
                          outputFolder=None, 
                          keep_discarded_files=False):
    """ 
    :param readsContainingTr a file in SAM/BAM/FASTQ format containing mapped and annotated reads
    :param idFile a tab delimited file (BARCODE - X - Y) containing all the barcodes
    :param miss_matches the number of allowed missmatches 
    :param kmer the kmer length
    :param start_positon the start position of the barcode
    :param over_hang the number of bases to allow for overhang
    :param outputFolder if set output files will be placed there
    :param keep_discarded_files if True files with the non demultiplexed reads will be generated
    This functions performs a demultiplexing using Taggd. Input annotated reads will be filter
    out looking at their barcodes
    """
    
    logger = logging.getLogger("STPipeline")
    logger.info("Start Mapping against the barcodes")
    
    outputFilePrefix = replaceExtension(getCleanFileName(readsContainingTr),'_demultiplexed')
    if outputFolder is not None and os.path.isdir(outputFolder): 
        outputFilePrefix = os.path.join(outputFolder, outputFilePrefix)

    #we know the output file from the prefix and suffix
    #we know the input is SAM so the output will be SAM as well
    outputFile = outputFilePrefix + "_matched.sam"

    # taggd options
    args = ['taggd_demultiplex.py',
            "--max-edit-distance", miss_matches,
            "--k", kmer,
            "--start-position", start_positon,
            "--overhang", over_hang]
    if not keep_discarded_files:
        args.append("--only-output-matched")
    args += [idFile, readsContainingTr, outputFilePrefix]

    try:
        proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdout, errmsg) = proc.communicate()
    except Exception as e:
        error = "Error demultiplexing: taggd execution failed"
        logger.info(error)
        logger.info(e)
        raise
    
    if not fileOk(outputFile):
        error = "Error demultiplexing: output file is not present " + outputFile
        logger.error(error)
        logger.error(stdout)
        logger.error(errmsg)
        raise RuntimeError(error + "\n")
    else:
        procOut = stdout.split("\n")
        logger.info("Barcode Mapping stats :")
        for line in procOut: 
            logger.info(str(line))

    logger.info("Finish Mapping against the barcodes")
    return outputFile