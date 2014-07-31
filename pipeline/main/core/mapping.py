#!/usr/bin/env python
"""
    Copyright (C) 2012  Spatial Transcriptomics AB,
    read LICENSE for licensing terms. 
    Contact : Jose Fernandez Navarro <jose.fernandez.navarro@scilifelab.se>

"""

""" This class contains wrappers to make systems calls for different aligners
most of the options can be passed as arguments
"""

import logging 
import os
import sys
import subprocess
from main.common.utils import *
from main.common.fastq_utils import *
import pysam
    
def bowtie2Map(fw, rv, ref_map, trim = 42, cores = 8, qual64 = False, discordant = False):  
    ''' maps pair end reads against a given genome using bowtie2 
    '''
    
    logger = logging.getLogger("STPipeline")
    
    if fw.endswith(".fastq") and rv.endswith(".fastq"):
        outputFile = replaceExtension(fw,".sam")
    else:
        logger.error("Error: Input format not recognized " + fw + " , " + rv)
        raise RuntimeError("Error: Input format not recognized")
    
    logger.info("Start Bowtie2 Mapping, output = " + outputFile)
    
    qual_flags = ["--phred64"] if qual64 else ["--phred33"] 
    core_flags = ["-p", str(cores)] if cores > 1 else []
    trim_flags = ["--trim5", trim] 
    io_flags   = ["-q", "-X", 2000, "--sensitive"] ##500 (default) is too selective
    io_flags  += ["--no-discordant"] if discordant else []

    args = ['bowtie2']
    args += trim_flags
    args += qual_flags
    args += core_flags
    args += io_flags
    args += ["-x", ref_map, "-1", fw, "-2", rv, "-S", outputFile] 
    
    proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, errmsg) = proc.communicate()

    if not fileOk(outputFile):
        logger.error("Error: output file is not present : " + outputFile)
        logger.error(stdout)
        logger.error(errmsg)
        raise RuntimeError("Error: output file is not present")
    else:
        procOut = [x for x in errmsg.split("\n") if x.find("Warning") == -1 and x.find("Error") == -1]
        logger.info('Mapping stats on paired end mode with 5-end trimming of ' + str(trim))
        for line in procOut:
            logger.info(str(line))

    logger.info("Finish Bowtie2 Mapping")

    return outputFile


def bowtie2_contamination_map(fastq, contaminant_index, trim=42, cores=8, qual64=False):
    """ Maps reads against contaminant genome index with Bowtie2 and returns
    the fastq of unaligned reads.
    """
    
    logger = logging.getLogger("STPipeline")

    if fastq.endswith(".fastq"):
        contaminated_file = replaceExtension(fastq,"_contaminated.sam")
        clean_fastq = replaceExtension(fastq,"_clean.fastq")
    else:
        logger.error("Error: Input format not recognized " + fastq)
        raise RuntimeError("Error: Input format not recognized")

    logger.info("Start Bowtie2 Mapping rRNA")
    
    qual_flags = ["--phred64"] if qual64 else ["--phred33"]
    core_flags = ["-p", str(cores)] if cores > 1 else []
    trim_flags = ["--trim5", trim]
    io_flags   = ["-q","--sensitive"]
    
    args = ['bowtie2']
    args += trim_flags
    args += qual_flags
    args += core_flags
    args += io_flags
    args += ["-x", contaminant_index, "-U", fastq, "-S", contaminated_file, "--un", clean_fastq]

    proc = subprocess.Popen([str(i) for i in args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, errmsg) = proc.communicate()

    if not fileOk(contaminated_file) or not fileOk(clean_fastq):
        logger.error("Error: output file is not present " + contaminated_file + " , " + clean_fastq)
        logger.error(stdout)
        logger.error(errmsg)
        raise RuntimeError("Error: output file is not present")
    else:
        procOut = [x for x in errmsg.split("\n") if x.find("Warning") == -1 and x.find("Error") == -1]
        logger.info('Contaminant mapping stats against ' +
                       contaminant_index + ' with 5-end trimming of ' +
                       str(trim))
        for line in procOut:
            logger.info(str(line))

    logger.info("Finish Bowtie2 Mapping rRNA")

    return clean_fastq, contaminated_file

def filterUnmapped(sam, discard_fw=False, discard_rw=False):
    ''' filter unmapped and discordant reads
    '''
    
    logger = logging.getLogger("STPipeline")

    if sam.endswith(".sam"):
        outputFileSam = replaceExtension(sam,'_filtered.sam')
    else:
        logger.error("Error: Input format not recognized " + sam)
        raise RuntimeError("Error: Input format not recognized")

    logger.info("Start converting Sam to Bam and filtering unmapped")
    
    # Remove found duplicates in the Forward Reads File
    input = pysam.Samfile(sam, "r")
    output = pysam.Samfile(outputFileSam, 'wh', header=input.header)

    for read in input:
        # filtering out not pair end reads
        if not read.is_paired:
            logger.error("Error: input Sam file contains not paired reads")
            raise RuntimeError("Error: input Sam file contains not paired reads")

        if read.is_proper_pair and not read.mate_is_unmapped:
            # if read is a concordant pair and it is mapped,
            output.write(read)
        elif not read.is_proper_pair and not read.is_unmapped:
            # if read is a discordant mapped pair or uniquely mapped (mixed mode bowtie2)
            if read.is_read1 and not discard_fw:
                # if read is forward and we dont want to discard it
                output.write(read)
            elif read.is_read2 and not discard_rw:
                # if read is reverse and we dont want to discard it
                output.write(read)
            else:
                # I want to discard both forward and reverse and unmapped
                pass
        else:
            # not mapped stuff discard
            pass  
            
    input.close()
    output.close()

    if not fileOk(outputFileSam):
        logger.error("Error: output file is not present " + outputFileSam)
        raise RuntimeError("Error: output file is not present")
    
    logger.info("End converting Sam to Bam and filtering unmapped")
    
    return outputFileSam

def getTrToIdMap(readsContainingTr, idFile, m, k, s, l, e):
    ''' barcode demultiplexing mapping with old findindexes 
    '''
    
    logger = logging.getLogger("STPipeline")
    
    if not fileOk(readsContainingTr) or not fileOk(idFile):
        logger.error("Error: Input files not present , transcript file = " + readsContainingTr + " ids = " + idFile)
        raise RuntimeError("Error: Input format not recognized")
    
    logger.info("Start Mapping against the barcodes")
    
    outputFile = replaceExtension(readsContainingTr,'_nameMap.txt')

    args = ['findIndexes',
            "-m", str(m), "-k", str(k), "-s", str(s),
            "-l", str(l), "-o", str(outputFile), idFile,
            readsContainingTr]
    
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, errmsg) = proc.communicate()
 
    if not os.path.isfile(outputFile) or os.path.getsize(outputFile) == 0:
        logger.error("Error: output file is not present " + outputFile)
        logger.error(stdout)
        logger.error(errmsg)
        raise Exception("Error: output file is not present")
    else:
        procOut = stdout.split("\n")
        logger.info('Barcode Mapping stats :')
        for line in procOut: 
            logger.info(str(line))

    logger.info("Finish Mapping against the barcodes")
    
    return outputFile