#!/usr/bin/env python
""" 
This module contains wrappers to make systems calls for different annotation tools
most of the options can be passed as arguments
"""

import logging
import subprocess
import os
import pysam
from stpipeline.common.utils import getExtension, Prepender, fileOk


def annotateReadsFeatureCounts(samFile, gtfFile, mode, outputFolder=None):
    #TODO
    # FeatureCounts options
    # -a (annotation file)
    # -o (output file)
    # -F GTF
    # -t exon
    # -g gene_id
    # -f (perform counting at exon level) NO
    # -O (align reads to all their overlapping features)
    # -s 0 (unstranded) 1 (stranded) 2 (reversely stranded) 
    # -M (multimapped reads are counted) NO
    # --largestOverlap (assign reasd to a meta-feature/feature that has the largest number of overlapping genes)
    # --minOverlap 1 (min number of overlapping bases required)
    # --primary (count primary alignments only)
    # --ignoreDup (ignore duplicated reads)
    # -p (count reads pairs instead of individual reads)
    # -P (check validity of paired-end distance)
    # -B (match only read pairs that have both successfully aligned)
    # --donotsort (useful when using pairs)
    # OUTPUT
    # (‘Geneid’, ‘Chr’, ‘Start’, ‘End’, ‘Strand’ and ‘Length’)
    return
     
def annotateReadsHTSeq(samFile, gtfFile, mode, outputFolder=None):
    """ 
    :param samFile sam file contained mapped reads
    :param gtfFile an annotation file in GTF format
    :param mode htseq-count overlapping mode
    :param outputFolder true if we want to place the output file in a given folder
    Annotate the reads using htseq-count tool
    """
    
    logger = logging.getLogger("STPipeline")
    
    sam_type = getExtension(samFile).lower()
    outputFile = 'annotated.' + sam_type
    if outputFolder is not None and os.path.isdir(outputFolder):
        outputFile = os.path.join(outputFolder, outputFile)
    
    # Get the same file header because HTSeq-count will remove it from the output
    samfile = pysam.AlignmentFile(samFile, "r")
    samfile_header = samfile.text
    samfile.close()
    
    #do not want to show the debug messages
    discard_output = open(os.devnull,"w")
    
    #-q (suppress warning reports)
    #-a (min quality)
    #-f (format)
    #-m (annotation mode)
    #-s (strandeness)
    #-i (attribute in GFF to be used as ID)
    #-t (feature type to be used in GFF)
    #-r (input sorted order : name - pos)
    args = ['htseq-count',"-r", "name", "-q", "-a", "0", "-f", sam_type, "-m" , mode, "-s", "no", "-t", 
            "exon", "-i","gene_id" , "-o", outputFile, samFile, gtfFile]
    try:
        subprocess.check_call(args, stdout=discard_output, stderr=subprocess.PIPE)
    except Exception as e:
        error = "Error annotation: HTSEQ execution failed"
        logger.error(error)
        logger.error(e)
        raise
    
    if not fileOk(outputFile):
        error = "Error: output file is not present " + outputFile
        logger.error(error)
        raise RuntimeError(error + "\n")
      
    # Attach back the header to the SAM file
    # Okay, the idea is to attach the SAM header that we extract
    # from the input file before annotation. 
    # Unfortunately, PySAM does not support to prepend a SAM header
    # So we need to parse the header and prepend the lines to the SAM
    # output from htseq-count.
    with Prepender(outputFile) as f:
        f.write_lines(samfile_header.split("\n"))

    return outputFile