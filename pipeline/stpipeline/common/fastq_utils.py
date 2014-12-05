#!/usr/bin/env python
""" 
This module contains some functions to deal with fastq files
"""

from stpipeline.common.utils import *
import logging 
from itertools import izip
from cutadapt.qualtrim import quality_trim_index

def coroutine(func):
    """ 
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

def trim_quality(record, trim_distance, min_qual=20, 
                 min_length=28, qual64=False):    
    """
    Perfoms a bwa-like quality trimming on the sequence and 
    quality in tuple record(name,seq,qual)
    """
    qscore = record[2]
    sequence = record[1]
    name = record[0]
    phred = 64 if qual64 else 33
    
    #get the position at which to trim
    cut_index = quality_trim_index(qscore, min_qual, phred)
    #get the number of bases suggested to trim
    nbases = len(qscore) - cut_index
    
    #check if the trimmed sequence would have min length accounting for trim_distance
    #offset 
    if (len(sequence) - (trim_distance + nbases)) >= min_length:
        new_seq = record[1][:(len(sequence) - nbases)]
        new_qual = record[2][:(len(sequence) - nbases)]
        return name, new_seq, new_qual
    else:
        return None
    
def getFake(record):
    """ 
    Generates a fake fastq record(name,seq,qual) from the record given as input
    """
    new_seq = "".join("N" for k in record[1])
    new_qual = "".join("B" for k in record[2])
    return (record[0],new_seq,new_qual)

def readfq(fp): # this is a generator function
    """ 
    Heng Li's fasta/fastq reader function.
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

@coroutine
def writefq(fp):  # This is a coroutine
    """ 
    Fastq writing generator sink.
    Send a (header, sequence, quality) triple to the instance to write it to
    the specified file pointer.
    """
    fq_format = '@{header}\n{sequence}\n+\n{quality}\n'
    try:
        while True:
            record = yield
            read = fq_format.format(header=record[0], sequence=record[1], quality=record[2])
            fp.write(read)

    except GeneratorExit:
        return

def reformatRawReads(fw, rw, trim_fw=42, trim_rw=5,
                     min_qual=20, min_length=28, qual64=False, outputFolder=None, keep_discarded_files=False):
    """ 
    Converts reads in rw file appending the first (distance - trim)
    bases of fw and also add FW or RW string to reads names
    It also performs a bwa quality trim of the fw and rw reads, when
    the trimmed read is below min lenght it will discarded.
    """
    logger = logging.getLogger("STPipeline")
    
    if fw.endswith(".fastq") and rw.endswith(".fastq"):
        out_rw = replaceExtension(getCleanFileName(rw),'_formated.fastq')
        if outputFolder is not None and os.path.isdir(outputFolder) : 
            out_rw = os.path.join(outputFolder, out_rw)
            
        out_fw = replaceExtension(getCleanFileName(fw),'_formated.fastq')
        if outputFolder is not None and os.path.isdir(outputFolder): 
            out_fw = os.path.join(outputFolder, out_fw)
        
        out_fw_discarded = replaceExtension(getCleanFileName(fw),'_formated_discarded.fastq')
        if outputFolder is not None and os.path.isdir(outputFolder): 
            out_fw_discarded = os.path.join(outputFolder, out_fw_discarded)
            
        out_rw_discarded = replaceExtension(getCleanFileName(rw),'_formated_discarded.fastq')
        if outputFolder is not None and os.path.isdir(outputFolder) : 
            out_rw_discarded = os.path.join(outputFolder, out_rw_discarded)
            
    else:
        error = "Error: Input format not recognized " + out_fw + " , " + out_rw
        logger.error(error)
        raise RuntimeError(error + "\n")

    logger.info("Start Reformatting and Filtering raw reads")
    
    out_fw_handle = safeOpenFile(out_fw, 'w')
    out_fw_writer = writefq(out_fw_handle)
    
    out_rw_handle = safeOpenFile(out_rw, 'w')
    out_rw_writer = writefq(out_rw_handle)

    if keep_discarded_files:
        out_rw_handle_discarded = safeOpenFile(out_rw_discarded, 'w')
        out_rw_writer_discarded = writefq(out_rw_handle_discarded)
        
        out_fw_handle_discarded = safeOpenFile(out_fw_discarded, 'w')
        out_fw_writer_discarded = writefq(out_fw_handle_discarded)
    
        
    fw_file = safeOpenFile(fw, "rU")
    rw_file = safeOpenFile(rw, "rU")

    total_reads = 0
    dropped_fw = 0
    dropped_rw = 0

    for line1, line2 in izip(readfq(fw_file), readfq(rw_file)):
        total_reads += 1

        # Trim rw
        line2_trimmed = trim_quality(line2, trim_rw, min_qual, min_length, qual64)
        # Trim fw
        line1_trimmed = trim_quality(line1, trim_fw, min_qual, min_length, qual64)

        if line1_trimmed is not None:
            out_fw_writer.send(line1_trimmed)
        else:
            # write fake sequence so bowtie wont fail for having rw and fw with different lenghts
            out_fw_writer.send(getFake(line1))
            dropped_fw += 1
            if keep_discarded_files:
                out_fw_writer_discarded.send(line1)

        if line2_trimmed is not None:
            # Add the barcode and polyTs from fw only if rw has not been completely trimmed
            # we do not include in the new seq/qual the number of reverse bases that we want
            # to trim so they are not included in the mapping
            new_seq = line1[1][:trim_fw] + line2_trimmed[1][trim_rw:]
            new_qual = line1[2][:trim_fw] + line2_trimmed[2][trim_rw:]
            record = (line2_trimmed[0], new_seq, new_qual)
            out_rw_writer.send(record)

        else:
            # write fake sequence so bowtie wont fail for having rw and fw with different lenghts
            out_rw_writer.send(getFake(line2))
            dropped_rw += 1  
            if keep_discarded_files:
                out_rw_writer_discarded.send(line2)
    
    out_fw_writer.close()
    out_rw_writer.close()
    out_fw_handle.close()
    out_rw_handle.close()
    if keep_discarded_files:
        out_fw_handle_discarded.close()
        out_rw_handle_discarded.close()
    fw_file.close()
    rw_file.close()
    
    if not fileOk(out_fw) or not fileOk(out_rw):
        error = "Error: output file is not present " + out_fw + " , " + out_rw
        logger.error(error)
        raise RuntimeError(error + "\n")
    else:
        logger.info("Trimming stats total reads : " + str(total_reads))
        logger.info("Trimming stats fw 1 : " + str(dropped_fw) + " reads have been dropped on the forward reads!")
        perc1 = '{percent:.2%}'.format(percent= float(dropped_fw) / float(total_reads) )
        logger.info("Trimming stats fw 2 : you just lost about " + perc1 + " of your data on the forward reads!")
        logger.info("Trimming stats rw 1 : " + str(dropped_rw) + " reads have been dropped on the reverse reads!") 
        perc2 = '{percent:.2%}'.format(percent= float(dropped_rw) / float(total_reads) )
        logger.info("Trimming stats rw 2 : you just lost about " + perc2 + " of your data on the reverse reads!")
        logger.info("Trimming stats reads remaining: " + str(total_reads - dropped_fw - dropped_rw))
        
    logger.info("Finish Reformatting and Filtering raw reads")
    
    return out_fw, out_rw


