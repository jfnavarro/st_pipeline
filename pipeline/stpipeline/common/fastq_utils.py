#!/usr/bin/env python
""" 
This module contains some functions to deal with fastq files
"""

from stpipeline.common.utils import *
from stpipeline.common.adaptors import removeAdaptor
import logging 
from itertools import izip

def coroutine(func):
    """ 
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

def quality_trim_index(qualities, cutoff, base=33):
    """
    NOTE : function snippet from CutAdapt 
    https://code.google.com/p/cutadapt/
    
    Find the position at which to trim a low-quality end from a nucleotide sequence.

    Qualities are assumed to be ASCII-encoded as chr(qual + base).

    The algorithm is the same as the one used by BWA within the function
    'bwa_trim_read':
    - Subtract the cutoff value from all qualities.
    - Compute partial sums from all indices to the end of the sequence.
    - Trim sequence at the index at which the sum is minimal.
    """
    s = 0
    max_qual = 0
    max_i = len(qualities)
    for i in reversed(xrange(max_i)):
        q = ord(qualities[i]) - base
        s += cutoff - q
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            max_i = i
    return max_i

def trim_quality(record, trim_distance, min_qual=20, 
                 min_length=28, qual64=False):    
    """
    :param record the fastq read (name,sequence,quality)
    :param trim_distance the number of bases to be trimmed (not considered)
    :param min_qual the quality threshold to trim
    :param min_length the min length of a valid read after trimming
    :param qual64 true of the qualities are in phred64 format
    Perfoms a bwa-like quality trimming on the sequence and 
    quality in tuple record(name,seq,qual)
    Returns the trimmed read or None if the read has to be discarded
    """
    qscore = record[2][trim_distance:]
    sequence = record[1][trim_distance:]
    name = record[0]
    phred = 64 if qual64 else 33
    
    #get the position at which to trim
    cut_index = quality_trim_index(qscore, min_qual, phred)
    #get the number of bases suggested to trim
    nbases = len(qscore) - cut_index
    
    #check if the trimmed sequence would have min length (at least)
    #if so return the trimmed read otherwise return None
    if (len(sequence) - nbases) >= min_length:
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
    return record[0], new_seq, new_qual

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
  
def reverse_complement(seq):
    """
    :param seq a FASTQ sequence
    This functions returns the reverse complement
    of the sequence given as input
    """
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}   
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases
  
def reformatRawReads(fw, rw, 
                     barcode_start=0, barcode_length=18,
                     filter_AT_content=90,
                     molecular_barcodes=False, mc_start=18, mc_end=27,
                     trim_fw=42, trim_rw=5,
                     min_qual=20, min_length=28,
                     polyA_min_distance=0, polyT_min_distance=0, 
                     polyG_min_distance=0, polyC_min_distance=0,
                     qual64=False, outputFolder=None, keep_discarded_files=False):
    """ 
    :param fw the fastq file with the forward reads
    :param rw the fastq file with the reverse reads
    :param barcode_start the base where the barcode sequence starts
    :param barcode_length the number of bases of the barcodes
    :param molecular_barcodes if True the forward reads contain molecular barcodes
    :param mc_start the start position of the molecular barcodes if any
    :param mc_end the end position of the molecular barcodes if any
    :param trim_fw how many bases we want to trim (not consider in the forward)
    :param trim_rw how many bases we want to trim (not consider in the reverse)
    :param min_qual the min quality value to use to trim quality
    :param min_length the min valid length for a read after trimming
    :param polyA_min_distance if >0 we remove PolyA adaptors from the reads
    :param polyT_min_distance if >0 we remove PolyT adaptors from the reads
    :param polyG_min_distance if >0 we remove PolyG adaptors from the reads
    :param qual64 true of qualities are in phred64 format
    :param outputFolder optional folder to output files
    :param keep_discarded_files when true files containing the discarded reads will be created
    This function does three things (all here for speed optimization)
      - It appends the barcode and the molecular barcode (if any)
        from forward reads to reverse reads
      - It performs a BWA quality trimming discarding very short reads
      - It removes adaptors from the reads (optional)
    """
    logger = logging.getLogger("STPipeline")
    logger.info("Start Reformatting and Filtering raw reads")
    
    out_rw = 'R1_trimmed_formated.fastq'
    out_fw = 'R2_trimmed_formated.fastq'
    out_fw_discarded = 'R1_trimmed_formated_discarded.fastq'
    out_rw_discarded = 'R2_trimmed_formated_discarded.fastq'
    
    if outputFolder is not None and os.path.isdir(outputFolder):
        out_rw = os.path.join(outputFolder, out_rw)
        out_fw = os.path.join(outputFolder, out_fw)
        out_fw_discarded = os.path.join(outputFolder, out_fw_discarded)
        out_rw_discarded = os.path.join(outputFolder, out_rw_discarded)
    
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
    
    adaptorA = "".join("A" for k in xrange(polyA_min_distance))
    adaptorT = "".join("T" for k in xrange(polyT_min_distance))
    adaptorG = "".join("G" for k in xrange(polyG_min_distance))
    adaptorC = "".join("C" for k in xrange(polyG_min_distance))
    
    for line1, line2 in izip(readfq(fw_file), readfq(rw_file)):

        if line1[0].split()[0] != line2[0].split()[0]:
            logger.warning("Pair raids found with different names " + line1[0] + " and " + line2[0])
        
        total_reads += 1
        
        #get the barcode and molecular barcode if any from the forward read
        #to be attached to the reverse read
        to_append_sequence = line1[1][barcode_start:barcode_length]
        to_append_sequence_quality = line1[2][barcode_start:barcode_length]
        if molecular_barcodes:
            to_append_sequence += line1[1][mc_start:mc_end]
            to_append_sequence_quality += line1[2][mc_start:mc_end]
        
        original_line1 = line1
        original_line2 = line2
        
        # If read - trimming is not long enough or has a high AT content discard...
        if (len(line1[1]) - trim_fw) < min_length or \
        ((line1[1].count("A") + line1[1].count("T")) / len(line1[1])) * 100 >= filter_AT_content:
            line1 = None
        if (len(line2[1]) - trim_rw) < min_length or \
        ((line2[1].count("A") + line2[1].count("T")) / len(line2[1])) * 100 >= filter_AT_content:
            line2 = None
              
        # if indicated we remove the adaptor PolyA from both reads
        if polyA_min_distance > 0:
            line1 = removeAdaptor(line1, adaptorA, trim_fw, "5")
            line2 = removeAdaptor(line2, adaptorA, trim_rw, "5")
            
        # if indicated we remove the adaptor PolyT from both reads
        if polyT_min_distance > 0:
            line1 = removeAdaptor(line1, adaptorT, trim_fw, "5")
            line2 = removeAdaptor(line2, adaptorT, trim_rw, "5")
       
        # if indicated we remove the adaptor PolyG from both reads
        if polyG_min_distance > 0:
            line1 = removeAdaptor(line1, adaptorG, trim_fw, "5")
            line2 = removeAdaptor(line2, adaptorG, trim_rw, "5")
        
        # if indicated we remove the adaptor PolyC from both reads
        if polyC_min_distance > 0:
            line1 = removeAdaptor(line1, adaptorC, trim_fw, "5")
            line2 = removeAdaptor(line2, adaptorC, trim_rw, "5")
          
        line2_trimmed = None
        line1_trimmed = None
             
        # Trim rw
        if line2 is not None:
            line2_trimmed = trim_quality(line2, trim_rw, min_qual, min_length, qual64)
            
        # Trim fw
        if line1 is not None:
            line1_trimmed = trim_quality(line1, trim_fw, min_qual, min_length, qual64)
        
        if line1_trimmed is not None:
            out_fw_writer.send(line1_trimmed)
        else:
            # write fake sequence so mapping wont fail for having rv and fw with different lengths
            out_fw_writer.send(getFake(original_line1))
            dropped_fw += 1
            if keep_discarded_files:
                out_fw_writer_discarded.send(original_line1)

        if line2_trimmed is not None:
            # Attach the barcode from fw only if rw has not been completely trimmed
            new_seq = to_append_sequence + line2_trimmed[1]
            new_qual = to_append_sequence_quality + line2_trimmed[2]
            out_rw_writer.send((line2_trimmed[0], new_seq, new_qual))
        else:
            # write fake sequence so mapping wont fail for having rv and fw with different lengths
            out_rw_writer.send(getFake(original_line2))
            dropped_rw += 1  
            if keep_discarded_files:
                out_rw_writer_discarded.send(original_line2)
    
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
        logger.info("Trimming stats total reads : " + str(total_reads * 2))
        logger.info("Trimming stats fw 1 : " + str(dropped_fw) + " reads have been dropped on the forward reads!")
        perc1 = '{percent:.2%}'.format(percent= float(dropped_fw) / float(total_reads) )
        logger.info("Trimming stats fw 2 : you just lost about " + perc1 + " of your data on the forward reads!")
        logger.info("Trimming stats rv 1 : " + str(dropped_rw) + " reads have been dropped on the reverse reads!") 
        perc2 = '{percent:.2%}'.format(percent= float(dropped_rw) / float(total_reads) )
        logger.info("Trimming stats rv 2 : you just lost about " + perc2 + " of your data on the reverse reads!")
        logger.info("Trimming stats reads remaining: " + str((total_reads * 2) - dropped_fw - dropped_rw))
        
    logger.info("Finish Reformatting and Filtering raw reads")
    return out_fw, out_rw


