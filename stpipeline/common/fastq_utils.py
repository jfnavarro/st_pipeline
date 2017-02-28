""" 
This module contains some specific functionalities for
ST fastq files, mainly quality filtering functions.
"""

from stpipeline.common.utils import safeOpenFile, fileOk
from stpipeline.common.adaptors import removeAdaptor
from stpipeline.common.stats import qa_stats
import logging 
from itertools import izip
from sqlitedict import SqliteDict
import os
import re

def coroutine(func):
    """ 
    Coroutine decorator, starts coroutines upon initialization.
    """
    def start(*args, **kwargs):
        cr = func(*args, **kwargs)
        cr.next()
        return cr
    return start

def readfq(fp): # this is a generator function
    """ 
    Heng Li's fasta/fastq reader function.
    # https://github.com/lh3/readfq/blob/master/readfq.py
    # Unlicensed. 
    Parses fastq records from a file using a generator approach.
    :param fp: opened file descriptor
    :returns an iterator over tuples (name,sequence,quality)
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
    
def quality_trim_index(bases, qualities, cutoff, base=33):
    """
    Function snippet and modified from CutAdapt 
    https://github.com/marcelm/cutadapt/
    
    Copyright (c) 2010-2016 Marcel Martin <marcel.martin@scilifelab.se>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN C

    Find the position at which to trim a low-quality end from a nucleotide sequence.

    Qualities are assumed to be ASCII-encoded as chr(qual + base).

    The algorithm is the same as the one used by BWA within the function
    'bwa_trim_read':
    - Subtract the cutoff value from all qualities.
    - Compute partial sums from all indices to the end of the sequence.
    - Trim sequence at the index at which the sum is minimal.
    
    This variant works on NextSeq data.
    With Illumina NextSeq, bases are encoded with two colors. 'No color' (a
    dark cycle) usually means that a 'G' was sequenced, but that also occurs
    when sequencing falls off the end of the fragment. The read then contains
    a run of high-quality G bases in the end.
    This routine works as the one above, but counts qualities belonging to 'G'
    bases as being equal to cutoff - 1.
    """
    s = 0
    max_qual = 0
    max_i = len(qualities)
    for i in reversed(xrange(max_i)):
        q = ord(qualities[i]) - base
        if bases[i] == 'G':
            q = cutoff - 1
        s += cutoff - q
        if s < 0:
            break
        if s > max_qual:
            max_qual = s
            max_i = i
    return max_i

def trim_quality(sequence,
                 quality,
                 min_qual=20, 
                 min_length=30, 
                 phred=33):    
    """
    Quality trims a fastq read using a BWA approach.
    It returns the trimmed record or None if the number of bases
    after trimming is below a minimum.
    :param sequence: the sequence of bases of the read
    :param quality: the quality scores of the read
    :param min_qual the quality threshold to trim (consider a base of bad quality)
    :param min_length: the minimum length of a valid read after trimming
    :param phred: the format of the quality string (33 or 64)
    :type sequence: str
    :type quality: str
    :type min_qual: integer
    :type min_length: integer
    :type phred: integer
    :return: A tuple (base, qualities) or (None,None)
    """
    if len(sequence) < min_length:
        return None, None
    # Get the position at which to trim (number of bases to trim)
    cut_index = quality_trim_index(sequence, quality, min_qual, phred)
    # Check if the trimmed sequence would have min length (at least)
    # if so return the trimmed read otherwise return None
    if (cut_index + 1) >= min_length:
        new_seq = sequence[:cut_index]
        new_qual = quality[:cut_index]
        return new_seq, new_qual
    else:
        return None, None
  
def check_umi_template(umi, template):
    """
    Checks that the UMI (molecular barcode) given as input complies
    with the pattern given in template.
    Returns True if the UMI complies
    :param umi: a molecular barcode
    :param template: a reg-based template with the same
                    distance of the UMI that should tell how the UMI should be formed
    :type umi: str
    :type template: str
    :return: True if the given molecular barcode fits the pattern given
    """
    p = re.compile(template)
    return p.match(umi) is not None

def filterInputReads(fw, 
                     rv,
                     out_fw,
                     out_rw,
                     out_rw_discarded,
                     filter_AT_content,
                     filter_GC_content,
                     umi_start, 
                     umi_end,
                     min_qual, 
                     min_length,
                     polyA_min_distance, 
                     polyT_min_distance, 
                     polyG_min_distance, 
                     polyC_min_distance,
                     qual64,
                     umi_filter,
                     umi_filter_template,
                     umi_quality_bases):
    """
    This function does few things (all done in one loop for performance reasons)
      - It performs a sanity check (forward and reverse reads same length and order)
      - It performs a BWA-based quality trimming discarding very short reads
      - It removes adaptors from the reads (optional)
      - It checks for AT and GC content (optional)
      - It performs a sanity check on the UMI (optional)
    Reads that do not pass the filters are discarded (both R1 and R2)
    :param fw: the fastq file with the forward reads
    :param rv: the fastq file with the reverse reads
    :param out_fw: the name of the output file for the forward reads
    :param out_rw: the name of the output file for the reverse reads
    :param out_rw_discarded: the name of the output file for discarded reverse reads
    :param filter_AT_content: % of A and T bases a read2 must have to be discarded
    :param filter_GC_content: % of G and C bases a read2 must have to be discarded
    :param umi_start: the start position of the UMI
    :param umi_end: the end position of the UMI
    :param min_qual: the min quality value to use in the trimming
    :param min_length: the min valid length for a read after trimming
    :param polyA_min_distance: if >10 we remove PolyA adaptors from the reads
    :param polyT_min_distance: if >10 we remove PolyT adaptors from the reads
    :param polyG_min_distance: if >10 we remove PolyG adaptors from the reads
    :param qual64: true of qualities are in phred64 format
    :param umi_filter: performs a UMI quality template filter when True
    :param umi_filter_template: the template to use for the UMI filter
    :param umi_quality_bases: the number of low quality bases allowed in an UMI
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(fw) or not os.path.isfile(rv):
        error = "Error doing quality trimming, input file/s not present {}\n{}\n".format(fw,rv)
        logger.error(error)
        raise RuntimeError(error)
    
    # Check if discarded files must be written out 
    keep_discarded_files = out_rw_discarded is not None
    
    # Create output file writers
    out_rv_handle = safeOpenFile(out_rw, 'w')
    out_rv_writer = writefq(out_rv_handle)
    out_fw_handle = safeOpenFile(out_fw, 'w')
    out_fw_writer = writefq(out_fw_handle)
    if keep_discarded_files:
        out_rv_handle_discarded = safeOpenFile(out_rw_discarded, 'w')
        out_rv_writer_discarded = writefq(out_rv_handle_discarded)
    
    # Some counters
    total_reads = 0
    dropped_rv = 0
    dropped_umi = 0
    dropped_umi_template = 0
    dropped_AT = 0
    dropped_GC = 0
    dropped_adaptor = 0
    
    # Build fake sequence adaptors with the parameters given
    adaptorA = "".join("A" for k in xrange(polyA_min_distance))
    adaptorT = "".join("T" for k in xrange(polyT_min_distance))
    adaptorG = "".join("G" for k in xrange(polyG_min_distance))
    adaptorC = "".join("C" for k in xrange(polyC_min_distance))
    # Not recommended to do adaptor trimming for adaptors smaller than 8
    do_adaptorA = polyA_min_distance >= 8
    do_adaptorT = polyT_min_distance >= 8
    do_adaptorG = polyG_min_distance >= 8
    do_adaptorC = polyC_min_distance >= 8
    # Allow for 2 miss-matches when removing adaptors
    adaptor_missmatches = 2
    
    # Quality format
    phred = 64 if qual64 else 33
    
    # Open fastq files with the fastq parser
    fw_file = safeOpenFile(fw, "rU")
    rv_file = safeOpenFile(rv, "rU")
    for (header_fw, sequence_fw, quality_fw), (header_rv, sequence_rv, quality_rv) \
    in izip(readfq(fw_file), readfq(rv_file)):
        
        if not sequence_fw or not sequence_rv:
            error = "Error doing quality trimming, Checks of raw reads.\n" \
            "The input files {},{} are not of the same length".format(fw,rv)
            logger.error(error)
            break
        
        if header_fw.split()[0] != header_rv.split()[0]:
            logger.warning("Pair reads found with different " \
                           "names {} and {}".format(header_fw,header_rv))
            
        # Increase reads counter
        total_reads += 1
        discard_read = False
        
        # If we want to check for UMI quality and the UMI is incorrect
        # then we discard the reads
        if umi_filter \
        and not check_umi_template(sequence_fw[umi_start:umi_end], umi_filter_template):
            dropped_umi_template += 1
            discard_read = True
        
        # Check if the UMI has many low quality bases
        if not discard_read and \
        len([b for b in quality_fw[umi_start:umi_end] if (ord(b) - phred) < min_qual]) > umi_quality_bases:
            dropped_umi += 1
            discard_read = True

        # If reverse read has a high AT content discard...
        if not discard_read and \
        ((sequence_rv.count("A") + sequence_rv.count("T")) / len(sequence_rv)) * 100 >= filter_AT_content:
            dropped_AT += 1
            discard_read = True

        # If reverse read has a high GC content discard...
        if not discard_read and \
        ((sequence_rv.count("G") + sequence_rv.count("C")) / len(sequence_rv)) * 100 >= filter_GC_content:
            dropped_GC += 1
            discard_read = True
                           
        # Store the original reads to write them to the discarded output if applies
        if keep_discarded_files:    
            orig_sequence_rv = sequence_rv
            orig_quality_rv = quality_rv 
            
        if not discard_read:
            # if indicated we remove the artifacts PolyA from reverse reads
            if do_adaptorA: 
                sequence_rv, quality_rv = removeAdaptor(sequence_rv, quality_rv, 
                                                        adaptorA, adaptor_missmatches) 
            # if indicated we remove the artifacts PolyT from reverse reads
            if do_adaptorT: 
                sequence_rv, quality_rv = removeAdaptor(sequence_rv, quality_rv, 
                                                        adaptorT, adaptor_missmatches) 
            # if indicated we remove the artifacts PolyG from reverse reads
            if do_adaptorG: 
                sequence_rv, quality_rv = removeAdaptor(sequence_rv, quality_rv, 
                                                        adaptorG, adaptor_missmatches) 
            # if indicated we remove the artifacts PolyC from reverse reads
            if do_adaptorC: 
                sequence_rv, quality_rv = removeAdaptor(sequence_rv, quality_rv, 
                                                        adaptorC, adaptor_missmatches)
            # Check if the read is smaller than the minimum after removing artifacts   
            if len(sequence_rv) < min_length:
                dropped_adaptor += 1
                discard_read = True
            else:              
                # Trim reverse read (will return None if length of trimmed sequence is less than min_length)
                sequence_rv, quality_rv = trim_quality(sequence_rv, quality_rv, 
                                                       min_qual, min_length, phred)
                if not sequence_rv or not quality_rv:
                    discard_read = True
            
        # Write reverse read to output
        if not discard_read:
            out_rv_writer.send((header_rv, sequence_rv, quality_rv))
            out_fw_writer.send((header_fw, sequence_fw, quality_fw))
        else:
            dropped_rv += 1  
            if keep_discarded_files:
                out_rv_writer_discarded.send((header_rv, orig_sequence_rv, orig_quality_rv))
    
    fw_file.close()
    rv_file.close()
    out_rv_handle.flush()
    out_rv_handle.close()
    out_rv_writer.close()
    out_fw_handle.flush()
    out_rv_writer.close()
    out_fw_writer.close()
    if keep_discarded_files:
        out_rv_handle_discarded.flush()
        out_rv_handle_discarded.close()
        out_rv_writer_discarded.close()
        
    # Write info to the log
    logger.info("Trimming stats total reads (pair): {}".format(total_reads))
    logger.info("Trimming stats {} reads have been dropped!".format(dropped_rv)) 
    perc2 = '{percent:.2%}'.format(percent= float(dropped_rv) / float(total_reads) )
    logger.info("Trimming stats you just lost about {} of your data".format(perc2))
    logger.info("Trimming stats reads remaining: {}".format(total_reads - dropped_rv))
    logger.info("Trimming stats dropped pairs due to incorrect UMI: {}".format(dropped_umi_template))
    logger.info("Trimming stats dropped pairs due to low quality UMI: {}".format(dropped_umi))
    logger.info("Trimming stats dropped pairs due to high AT content: {}".format(dropped_AT))
    logger.info("Trimming stats dropped pairs due to high GC content: {}".format(dropped_GC))
    logger.info("Trimming stats dropped pairs due to presence of artifacts: {}".format(dropped_adaptor))
    
    # Check that output file was written ok
    if not fileOk(out_rw):
        error = "Error doing quality trimming checks of raw reads." \
        "\nOutput file not present {}\n".format(out_rw)
        logger.error(error)
        raise RuntimeError(error)
    
    # Adding stats to QA Stats object
    qa_stats.input_reads_forward = total_reads
    qa_stats.input_reads_reverse = total_reads
    qa_stats.reads_after_trimming_forward = total_reads
    qa_stats.reads_after_trimming_reverse = (total_reads - dropped_rv)

#TODO this approach uses too much memory
#     find a better solution (maybe Cython or C++)
def hashDemultiplexedReads(reads,
                           umi_start,
                           umi_end,
                           low_memory):
    """
    This function extracts the read name, the UMI and the x,y coordinates
    from the reads given as input (output of TaggD) and returns a dictionary
    with the clean read name as key and (x,y,umi) as value. X and Y correspond
    to the array coordinates of the barcode of the read and UMI is extracted from the read
    sequence.
    :param reads: path to a file with the fastq reads after demultiplexing (TaggD)
    :param umi_start: the start position of the UMI
    :param umi_end: the end position of the UMI
    :param low_memory: True to use a key-value db instead of dict
    :type reads: str
    :type umi_start: integer
    :type umi_end: integer
    :type low_memory: boolean
    :return: a dictionary of read_name -> (x,y,umi) tags where umi is optional
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(reads):
        error = "Error parsing TaggD output, input file not present {}\n".format(reads)
        logger.error(error)
        raise RuntimeError(error)
    
    if low_memory:
        hash_reads = SqliteDict(autocommit=False, flag='c', journal_mode='OFF')
    else:
        hash_reads = dict()
    
    fastq_file = safeOpenFile(reads, "rU")
    for name, sequence, _ in readfq(fastq_file):
        # Assumes the header is like this
        # @NS500688:111:HNYW7BGXX:1:11101:13291:1099 1:N:0:TGCCCA B0:Z:GTCCCACTGGAACGACTGTCCCGCATC B1:Z:678 B2:Z:678
        header_tokens = name.split()
        assert(len(header_tokens) > 3)
        assert(len(sequence) >= umi_end)
        # Get the X and Y tags from the header of the read
        x = header_tokens[-2]
        y = header_tokens[-1]
        # The UMI is retrieved from the sequence
        umi = sequence[umi_start:umi_end]
        # We keep the same naming convention for the UMI attribute
        tags = (x,y,"B3:Z:{}".format(umi))
        # In order to save memory we truncate the read
        # name to only keep the unique part (lane, tile, x_pos, y_pos)
        # TODO this procedure is specific to only Illumina technology
        key = "".join(header_tokens[0].split(":")[-4:])
        hash_reads[key] = tags
        
    if low_memory: hash_reads.commit()
    fastq_file.close()    
    return hash_reads