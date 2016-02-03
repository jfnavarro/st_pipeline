#!/usr/bin/env python
""" 
This module contains some functions to deal with fastq files
"""

from stpipeline.common.utils import *
from stpipeline.common.adaptors import removeAdaptor
import logging 
from stpipeline.common.stats import Stats
from itertools import izip
from blist import sorteddict
import shelve
try:
    import gdbm
    found_gdbm = True
except ImportError:
    found_gdbm = False
    
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
    num_bases_trimmed = len(sequence)
    name = record[0]
    phred = 64 if qual64 else 33
    
    # Get the position at which to trim
    cut_index = quality_trim_index(qscore, min_qual, phred)
    # Get the number of bases suggested to trim
    nbases = num_bases_trimmed - cut_index
    
    # Check if the trimmed sequence would have min length (at least)
    # if so return the trimmed read otherwise return None
    if (num_bases_trimmed - nbases) >= min_length:
        new_seq = record[1][:(num_bases_trimmed - nbases)]
        new_qual = record[2][:(num_bases_trimmed - nbases)]
        return name, new_seq, new_qual
    else:
        return None
    
def getFake(header, num_bases):
    """ 
    Generates a fake fastq record(name,seq,qual) from header and length
    given as input
    """
    new_seq = "".join("N" for k in xrange(num_bases))
    new_qual = "".join("B" for k in xrange(num_bases))
    return header, new_seq, new_qual

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
  
#TODO optimize an refactor this (maybe use reg-exp)
def check_umi_template(umi, template):
    """
    Checks that the UMI given as input complies
    with the pattern given in template
    :param umi a UMI from a read
    :param template a reg-based template with the same
    distance of the UMI that should tell how the UMI should
    look.
    The functions returns true if the UMI complies
    """
    assert(len(umi) == len(template))

    for i,ele in enumerate(template):
        umi_base = umi[i]
        if ele == "W":
            if umi_base != "A" and umi_base != "T":
                return False
        elif ele == "S":
            if umi_base != "C" and umi_base != "G":
                return False
        elif ele == "N":
            if umi_base != "A" and umi_base != "T" and umi_base != "C" and umi_base != "G":
                return False
        elif ele == "V":
            if umi_base == "T":
                return False
        elif ele == "A":
            if umi_base != "A":
                return False
        elif ele == "C":
            if umi_base != "C":
                return False
        elif ele == "G":
            if umi_base != "G":
                return False
        elif ele == "T":
            if umi_base != "T":
                return False
        elif ele == "U":
            if umi_base != "U":
                return False
        elif ele == "R":
            if umi_base != "A" and umi_base != "G":
                return False
        elif ele == "Y":
            if umi_base != "C" and umi_base != "T":
                return False
        elif ele == "K":
            if umi_base != "G" and umi_base != "T":
                return False
        elif ele == "M":
            if umi_base != "A" and umi_base != "C":
                return False
        elif ele == "B":
            if umi_base == "A":
                return False
        elif ele == "D":
            if umi_base == "C":
                return False
        elif ele == "H":
            if umi_base == "G":
                return False
        else:
            return False
            
    return True

#TODO Cythonize this 
def reformatRawReads(fw, 
                     rw,
                     qa_stats,
                     barcode_start=0, 
                     barcode_length=18,
                     filter_AT_content=90,
                     molecular_barcodes=False, 
                     mc_start=18, 
                     mc_end=27,
                     trim_rw=5,
                     min_qual=20, 
                     min_length=28,
                     polyA_min_distance=0, 
                     polyT_min_distance=0, 
                     polyG_min_distance=0, 
                     polyC_min_distance=0,
                     qual64=False, 
                     outputFolder=None, 
                     keep_discarded_files=False,
                     umi_filter=False,
                     umi_filter_template="WSNNWSNNV"):
    """ 
    :param fw the fastq file with the forward reads
    :param rw the fastq file with the reverse reads
    :param barcode_start the base where the barcode sequence starts
    :param barcode_length the number of bases of the barcodes
    :param molecular_barcodes if True the forward reads contain molecular barcodes
    :param mc_start the start position of the molecular barcodes if any
    :param mc_end the end position of the molecular barcodes if any
    :param trim_rw how many bases we want to trim (not consider in the reverse)
    :param min_qual the min quality value to use to trim quality
    :param min_length the min valid length for a read after trimming
    :param polyA_min_distance if >0 we remove PolyA adaptors from the reads
    :param polyT_min_distance if >0 we remove PolyT adaptors from the reads
    :param polyG_min_distance if >0 we remove PolyG adaptors from the reads
    :param qual64 true of qualities are in phred64 format
    :param outputFolder optional folder to output files
    :param keep_discarded_files when true files containing the discarded reads will be created
    :param umi_filter performs a UMI quality filter when True
    :param umi_filter_template the template to use for the UMI filter
    This function does three things (all here for speed optimization)
      - It performs a sanity check
      - It performs a BWA quality trimming discarding very short reads
      - It removes adaptors from the reads (optional)
    """
    logger = logging.getLogger("STPipeline")
    # Create file handlers
    out_rw = 'R2_trimmed_formated.fastq'
    out_rw_discarded = 'R2_trimmed_formated_discarded.fastq'
    if outputFolder and os.path.isdir(outputFolder):
        out_rw = os.path.join(outputFolder, out_rw)
        out_rw_discarded = os.path.join(outputFolder, out_rw_discarded)
    out_rw_handle = safeOpenFile(out_rw, 'w')
    out_rw_writer = writefq(out_rw_handle)
    if keep_discarded_files:
        out_rw_handle_discarded = safeOpenFile(out_rw_discarded, 'w')
        out_rw_writer_discarded = writefq(out_rw_handle_discarded)
    
    # Some counters
    total_reads = 0
    dropped_rw = 0
    dropped_umi = 0
    
    # Build fake sequence adaptors with the parameters given
    adaptorA = "".join("A" for k in xrange(polyA_min_distance))
    adaptorT = "".join("T" for k in xrange(polyT_min_distance))
    adaptorG = "".join("G" for k in xrange(polyG_min_distance))
    adaptorC = "".join("C" for k in xrange(polyG_min_distance))
    
    # Check if barcode settings are correct
    iscorrect_mc = molecular_barcodes
    if mc_start < (barcode_start + barcode_length) \
    or mc_end < (barcode_start + barcode_length):
        logger.warning("Your UMI sequences overlap with the barcodes sequences")
        iscorrect_mc = False
    
    # Open fastq files with the fastq parser
    fw_file = safeOpenFile(fw, "rU")
    rw_file = safeOpenFile(rw, "rU")
    # Iterate records in paralel (Assumes they are ordered by name)
    for line1, line2 in izip(readfq(fw_file), readfq(rw_file)):
        
        if not line1 or not line2:
            logger.error("The input files %s,%s are not of the same length" % (fw,rw))
            break
        
        header_fw = line1[0]
        header_rv = line2[0]
        sequence_fw = line1[1]
        sequence_rv = line2[1]
        num_bases_rv = len(sequence_rv)
        quality_rv = line2[2]
        
        if header_fw.split()[0] != header_rv.split()[0]:
            logger.warning("Pair reads found with different names %s and %s" % (header_fw,header_rv))
            
        # Increase reads counter
        total_reads += 1
        
        # If we want to check for UMI quality and the UMI is incorrect
        # we discard the reads
        if iscorrect_mc and umi_filter \
        and not check_umi_template(sequence_fw[mc_start:mc_end], umi_filter_template):
            dropped_umi += 1
            line2 = None
                                                      
        # If read - trimming is not long enough or has a high AT content discard...
        if (num_bases_rv - trim_rw) < min_length or \
        ((sequence_rv.count("A") + sequence_rv.count("T")) / num_bases_rv) * 100 >= filter_AT_content:
            line2 = None
              
        # if indicated we remove the adaptor PolyA from reverse reads
        if polyA_min_distance > 0 and line2: line2 = removeAdaptor(line2, adaptorA, trim_rw, "5") 
        # if indicated we remove the adaptor PolyT from reverse reads
        if polyT_min_distance > 0 and line2: line2 = removeAdaptor(line2, adaptorT, trim_rw, "5") 
        # if indicated we remove the adaptor PolyG from reverse reads
        if polyG_min_distance > 0 and line2: line2 = removeAdaptor(line2, adaptorG, trim_rw, "5") 
        # if indicated we remove the adaptor PolyC from reverse reads
        if polyC_min_distance > 0 and line2: line2 = removeAdaptor(line2, adaptorC, trim_rw, "5") 
            
        # Trim reverse read
        line2 = trim_quality(line2, trim_rw, min_qual, min_length, qual64) if line2 else None
        
        # Write reverse read to output
        if line2:
            out_rw_writer.send(line2)
        else:
            dropped_rw += 1  
            if keep_discarded_files:
                # Write to discarded the original record
                out_rw_writer_discarded.send((header_rv, sequence_rv, quality_rv))
    
    out_rw_handle.close()
    if keep_discarded_files: out_rw_handle_discarded.close()
    fw_file.close()
    rw_file.close()
    
    if not fileOk(out_rw):
        error = "Error reformatting raw reads: output file not present %s" % (out_rw)
        logger.error(error)
        raise RuntimeError(error + "\n")
    else:
        logger.info("Trimming stats total reads (pair): %s" % (str(total_reads)))
        logger.info("Trimming stats reverse: %s reads have been dropped!" % (str(dropped_rw))) 
        perc2 = '{percent:.2%}'.format(percent= float(dropped_rw) / float(total_reads) )
        logger.info("Trimming stats reverse: you just lost about %s of your data" % (perc2))
        logger.info("Trimming stats reads (reverse) remaining: %s" % (str(total_reads - dropped_rw)))
        if umi_filter:
            logger.info("Trimming stats dropped pairs due to incorrect UMI: %s" % (str(dropped_umi)))
            
        # Adding stats to QA Stats object
        qa_stats.input_reads_forward = total_reads
        qa_stats.input_reads_reverse = total_reads
        qa_stats.reads_after_trimming_forward = 0
        qa_stats.reads_after_trimming_reverse = total_reads - dropped_rw
        
    return out_rw

def hashDemultiplexedReads(reads,
                           has_umi,
                           umi_start,
                           umi_end,
                           low_memory):
    """
    @param reads fastq reads after demultiplexing
    @param has_umi True if the read sequence contains UMI
    @param umi_start the start position of the UMI
    @param umi_end the end position of the UMI
    @param low_memory True to use a key-value db instead of dict
    This function extracts the read name and the barcode
    from the reads given as input and returns a hash
    with the clean read name as key and (barcode,x,y,umi) as
    values (umi is optional)
    """
    if low_memory and found_gdbm:
        #TODO use a global static name variable
        hash_reads = shelve.Shelf(gdbm.open("st_temp_hash_demux", "n")) 
    else:
        hash_reads = sorteddict()

    fastq_file = safeOpenFile(reads, "rU")
    for line in readfq(fastq_file):
        # Assumes the header ends like this B0:Z:GTCCCACTGGAACGACTGTCCCGCATC B1:Z:678 B2:Z:678
        header_tokens = line[0].split()
        assert(len(header_tokens) == 5) 
        barcode = header_tokens[-3]
        # TODO double check that the order of X and Y is correct
        x = header_tokens[-2]
        y = header_tokens[-1]
        # Assumes STAR will only output the first token of the read name
        # We keep the same naming for the extra attributes
        tags = [barcode,x,y]
        if has_umi:
            # Add the UMI as an extra tag
            umi = line[1][umi_start:umi_end]
            tags.append("B3:Z:%s" % umi)
        hash_reads[header_tokens[0]] = tags
    fastq_file.close()
    return hash_reads