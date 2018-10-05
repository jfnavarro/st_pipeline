import sys
import os
import re
import pysam
import ctypes
import logging
from stpipeline.common.utils import safeOpenFile, fileOk, is_fifo
from stpipeline.common.fastq_utils import *
from stpipeline.common.sam_utils import convert_to_AlignedSegment
from stpipeline.common.adaptors import removeAdaptor
from stpipeline.common.stats import qa_stats
from itertools import izip

bam_header = {
        'HD': {'VN': '1.5', 'SO':'unsorted'},
        'RG': [{'ID': '0', 'SM' : 'unknown_sample', 'PL' : 'ILLUMINA' }]
    }
        
def InputReadsFilter(fw,
                     rv,
                     out_rv,
                     out_rv_discarded,
                     barcode_length,
                     start_position,
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
                     polyN_min_distance,
                     qual64,
                     umi_filter,
                     umi_filter_template,
                     umi_quality_bases,
                     adaptor_missmatches,
                     overhang,
                     disable_umi,
                     disable_barcode):
    """
    This class handles the input read filtering and quality trimming
      - It performs a sanity check (forward and reverse reads same length and order)
      - It performs a BWA-based quality trimming discarding very short reads
      - It removes adaptors from the reads (optional)
      - It checks for AT and GC content (optional)
      - It performs a sanity check on the UMI (optional)
    Reads that do not pass the filters are discarded (both R1 and R2)
    Reads that pass the filter are written as BAM (R2)
    :param fw: the FASTQ file with the forward reads (R1)
    :param rv: the FASTQ file with the reverse reads (R2)
    :param out_rv: the name of the output file for reverse (BAM)
    :param out_rv_discarded: the name of the output file for discarded reverse reads (FASTQ)
    :param barcode_length: length of the barcode sequence (integer)
    :param start_position: the start position of the barcode
    :param filter_AT_content: % of A and T bases a read2 must have to be discarded
    :param filter_GC_content: % of G and C bases a read2 must have to be discarded
    :param umi_start: the start position of the UMI
    :param umi_end: the end position of the UMI
    :param min_qual: the min quality value to use in the trimming
    :param min_length: the min valid length for a read after trimming
    :param polyA_min_distance: if >5 remove PolyA adaptors from the reads
    :param polyT_min_distance: if >5 remove PolyT adaptors from the reads
    :param polyG_min_distance: if >5 remove PolyG adaptors from the reads
    :param polyC_min_distance: if >5 remove PolyC adaptors from the reads
    :param polyN_min_distance: if >5 remove PolyN adaptors from the reads
    :param qual64: true of qualities are in phred64 format
    :param umi_filter: performs a UMI quality template filter when True
    :param umi_filter_template: the template to use for the UMI filter
    :param umi_quality_bases: the number of low quality bases allowed in an UMI
    :param adaptor_missmatches: number of miss-matches allowed when removing adaptors
    :param overhang: overhang to be used for the barcodes (integer)
    :param disable_umi: true if the reads do not contain UMIs
    :param disable_barcode: true if the reads do not contain barcodes
    """
    logger = logging.getLogger("STPipeline")
    if not (os.path.isfile(fw) or is_fifo(fw)) or not (os.path.isfile(rv) or is_fifo(rv)):
        error = "Error doing quality trimming, input file/s not present {}\n{}\n".format(fw,rv)
        logger.error(error)
        raise RuntimeError(error)

    # Check if discarded files must be written out
    cdef bint keep_discarded_files = out_rv_discarded is not None

    # Build fake sequence adaptors with the parameters given
    cdef str adaptorA = "".join("A" for k in xrange(polyA_min_distance))
    cdef str adaptorT = "".join("T" for k in xrange(polyT_min_distance))
    cdef str adaptorG = "".join("G" for k in xrange(polyG_min_distance))
    cdef str adaptorC = "".join("C" for k in xrange(polyC_min_distance))
    cdef str adaptorN = "".join("N" for k in xrange(polyN_min_distance))

    # Not recommended to do adaptor trimming for adaptors smaller than 5
    cdef bint do_adaptorA = polyA_min_distance >= 5
    cdef bint do_adaptorT = polyT_min_distance >= 5
    cdef bint do_adaptorG = polyG_min_distance >= 5
    cdef bint do_adaptorC = polyC_min_distance >= 5
    cdef bint do_adaptorN = polyN_min_distance >= 5
    cdef bint do_AT_filter = filter_AT_content > 0
    cdef bint do_GC_filter = filter_GC_content > 0

    # Quality format
    cdef int phred = 64 if qual64 else 33
      
    # Some counters
    cdef int total_reads = 0
    cdef int dropped_umi = 0
    cdef int dropped_umi_template = 0
    cdef int dropped_AT = 0
    cdef int dropped_GC = 0
    cdef int dropped_adaptor = 0
    cdef int too_short_after_trimming = 0
          
    # Some variables to avoid overhead in the loop
    cdef str header_fw
    cdef str sequence_fw
    cdef str quality_fw
    cdef str header_rv
    cdef str sequence_rv
    cdef str quality_rv
    cdef str orig_sequence_rv
    cdef str orig_quality_rv
    cdef bint discard_read
    
    # Create output file writers
    bam_file = pysam.AlignmentFile(out_rv, "wbu", header=bam_header)
    fw_file = safeOpenFile(fw, "rU")
    rv_file = safeOpenFile(rv, "rU")
    if keep_discarded_files:
        out_rv_handle_discarded = safeOpenFile(out_rv_discarded, 'w')
        out_rv_writer_discarded = writefq(out_rv_handle_discarded)
        
    for (header_fw, sequence_fw, quality_fw), \
    (header_rv, sequence_rv, quality_rv) in izip(readfq(fw_file), readfq(rv_file)):
        
        discard_read = False
        orig_sequence_rv, orig_quality_rv = sequence_rv, quality_rv
        total_reads += 1
        
        if not sequence_fw or not sequence_rv:
            error = "Error doing quality trimming.\n" \
            "The input files are not of the same length"
            logger.error(error)
            break

        if header_fw.split()[0] != header_rv.split()[0]:
            logger.warning("Pair reads found with different " \
                                "names {} and {}".format(header_fw,header_rv))

        # get the barcode sequence
        if disable_barcode:
            barcode = None
        else:
            barcode = sequence_fw[max(0,start_position-overhang):(start_position+barcode_length+overhang)]

        if not disable_umi:
            # If we want to check for UMI quality and the UMI is incorrect
            # then we discard the reads
            umi_seq = sequence_fw[umi_start:umi_end]
            if umi_filter \
            and not check_umi_template(umi_seq, umi_filter_template):
                dropped_umi_template += 1
                discard_read = True

            # Check if the UMI has many low quality bases
            umi_qual = quality_fw[umi_start:umi_end]
            if not discard_read and (umi_end - umi_start) >= umi_quality_bases and \
            len([b for b in umi_qual if (ord(b) - phred) < min_qual]) > umi_quality_bases:
                dropped_umi += 1
                discard_read = True
        else:
            umi_seq = None

        # If reverse read has a high AT content discard...
        if not discard_read and do_AT_filter and \
        ((sequence_rv.count("A") + sequence_rv.count("T")) / len(sequence_rv)) * 100 >= filter_AT_content:
            dropped_AT += 1
            discard_read = True

        # If reverse read has a high GC content discard...
        if not discard_read and do_GC_filter and \
        ((sequence_rv.count("G") + sequence_rv.count("C")) / len(sequence_rv)) * 100 >= filter_GC_content:
            dropped_GC += 1
            discard_read = True

        if not discard_read:
            # Perform adaptor/homopolymer filters
            if do_adaptorA and len(sequence_rv) > min_length:
                sequence_rv, quality_rv = removeAdaptor(
                    sequence_rv, quality_rv, adaptorA, adaptor_missmatches)
            if do_adaptorT and len(sequence_rv) > min_length:
                sequence_rv, quality_rv = removeAdaptor(
                    sequence_rv, quality_rv, adaptorT, adaptor_missmatches)
            if do_adaptorG and len(sequence_rv) > min_length:
                sequence_rv, quality_rv = removeAdaptor(
                    sequence_rv, quality_rv, adaptorG, adaptor_missmatches)
            if do_adaptorC and len(sequence_rv) > min_length:
                sequence_rv, quality_rv = removeAdaptor(
                    sequence_rv, quality_rv, adaptorC, adaptor_missmatches)
            if do_adaptorN and len(sequence_rv) > min_length:
                sequence_rv, quality_rv = removeAdaptor(
                    sequence_rv, quality_rv, adaptorN, adaptor_missmatches)

            # Check if the read is smaller than the minimum after removing artifacts
            if len(sequence_rv) < min_length:
                dropped_adaptor += 1
                discard_read = True

        if not discard_read:
            # Trim reverse read (will return None if length of trimmed sequence is less than min_length)
            sequence_rv, quality_rv = trim_quality(
                sequence_rv,
                quality_rv,
                min_qual,
                min_length,
                phred)
            if not sequence_rv or not quality_rv:
                too_short_after_trimming += 1
                discard_read = True

        if not discard_read:
            bam_file.write(
                convert_to_AlignedSegment(
                    header_rv,
                    sequence_rv,
                    quality_rv,
                    barcode,
                    umi_seq))
        elif keep_discarded_files:
            out_rv_writer_discarded.send((header_rv, orig_sequence_rv, orig_quality_rv))

    bam_file.close()                    
    fw_file.close()
    rv_file.close()
    if keep_discarded_files:
        out_rv_handle_discarded.flush()
        out_rv_handle_discarded.close()
        out_rv_writer_discarded.close()
                
    # Write info to the log
    cdef int dropped_rv = dropped_umi + dropped_umi_template + \
                          dropped_AT + dropped_GC + dropped_adaptor + \
                          too_short_after_trimming
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
    logger.info("Trimming stats dropped pairs due to being too short: {}".format(too_short_after_trimming))
    
    # Check that output file was written ok
    cdef int remaining_reads = (total_reads - dropped_rv)
    if not fileOk(out_rv) or remaining_reads == 0:
        error = "Error doing quality trimming checks of raw reads." \
        "\nOutput file not present {}\n".format(out_rw)
        logger.error(error)
        raise RuntimeError(error)
    
    # Adding stats to QA Stats object
    qa_stats.input_reads_forward = total_reads
    qa_stats.input_reads_reverse = total_reads
    qa_stats.reads_after_trimming_forward = remaining_reads
    qa_stats.reads_after_trimming_reverse = remaining_reads
