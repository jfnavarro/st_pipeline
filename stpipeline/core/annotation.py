""" 
This module contains a modified version of htseq-count
with slight modifications to perform annotation
of ST mapped reads against a reference.
"""
import logging
import os
import pysam
from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import itertools
import HTSeq
import multiprocessing
import operator
import math
import time


class UnknownChrom( Exception ):
    pass

def invert_strand( iv ):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError, "Illegal strand"  
    return iv2

def count_reads_in_features(sam_filename, 
                            gff_filename, 
                            samtype, 
                            order, 
                            stranded,
                            overlap_mode, 
                            feature_type, 
                            id_attribute, 
                            quiet, 
                            minaqual, 
                            samout, 
                            include_non_annotated, 
                            htseq_no_ambiguous,
                            outputDiscarded):
    """
    This is taken from the function count_reads_in_features() from the 
    script htseq-count in the HTSeq package version 0.61.p2 
    The reason to do so is to fix two really small bugs related to the SAM output.
    The code of the function is small and simple so for now we
    will use the patched function here. A patch request has been sent
    to the HTSeq team.
    The description of the parameters are the same as htseq-count.
    Two parameters were added to filter out what to write in the sam output
    
    The HTSEQ License
    HTSeq is free software: you can redistribute it and/or modify it under the terms of 
    the GNU General Public License as published by the Free Software Foundation, 
    either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    The full text of the GNU General Public License, version 3, 
    can be found here: http://www.gnu.org/licenses/gpl-3.0-standalone.html
    """
    # Set up the filters
    count_reads_in_features.filter_htseq = \
    ["__too_low_aQual", "__not_aligned", "__alignment_not_unique"]
    if not include_non_annotated: 
        count_reads_in_features.filter_htseq.append("__no_feature")
    count_reads_in_features.filter_htseq_no_ambiguous = htseq_no_ambiguous

    # Open SAM/BAM output file
    flag_write = "wb" if samtype == "bam" else "wh"
    flag_read = "rb" if samtype == "bam" else "r"
    saminfile = pysam.AlignmentFile(sam_filename, flag_read)
    count_reads_in_features.samoutfile = pysam.AlignmentFile(samout, flag_write, template=saminfile)
    if outputDiscarded is not None:
        count_reads_in_features.samdiscarded = pysam.AlignmentFile(outputDiscarded, flag_write, template=saminfile)
    saminfile.close()
    
    # Counter of annotated records
    count_reads_in_features.annotated = 0
    
    # Function to write to SAM output
    def write_to_samout(read, assignment):
        # Creates the PySAM record
        # to_pysam_AlignedSegment is the new method in HTSeq>=0.7.0 that
        # uses the latest Pysam API and reports the correct sequences
        sam_record = read.to_pysam_AlignedSegment(count_reads_in_features.samoutfile)
        sam_record.set_tag("XF", assignment, "Z")
        if read is not None and assignment not in count_reads_in_features.filter_htseq \
        and not (count_reads_in_features.filter_htseq_no_ambiguous and assignment.find("__ambiguous") != -1):
            count_reads_in_features.samoutfile.write(sam_record)
            count_reads_in_features.annotated += 1
        elif outputDiscarded is not None:
            count_reads_in_features.samdiscarded.write(sam_record)
                
    # Annotation objects
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no")
    counts = {}
    gff = HTSeq.GFF_Reader(gff_filename)   

    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError, ("Feature %s does not contain a '%s' attribute" \
                                       % (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError, ("Feature %s at %s does not have strand information but you are " \
                                       "running htseq-count in stranded mode. Use '--stranded=no'." % 
                                       (f.name, f.iv))
                features[f.iv] += feature_id
                counts[f.attr[id_attribute]] = 0
    except:
        raise
    
    if len(counts) == 0:
        raise RuntimeError, "No features of type '%s' found.\n" % feature_type
        
    if samtype == "sam":
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
    elif samtype == "bam":
        SAM_or_BAM_Reader = HTSeq.BAM_Reader
    else:
        raise ValueError, "Unknown input format %s specified." % samtype

    try:
        read_seq = SAM_or_BAM_Reader(sam_filename)
    except:
        raise RuntimeError, "Error occurred when reading beginning of SAM/BAM file."

    try:
        
        for r in read_seq:
            if not r.aligned:
                write_to_samout(r, "__not_aligned")
                continue
            try:
                if r.optional_field("NH") > 1:
                    write_to_samout(r, "__alignment_not_unique")
                    # Should these reads potentially be printed twice?
                    # should there not be a continue statement here?
                    # otherwise the read will move on through the if statemets
                    # until it gets a gene id annotation and will be printed again?
            except KeyError:
                pass
            if r.aQual < minaqual:
                write_to_samout(r, "__too_low_aQual")
                continue
            if stranded != "reverse":
                iv_seq = (co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0)
            else:
                iv_seq = (invert_strand(co.ref_iv) 
                          for co in r.cigar 
                          if co.type == "M" and co.size > 0)
            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            fs = fs.union(fs2)
                elif overlap_mode == "intersection-strict" or overlap_mode == "intersection-nonempty":
                    fs = None
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            if len(fs2) > 0 or overlap_mode == "intersection-strict":
                                if fs is None:
                                    fs = fs2.copy()
                                else:
                                    fs = fs.intersection(fs2)
                else:
                    raise RuntimeError, "Illegal overlap mode."
                
                if fs is None:
                    continue
                elif len(fs) == 0:
                    write_to_samout(r, "__no_feature")
                elif len(fs) > 1:
                    write_to_samout(r, "__ambiguous[" + '+'.join(fs) + "]")
                else:
                    write_to_samout(r, list(fs)[0])
                    
            except UnknownChrom:
                pass

    except:
        count_reads_in_features.samoutfile.close()
        if outputDiscarded is not None:
            count_reads_in_features.samdiscarded.close()
        raise

    count_reads_in_features.samoutfile.close()
    if outputDiscarded is not None:
        count_reads_in_features.samdiscarded.close()
    return count_reads_in_features.annotated

def annotateReads(mappedReads, 
                  gtfFile,
                  outputFile,
                  outputDiscarded,
                  mode,
                  strandness,
                  htseq_no_ambiguous, 
                  include_non_annotated):
    """
    Annotates a file with mapped reads (BAM) using a modified 
    version of the htseq-count tool. It writes the annotated records to a file.
    It assumes the input reads (BAM) are single end and do not contain
    multiple alignments or un-annotated reads.
    :param mappedReads: path to a BAM file with mapped reads sorted by coordinate
    :param gtfFile: path to an annotation file in GTF format
    :param outputFile: where to write the annotated records (BAM)
    :param outputDiscarded: where to write the non-annotated records (BAM)
    :param mode: htseq-count overlapping mode (see htseq-count documentation)
    :param strandness: the type of strandness to use when annotating (yes, no or reverse)
    :param htseq_no_ambiguous: true if we want to discard ambiguous annotations
    :param include_non_annotated: true if we want to include 
    non annotated reads as __no_feature in the output
    :param outputFile: the name/path to the output file
    :type mappedReads: str
    :type gtfFile: str
    :type outputFile: str
    :type outputDiscarded: str
    :type mode: str
    :type strandness: str
    :type htseq_no_ambiguos: boolean
    :type include_non_annotated: str
    :type outputFile: str
    :raises: RuntimeError, ValueError
    """
    
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(mappedReads):
        error = "Error during annotation, input file not present {}\n".format(mappedReads)
        logger.error(error)
        raise RuntimeError(error)
    
    try:

        # split bam file
        partial_bam_files = split_bam(mappedReads)

        # create an annotation subprocess for each partial bam
        annotated = 0
        subprocesses = [ multiprocessing.Process(
            target=count_reads_in_features,
            args=(
                bamfile_name,
                gtfFile,
                "bam", # Type BAM for files
                "pos", # Order pos or name
                strandness, # Strand yes/no/reverse
                mode, # intersection_nonempty, union, intersection_strict
                "exon", # feature type in GFF
                "gene_id", # gene_id or gene_name
                True, # Quiet mode
                0, # Min quality score
                "{0}.part_{1}.bam".format(outputFile,part),
                include_non_annotated,
                htseq_no_ambiguous,
                "{0}.part_{1}.bam".format(outputDiscarded,part) if outputDiscarded else outputDiscarded
            )) for part, bamfile_name in partial_bam_files.iteritems() ]

        # start work in child processes
        for p in subprocesses: p.start()

        # wait for all processes to finish
        while True in [p.is_alive() for p in subprocesses]: time.sleep(0.1)

        # join the children
        for p in subprocesses:
            assert p.exitcode == 0, "ERROR:: subprocess error during annotation."
            p.join()

        # merge the annotated bam files and summmarize the stats
        annotations=merge_bam( outputFile, ["{0}.part_{1}.bam".format(outputFile,part) for part in partial_bam_files] )
        if outputDiscarded:
            discarded_annotations=merge_bam( outputDiscarded, ["{0}.part_{1}.bam".format(outputDiscarded,part) for part in partial_bam_files] )

        annotated = sum(annotations.values())

    except Exception as e:
        error = "Error during annotation. HTSEQ execution failed\n"
        logger.error(error)
        raise e

    if not fileOk(outputFile):
        error = "Error during annotation. Output file not present {}\n".format(outputFile)
        logger.error(error)
        raise RuntimeError(error)

    logger.info("Annotated reads: {}".format(annotated))
    qa_stats.reads_after_annotation = annotated

def split_bam(input_bamfile_name):

    # set the number of output parts to equal the number of cores on the machine
    number_of_output_parts = multiprocessing.cpu_count()
    parts = { part:[] for part in range(number_of_output_parts) }
    part_read_counts = { part:0 for part in parts }
    output_file_names = {part:"{0}.part_{1}.bam".format(input_bamfile_name,part) for part in parts}

    # index and open the input_bamfile
    pysam.index(input_bamfile_name,'{0}.bai'.format(input_bamfile_name))
    input_bamfile = pysam.AlignmentFile( input_bamfile_name, mode='rb' )
    assert input_bamfile.check_index()

    # open the output bam files
    output_bamfiles = {
        part:pysam.AlignmentFile(file_name, mode="wbu", template=input_bamfile) \
        for part, file_name in output_file_names.iteritems()
        }

    # check the number of contigs
    index_stats = input_bamfile.get_index_statistics()

    # if the number of contigs exceeds the number of cores save each contig to a file
    # this is faster than splitting based on read count
    if len(index_stats) > number_of_output_parts and False:
        # added "and False" to the if-statement to always split based on equal read counts
        # this due to that the mitochondria sometimes seem to get a lot of reads

        # place the contigs to parts ordered by read count to get an as even distribution of reads as possible
        for contig, contig_read_count in \
        sorted([ (stat[0],stat[3]) for stat in index_stats], key=operator.itemgetter(1), reverse=True):
            part, part_read_count = sorted(
                [ (part,part_read_counts[part]) for part in parts], key=operator.itemgetter(1)
                )[0]
            parts[part].append( contig )
            part_read_counts[part] += contig_read_count

        for part, contigs in parts.iteritems():
            print parts, part_read_counts[part]
            for contig in contigs:
                for record in input_bamfile.fetch(contig=contig): output_bamfiles[part].write( record )

    # if number of contigs are less than cores split on equal read count instead
    else:
        total_read_count = input_bamfile.mapped + input_bamfile.unmapped
        reads_per_part = math.ceil(float(total_read_count)/number_of_output_parts)

        _tmp_read_counter = 0
        part = 0
        for record in input_bamfile.fetch(until_eof=True):
            output_bamfiles[part].write(record)
            _tmp_read_counter += 1
            if _tmp_read_counter == reads_per_part:
                part += 1
                _tmp_read_counter = 0

    return output_file_names

def merge_bam(merged_file_name, files_to_merge):

    annotations = {}

    with pysam.AlignmentFile( files_to_merge[0], mode='rb' ) as input_bamfile:
        merged_file = pysam.AlignmentFile(merged_file_name, mode="wb", template=input_bamfile)

    for file_name in files_to_merge:
        input_bamfile = pysam.AlignmentFile( file_name, mode='rb' )
        for record in input_bamfile.fetch(until_eof=True):
            merged_file.write(record)
            try:
                annotations[record.get_tag("XF")] += 1
            except KeyError:
                annotations[record.get_tag("XF")] = 1

    return annotations