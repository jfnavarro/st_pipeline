import os
import sys
import time
import re
import signal
import multiprocessing
import logging
import pysam
import operator
from collections import defaultdict
from pympler.asizeof import asizeof
from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
from stpipeline.common.gff_reader import gff_lines

class geneBuffer():
    """
    This object defines a buffer by holding a dictionary 
    of genes, spots coordinates and transcripts
    it assumes the transcripts to be added in a coordinate ordered fashion
    """
    def __init__(self, gff_filename):

        # defining variables
        cdef dict buffer = dict()
        cdef int last_position = 0
        cdef str last_chromosome = 'chrom'

        # binding variables to self
        self.buffer = buffer
        self.last_position = last_position
        self.last_chromosome = last_chromosome
        self.__compute_gene_end_coordinates(gff_filename)
        
    def __compute_gene_end_coordinates(self, gff_filename):
        """
        function that reads the end coordinate and chromosomes of all the genes present
        in the GFF file and save them as values of a dictionary with the gene ID as key
        dict: [GENE_id] => (chromosome, gene_end_coordinate)
        """
        # Create a dict with end coordinates for all genes
        # to be able to tell when a gene is "processed"
        # and thereby ready to be returned to the parent process
        cdef dict gene_end_coordinates = dict()
        cdef dict line
        cdef str gene_id
        cdef str seqname
        cdef int end

        # parse GTF file
        for line in gff_lines(gff_filename):
            seqname = line['seqname']
            end = int(line['end'])
            gene_id = line['gene_id']
            # save gene_id and rightmost genomic coordinate of each gene to dictionary
            try:
                if gene_id[0] == '"' and gene_id[-1] == '"': 
                    gene_id = gene_id[1:-1]
            except KeyError:
                raise ValueError(
                    'The gene_id attribute is missing in the annotation file ({0})\n'.format(self.gff_filename)
                    )
            try:
                if end > gene_end_coordinates[gene_id][1]:
                    gene_end_coordinates[gene_id] = (seqname, end)
            except KeyError:
                gene_end_coordinates[gene_id] = (seqname, end)

        # A fix to include any "__no_feature" annotations
        gene_end_coordinates['__no_feature'] = (None,-1)
        self.gene_end_coordinates = gene_end_coordinates

    def get_gene_end_position(self, gene):
        """
        Returns the genomic end coordinate and chromosome
        of the gene given as input
        """
        cdef tuple gene_end_coordinate
        cdef list ambiguous_genes
        try:
            gene_end_coordinate = self.gene_end_coordinates[gene]
        # if the gene annotation is not found in the GTF file
        # it can either be an ambiguous annotation from htseq or we should raise a error
        # So we try to get the end coordinate for HTSeq ambiguous annotation
        except KeyError:
            # trying to handle HTSeq ambiguous gene annotations correctly,
            # they will have an annotation with the following style __ambiguous[GENE1+GENE2]
            if gene[0:len('__ambiguous[')] == '__ambiguous[':
                try: # to get the right most right most coordinate ;)
                    ambiguous_genes = gene[len('__ambiguous['):-1].split('+')
                    gene_end_coordinate = max([self.gene_end_coordinates[amb_gene]
                                               for amb_gene in ambiguous_genes],
                                               key=operator.itemgetter(1))
                except KeyError:
                    raise ValueError("Error parsing unique events gene with id {0}" \
                                     "is not found in annotation file\n".format(gene))
            else:
                raise ValueError("Error parsing unique events gene with id {0}" \
                                 "was not found in the annotation file\n".format(gene))
        return gene_end_coordinate
                
    def add_transcript(self, gene, spot_coordinates, transcript, position):
        """
        Adds a transcript to the gene buffer
        Parameters:
        :param gene: the name of the gene
        :param spot_coordinates: the spot coordinates as a (x,y) tuple
        :param position: the transcript's lest most genomic coordinate
        :param transcript: the transcript information
            as a (chrom, start, end, clear_name, mapping_quality, strand, umi) tuple
            (i.e. AlignedSegment.reference_start)
        """
        self.last_position = position
        self.last_chromosome = transcript[0]
        cdef list new_reads_list
        # First try to add the "transcript" to the existing buffer
        try:
            self.buffer[gene][spot_coordinates].append(transcript)
        # In case the spot coordinates are not found make a new reads list,
        # add the "transcript" to it and add the spot coordinates to the specified gene
        except KeyError:
            new_reads_list = [transcript]
            try:
                self.buffer[gene][spot_coordinates] = new_reads_list
            except KeyError:
                # In case the gene is not in the buffer make a new gene dictionary 
                # with spot coordinates and read lists and add it to the gene buffer dictionary
                self.buffer[gene] = {spot_coordinates:new_reads_list}

    def check_and_clear_buffer(self, empty=False):
        """
        Iterates the genes in the buffer to check
        if any of the genes are outside the last inserted coordinates/chromosomes
        If so the gene will be returned in a list and deleted from the buffer
        :param empty: when True if forces to empty the buffer
        """
        cdef list _tmp = self.buffer.keys()
        cdef gene_transcripts = list()
        cdef str chrom
        cdef int end_position
        #NOTE seems like we have to go trough all the genes
        #in the buffer every time we want to check and clear
        #we cannot rely on the gene's genomic coordinates only
        #even though the BAM file is sorted by position
        for gene in _tmp:
            # fix to include any "__no_feature" annotations
            if gene == '__no_feature' and not empty: 
                continue
            #TODO better to cache these in a dictionary
            (chrom, end_position) = self.get_gene_end_position(gene)
            # check if the current position is past the gene end coordinate
            if empty or self.last_position > end_position or self.last_chromosome != chrom:
                yield (gene, self.buffer[gene])
                # Remove the gene from the buffer
                del self.buffer[gene]
        return gene_transcripts
                
def parse_unique_events(input_file, gff_filename=None):
    """
    This function parses the transcripts present in the filename given as input.
    It expects a coordinate sorted BAM file where the spot coordinates,
    gene and UMI are present as extra tags (optionally)
    Will yield a dictionary per gene with a spot coordinate tuple as keys
    foreach gene yield: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ]
    :param filename: the input file containing the annotated BAM records
    :param gff_filename: the gff file containing the gene coordinates (optional)
    """
    cdef object genes_buffer = geneBuffer(gff_filename) if gff_filename is not None else None
    cdef object genes_dict = dict()
    cdef object rec
    cdef str clear_name
    cdef int mapping_quality
    cdef int start
    cdef int end
    cdef str chrom
    cdef str strand
    cdef tuple transcript
    cdef tuple spot_coordinates
    cdef int x
    cdef int y
    cdef str gene

    # Open the log file and open the bam file for reading
    sam_file = pysam.AlignmentFile(input_file, "rb")
    
    # Parse the coordinate sorted bamfile record by record i.e. by genome 
    # coordinate from first chromosome to last
    for rec in sam_file.fetch(until_eof=True):

        # get the info about the record from the bam file
        clear_name = rec.query_name
        mapping_quality = rec.mapping_quality

        # Account for soft-clipped bases when retrieving the start/end coordinates
        start = int(rec.reference_start - rec.query_alignment_start)
        end = int(rec.reference_end + (rec.query_length - rec.query_alignment_end))
        chrom = sam_file.getrname(rec.reference_id)
        strand = "+"
        if rec.is_reverse:
            # We swap start and end if the transcript mapped to the reverse strand
            strand = "-"
            start, end = end, start

        # Get TAGGD tags from the bam file
        x,y,gene,umi = (-1,-1,'None','None')
        for (k, v) in rec.tags:
            if k == "B1":
                x = int(v) ## The X coordinate
            elif k == "B2":
                y = int(v) ## The Y coordinate
            elif k == "XF":
                gene = str(v) ## The gene name
            elif k == "B3":
                umi = str(v) ## The UMI
            else:
                continue

        # Create a new transcript and add it to the in memory gene_buffer dictionary
        transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
        if gff_filename is not None:
            genes_buffer.add_transcript(gene, (x,y), transcript, rec.reference_start)
            for g, t in genes_buffer.check_and_clear_buffer():
                yield (g, t)
        else:
            try:
                genes_dict[gene].append(transcript)
            except KeyError:
                genes_dict[gene] = [transcript]
        
    # Close the bam file and yield the last gene(s)
    sam_file.close()
    if gff_filename is not None:
        for g, t in genes_buffer.check_and_clear_buffer(True):
            yield (g, t)
    else:
        return genes_dict
