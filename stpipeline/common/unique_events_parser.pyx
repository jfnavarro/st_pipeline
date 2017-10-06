import os
import sys
import time
import re
import signal
import multiprocessing
import logging
import pysam
from collections import defaultdict
from pympler.asizeof import asizeof
from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
from stpipeline.common.gff_reader import gff_lines

class uniqueEventsParser():
    """
    The UniqueEventsParser class is a more memory efficient implementations
    of the stpipeline.comon.sam_utils.parseUniqueEvents-function.
    The class takes two values as input filename and gff_filename,
    as well as two keyword arguments verbose and max_genes_in_memory.

    The max_genes_in_memory variable defines the size of the underlying queue
    used to send genes from the file parsing function to the main process.
    
    The work done by the parser is defined by the self._worker_function function.
    In this function a annotated coordinate sorted bam file is parsed
    and genes are returned to the parent process through a multiprocessing.Queue self.q.

        Parses the transcripts present in the filename given as input.
        It expects a coordinate sorted BAM file where the spot coordinates,
        gene and UMI are present as extra tags.
        And a gtf file defining the coordinates of the genes in the reference genome.
        Will yield put a dictionary per gene with a spot coordinate tuple as keys into the self.q multiprocessing.Queue
        foreach gene put dict: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ]

    :param filename: the input file containing the annotated BAM records
    :param gff_filename: the gff file containing the gene coordinates
    :param verbose: boolean used to determine if info should be written to stdout (default False)
    :param max_genes_in_memory: integer number of genes to keep in memory, i.e. the size of the multiprocessing.Queue
    :return: a instance of the UniqueEventsParser
    """

    def __init__(self, filename, gff_filename, verbose=False, max_genes_in_memory=50):
        self.filename = filename
        self.gff_filename = gff_filename
        self.verbose = verbose
        self.max_genes_in_memory = max_genes_in_memory

    def run(self, ):
        """
        Function that starts the subprocess reading the file(s)
        """

        # create output queue
        self.q = multiprocessing.Queue(self.max_genes_in_memory)

        # start worker subprocess
        self.p = multiprocessing.Process(target=self._worker_function)
        self.p.start()

        if self.verbose: sys.stdout.write('Process={0} starting to parse unique events.\n'.format(str(self.p.pid)))

    def check_running(self, ):
        """
        Function that checks if the subprocess is running
        and return the status as as a string with value "RUNNING" or "COMPLETE"
        """

        # check for completion of worker process
        if self.p.is_alive(): return 'RUNNING'
        else:
            if self.p.exitcode == 0:
                self.p.join()
                if self.verbose: sys.stdout.write(
                    'Process={0} finished parsing unique events.\n'.format(str(self.p.pid))
                    )
                return 'COMPLETE'
            else:
                msg = 'ERROR:: Process={0} FAILED to parse unique events!!!.\n'.format(str(self.p.pid))
                if self.verbose: sys.stdout.write(msg)
                raise RuntimeError(msg)

    def _worker_function(self, ):
        """
        This function is a modified version of the stpipeline.comon.sam_utils.parseUniqueEvents-function:

        Parses the transcripts present in the filename given as input.
        It expects a coordinate sorted BAM file where the spot coordinates,
        gene and UMI are present as extra tags
        Will yield a dictionary per gene with a spot coordinate tuple as keys
        foreach gene yield: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ]
        :param filename: the input file containing the annotated BAM records
        :param gff_filename: the gff file containing the gene coordinates
        """

        if self.verbose: sys.stdout.write('INFO:: ENTERING => parseUniqueEvents_byCoordinate\n')

        # Open the log file and open the bamfile for reading
        logger = logging.getLogger("STPipeline")
        sam_file = pysam.AlignmentFile(self.filename, "rb")

        # define variables cython style for more efficient looping and in memory storage
        cdef object genes_buffer = geneBuffer( verbose=self.verbose )
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

        genes_buffer.get_gene_end_coordinates(self.gff_filename)

        # Parse the coordinate sorted bamfile record by record i.e. by genome coordinate from first chromosome to last
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
            # Check that all tags are present
            if 'None' in [gene,umi] or -1 in [x,y]:
                logger.warning("Warning parsing annotated reads.\n" \
                               "Missing attributes for record {}\n".format(clear_name))
                continue

            # Create a new transcript and add it to the in memory gene_buffer dictionary
            transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
            spot_coordinates = (x,y)
            genes_buffer.add_transcript( gene,spot_coordinates,transcript,rec.reference_start )
            genes_buffer.put_finished_genes_in_queue(self.q)

        # Close the bam file and yield the last gene(s)
        sam_file.close()
        genes_buffer.put_finished_genes_in_queue(self.q, empty=True)
        self.q.put( "COMPLETED" ) # send a last keyword to the parent process to define the end of the data

        # cleanup and exit
        while not self.q.empty():
            if self.verbose:
                sys.stdout.write(
                    'Process={0} Done processing input wainting for queue to empty.\n'.format(str(self.p.pid))
                    )
            time.sleep(1)
        self.q.close()
        self.q.join_thread()
        if self.verbose: sys.stdout.write('Process={0} returning.\n'.format(str(self.p.pid)))

    def all_unique_events(self,):
        """
        A generator function living in the parent process fetching genes from the queue
        and yelding the genes sent back.
        """
        
        # start the worker if not already done
        if self.check_running != 'RUNNING': self.run()
        
        # yield genes until the bam file has been parsed and the uniqueEventsParser shuts down
        while True:
            data = self.q.get()
            if self.check_running != 'COMPLETE' and data == 'COMPLETED':
                if self.verbose: sys.stdout.write('INFO:: got signal {0} from uep.\n'.format(data))
                while self.check_running() != 'COMPLETE':
                    if self.verbose: sys.stdout.write('INFO:: waiting for uep subprocesses to finish.\n')
                    time.sleep(0.1)
                break
            if self.verbose: sys.stdout.write('INFO:: got gene {0} from child process.\n'.format(data[0]))
            yield data

class geneBuffer():
    """
    This object defines a buffer by holding a dictionary of genes,spots coordinates and transcripts
    it assumes the transcripts to be added in a coordinate ordered fashion

    Basic usage:
      First create a insatance of the GeneBuffer
      Then get the end coordinates for each gene using the instance.get_gene_end_coordinates function and a gff file
      Then add transcript(s) using the instance.add_transcript function
      To send "finished" genes to a queue-object use the instance.put_finished_genes_in_queue function
    """

    def __init__(self, verbose=True):

        # defining varibles
        cdef dict _buffer = dict()
        cdef dict processed_genes = dict()
        cdef list deleted_genes = list()
        cdef int _total_transcript_counter = 0
        cdef int _transcripts_sent_counter = 0
        cdef double _speed_of_last_100k_transcripts = 0
        cdef double _start_time=time.time()
        cdef double _time_for_last_100k_print = _start_time
        cdef int _current_position = 0
        cdef str _current_chromosome = 'chrom'

        # binding variables to self
        self.buffer = _buffer
        self.processed_genes = processed_genes
        self.deleted_genes = deleted_genes
        self.current_position = _current_position
        self.current_chromosome = _current_chromosome
        self.total_transcript_counter = _total_transcript_counter
        self.transcripts_sent_counter = _transcripts_sent_counter
        self.speed_of_last_100k_transcripts = _speed_of_last_100k_transcripts
        self.start_time=_start_time
        self.time_for_last_100k_print = _time_for_last_100k_print

        # print options
        self.verbose = verbose
        if self.verbose: self.print_stats( header=True )

    def get_gene_end_coordinates(self, gff_filename):
        """
        function that reads the end coordinates of genes
        and save them as values of a dictionary with the gene ID as key
        dict: [GENE_id] => gene_end_coordinate
        """

        # Get dict with end coordinates for all genes
        # to be able to tell when a gene is "done"
        # and thereby ready to be returned to the parent process


        cdef dict gene_end_coordinates
        cdef dict line
        cdef str gene_id
        cdef str seqname
        cdef int end

        if self.verbose: sys.stdout.write('INFO:: LOADING gtf file...\n')
        gene_end_coordinates = dict()

        # parse gtf file
        for line in gff_lines(gff_filename):

            seqname = line['seqname']
            end = int(line['end'])
            gene_id = line['gene_id']

            # save gene_id and rigtmost genomic coordinate of each gene to dictionary
            try:
                if gene_id[0] == '"' and gene_id[-1] == '"': gene_id=gene_id[1:-1]
            except KeyError:
                raise ValueError(
                    'The gene_id attribute is missing in the annotation file ({0})\n'.format(self.gff_filename)
                    )
            try:
                if int(end) > gene_end_coordinates[ gene_id ][1]:
                    gene_end_coordinates[ gene_id ] = (seqname,int(end))
            except KeyError:
                gene_end_coordinates[ gene_id ] = (seqname,int(end))

        self.gene_end_coordinates = gene_end_coordinates

    def add_transcript(self,_gene,_spot_coordinates,_transcript,_leftmost_transcript_position ):
        """
        Adds a transcript to the gene buffer
        Parameters:
        :param _gene: the name of the gene
        :param _spot_coordinates: the spot cordinates as a (x,y) tuple
        :param _transcript: the transcript information
            as a (chrom, start, end, clear_name, mapping_quality, strand, umi) tuple
        :param _leftmost_transcript_position: the leftmost coordinate of the  original pysam record
            (i.e. AlignedSegment.reference_start)
        """

        cdef str gene = _gene
        cdef tuple spot_coordinates = _spot_coordinates
        cdef tuple transcript = _transcript
        cdef int leftmost_transcript_position = _leftmost_transcript_position

        self.current_position = leftmost_transcript_position
        self.current_chromosome = transcript[0]

        cdef dict new_gene_dict
        cdef dict new_spot_dict
        cdef list new_reads_list
        cdef tuple gene_end_coordinate
        cdef list ambiguous_genes

        # First try to add the "transcript" to the existing buffer
        try:
            self.buffer[gene]['spots'][spot_coordinates].append(transcript)

        # In case the spot coordinates are not found make a new reads list,
        # add the "transcript" to it and add the spotcordinates to the specified gene
        except KeyError:
            new_reads_list = [transcript]
            try:
                self.buffer[gene]['spots'][spot_coordinates] = new_reads_list
            except KeyError:

                # In case the gene is not in the buffer make a new gene dictionary with spot coordinates
                # and read lists and add it to the gene buffer dictionary
                # first try to fetch the genomic right most coordinate of the gene
                # from the dictionary with information from the gtf file
                # (the rightmost coordinate is the "end" coordinate of the gene
                # when parsing the bamfile by genome coordinate)
                try:
                    gene_end_coordinate = self.gene_end_coordinates[gene]

                # if the gene annotation is not found in the gtf file
                # it can either be an amgiuous annotation from htseq or we should raise a error
                # try to get the end coordinate for HTSeq ambiguous annotatios
                except KeyError:
                    # trying to handle HTSeq ambigiuos gene annotations correctly,
                    # they will have an annotation with the following style __ambiguous[GENE1+GENE2]
                    if gene[0:len('__ambiguous[')] == '__ambiguous[':#'GENE1+GENE2]'
                        try: # to get the right most right most coordinate ;)
                            ambiguous_genes = gene[len('__ambiguous['):-1].split('+')
                            gene_end_coordinate = max(
                                [ self.gene_end_coordinates[amb_gene] for amb_gene in ambiguous_genes ]
                                )
                        except KeyError:
                            raise ValueError('ERROR:: gene with id {0} is not found in gtf file\n'.format(gene))
                    else:
                        raise ValueError('ERROR:: gene with id {0} is not found in gtf file\n'.format(gene))

                # Create the new gene entry and add it to the buffer
                new_spot_dict = {spot_coordinates:new_reads_list}
                new_gene_dict = {
                        'spots':new_spot_dict,
                        'gene_end_coordinate':gene_end_coordinate
                    }
                self.buffer[gene] = new_gene_dict

        # Lastly update the record counter and write info if verbose is true
        self.total_transcript_counter += 1
        if self.verbose and self.total_transcript_counter%100000==0:
            self.speed_of_last_100k_transcripts = round(100000.0/(time.time()-self.time_for_last_100k_print),2)
            self.time_for_last_100k_print = time.time()
            self.print_stats()

    def put_finished_genes_in_queue(self,queue, empty=False):
        """
        Function for sending genes from the buffer to a queue object.
        if the empty parameter is set to true no further reads are
        expected to be added to the buffer and all genes currently
        in the buffer will be put in the output queue disregarding the
        genomic coordinate of the last transcript added.
        :param queue: a queue object to which genes will be sent
        :param empty: flag for sending all genes presently in buffer (default False)
        """

        # check the gene buffer to see if we passed the end coordinate of any of the genes in the buffer
        # if that is the case no more reads will be added to the gene
        # and we can send the completed genes to the parent process by puting it in the multiprocessing.Queue

        cdef list _tmp
        _tmp = list(self.buffer.keys()) # temporary list used to avoid changeing the size of the dict during iteration
        cdef list deleted_genes = [] #list used to mark genes for removal from the dict once sent to the parent process

        # For each gene in the buffer
        if empty and self.verbose: self.print_stats()
        for gene in _tmp:
            
            # fix to include any "__no_feature" annotations
            if gene == '__no_feature' and not empty: continue

            # check if the current position is past the gene end coordinate
            if empty \
            or self.current_position > self.buffer[gene]['gene_end_coordinate'][1] \
            or self.current_chromosome != self.buffer[gene]['gene_end_coordinate'][0]:

                # print stats
                if self.verbose: self.print_stats( )
                if empty and self.verbose: sys.stdout.write('INFO:: yielding last gene(s) {0} ... \n'.format(gene))

                # HERE we send the gene back to the parent process
                # by putting the spot dictionary and gene name in the multiprocessing.Queue
                queue.put( (gene, self.buffer[gene]['spots']) )

                # Check that the gene has not beeen processed be fore just to be sure nothing funky is going on
                assert gene not in self.processed_genes, 'ERROR: the gene {0} cannot be present twice in the genome.\n'.format(gene)
                self.processed_genes[gene] = True

                # Mark the gene for deletion from the buffer
                deleted_genes.append(gene)

        # Free some memory by deleting the genes from the buffer that were just sent to the parent process
        for gene in deleted_genes:
            for read_list in self.buffer[gene]['spots'].values(): self.transcripts_sent_counter += len(read_list)
            self.buffer[gene]['spots'] = None
            self.buffer[gene] = {}
            del self.buffer[gene]
        if empty and self.verbose: self.print_stats()

    def print_stats(self, header=False):
        """
        Function that prints a "current status" row to stdout
        """

        if header:
            header_str = '\t'.join( [
                'SIZE (MB)',
                'TOTREADS',
                'READSINBUF',
                'GENESINBUF',
                'TIME (s)',
                'AV_SPEED (reads/s)',
                'CU_SPEED (reads/s)',
                'POSITION'
            ])
            header_str+='\n'
            sys.stdout.write(header_str)

        sys.stdout.write('\t'.join([
                str(round(asizeof(self)/(1000.0*1000.0),2)),
                str(self.total_transcript_counter),
                str(self.total_transcript_counter-self.transcripts_sent_counter),
                str(len(self.buffer)),
                str(round(time.time()-self.start_time,3)),
                str(round(self.total_transcript_counter/(time.time()-self.start_time),2)),
                str(self.speed_of_last_100k_transcripts),
                '{0}:{1}\n'.format(str(self.current_chromosome),str(self.current_position))
            ]))