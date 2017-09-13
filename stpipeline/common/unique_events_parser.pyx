class UniqueEventsParser():
    """
    The UniqueEventsParser class is a more memory efficient implementations of the stpipeline.comon.sam_utils.parseUniqueEvents-function.
    The class takes two values as input filename and gff_filename, as well as two keyword arguments verbose and max_genes_in_memory.
    max_genes_in_memory defines the size of the underlying queue used to send genes from the file parsing function to the main process.
    The work done by the parser is defined by the self._worker_function function.
    In this function a annotated coordinate sorted bam file is parsed and genes are returned to the parent process through a multiprocessing.Queue self.q.
    
        Parses the transcripts present in the filename given as input.
        It expects a coordinate sorted BAM file where the spot coordinates, 
        gene and UMI are present as extra tags.
        And a gtf file defining the coordinates of the genes in the reference genome.
        Will yield put a dictionary per gene with a spot coordinate tuple as keys into the self.q multiprocessing.Queue
        foreach gene put dict: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ] 
    
    :param filename: the input file containing the annotated BAM records
    :param gff_filename: the gff file containing the gene coordinates
    :return: a instance of the UniqueEventsParser
    """
    
    def __init__(self, filename, gff_filename, verbose=False, max_genes_in_memory=100):
        self.filename = filename
        self.gff_filename = gff_filename
        self.verbose = verbose
        self.aborted = False
        self.max_genes_in_memory = max_genes_in_memory
    
    def stop(self,):
        """
        Function for stopping and aborting a running procsses.
        Mostly used for debugging and development.
        """

        if self.check_running() == 'RUNNING':
            
            import sys
            import os
            import signal

            if self.verbose: sys.stderr.write('KILLING Process='+str(self.p.pid)+'.\n')
            os.kill(self.p.pid,signal.SIGTERM)
            
            # empty the queue and join all processes!!
            self.q.close()
            while not self.q.empty(): self.q.get()
            self.q.join_thread()
            self.aborted = True
    
    def run(self, ):
        """
        Function that starts the subprocess reading the file(s)
        """
        
        # imports
        import multiprocessing
        import sys
        import time
        import re
        
        # create output queue
        self.q = multiprocessing.Queue(self.max_genes_in_memory)
        
        # start worker subprocess
        self.p = multiprocessing.Process(target=self._worker_function)
        self.aborted = False
        self.p.start()
        
        if self.verbose: sys.stderr.write('Process='+str(self.p.pid)+' starting to parse unique events.\n')
        
    def check_running(self, ):
        """
        Function that checks if the subprocess is running and return the status as as a string with value "RUNNING", "ABORTED" or "COMPLETE"
        """
        
        import sys
        
        # check for completion of worker process
        if self.p.is_alive(): return 'RUNNING'
        else:
            if self.aborted: return 'ABORTED'
            if self.p.exitcode == 0:
                self.p.join()
                if self.verbose: sys.stderr.write('Process='+str(self.p.pid)+' finished parsing unique events.\n')
                return 'COMPLETE'
            else:
                msg = 'ERROR:: Process='+str(self.p.pid)+' FAILED to parse unique events!!!.\n'
                if self.verbose: sys.stderr.write(msg)
                raise RuntimeError(msg)

    def read_gff_file(self, ):
        """
        function that reads the end coordinates of genes and returns them as values of a dictionary with the gene ID as key
        returns: dict: [GENE_id] => gene_end_coordinate
        """
        import sys

        cdef dict gene_end_coordinates
        cdef dict line_dict
        cdef str line
        cdef str key
        cdef str gene_id
        cdef list split_line
        cdef str seqname
        cdef int end
        cdef list attributes
        cdef str attribute
        
        if self.verbose: sys.stderr.write('INFO:: LOADING gtf file...\n')
        gene_end_coordinates = dict()
        
        # parse gtf file
        for line in open(self.gff_filename):
            
            # remove header
            if line[:4] == 'track' or line[0] == '#': continue
            
            # parse the line according to gtf format desc http://www.ensembl.org/info/website/upload/gff.html
            split_line = line.lstrip().rstrip().split('\t')
            line_dict = {}
            seqname = split_line[0]
            end = int(split_line[4])
            attributes = split_line[-1].rstrip().split(';')
            for attribute in attributes:
                if not attribute.lstrip().rstrip(): continue
                key = attribute.lstrip().split(' ')[0]
                value = attribute.lstrip().split(' ')[1]
                line_dict[key]=value 
            
            # save gene_id and rigtmost genomic coordinate of each gene to dictionary
            try:
                gene_id = line_dict['gene_id']
                if gene_id[0] == '"' and gene_id[-1] == '"': gene_id=gene_id[1:-1]
            except KeyError:
                raise ValueError('The gene_id attribute is missing in the annotation file ('+ self.gff_filename+')\nORIGINAL gff line: '+line+'\n')
            try:
                if int(end) > gene_end_coordinates[ gene_id ][1]:
                    gene_end_coordinates[ gene_id ] = (seqname,int(end))
            except KeyError:
                gene_end_coordinates[ gene_id ] = (seqname,int(end))
        
        # return the dict
        return gene_end_coordinates
    
    def print_stat_line(self, data, header=False):
        """
        Function that prints a current status row to stderr
        """
        
        import sys
        import time
        from pympler.asizeof import asizeof
        
        genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, chrom, start = data
        
        if header: sys.stderr.write('SIZE (MB)\tTOTREADS\tREADSINBUF\tGENESINBUF\tTIME (s)\tAV_SPEED (reads/s)\tCU_SPEED (reads/s)\tPOSITION\n')
        
        sys.stderr.write(
                str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
                str(tmp_counter_0)+'\t'+\
                str(tmp_counter_0-tmp_counter_1)+'\t'+\
                str(len(genes_buffer))+'\t'+\
                str(round(time.time()-start_time,3))+'\t'+\
                str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
                str(speed_last_100k)+'\t'+\
                str(chrom)+':'+str(start)+'\n'
            )

    def _worker_function(self, ):
        """
        Modified version of the stpipeline.comon.sam_utils.parseUniqueEvents-function:
        
        Parses the transcripts present in the filename given as input.
        It expects a coordinate sorted BAM file where the spot coordinates, 
        gene and UMI are present as extra tags
        Will yield a dictionary per gene with a spot coordinate tuple as keys
        foreach gene yield: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ] 
        :param filename: the input file containing the annotated BAM records
        :param gff_filename: the gff file containing the gene coordinates
        """

        #
        # Imports
        #
        from stpipeline.common.utils import fileOk
        from stpipeline.common.stats import qa_stats
        import os
        import logging 
        import pysam
        from collections import defaultdict
        import sys
        import time

        if self.verbose: sys.stderr.write('INFO:: ENTERING => parseUniqueEvents_byCoordinate\n')

        #
        # Get dict with end coordinates for all genes
        # to be able to tell when a gene is "done"
        # and thereby ready to be returned to the parent process
        #
        cdef dict gene_end_coordinates = {}
        gene_end_coordinates = self.read_gff_file()

        # Get conection to log and open the bamfile for reading
        logger = logging.getLogger("STPipeline")
        sam_file = pysam.AlignmentFile(self.filename, "rb")
   
        #
        # define variables cython style for more efficient looping and in memory storage
        #
        cdef dict genes_buffer = dict()
        cdef dict processed_genes = dict()
        cdef int tmp_counter_0 = 0
        cdef int tmp_counter_1 = 0
        cdef double speed_last_100k = 0
        cdef double start_time=time.time()
        cdef double time_last_100k = start_time
        cdef object rec
        cdef str clear_name
        cdef int mapping_quality
        cdef int start 
        cdef int end
        cdef str chrom 
        cdef str strand
        cdef tuple transcript
        cdef list deleted_genes
        cdef int x
        cdef int y
        cdef str gene
        cdef str seq
        cdef list _tmp
        cdef dict new_gene_dict
        cdef dict new_spot_dict
        cdef list new_reads_list
        cdef tuple spot_coordinates
        cdef tuple gene_end_coordinate
        
        if self.verbose: sys.stderr.write('INFO:: parsing bamfile ... \n')
        if self.verbose: self.print_stat_line( (genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, 'chrom', 'start'), header=True )

        # parse the bamfile record by record by  genome coordinate from first chromosome to last
        for rec in sam_file.fetch(until_eof=True):

            #
            # get the info about the record from bam
            #
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

            # Get TAGGD tags from bam
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
            
            # Skip reads that don't get a gene from HtSeq This should probably be handled in a nicer way later on!!
            # if gene[0:2] == '__': continue
            if gene == '__no_feature': continue
            
            # Create a new transcript and add it to the in memory gene_buffer dictionary
            transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
            spot_coordinates = (x,y)
            try:
                genes_buffer[gene]['spots'][spot_coordinates].append(transcript)
            except KeyError:
                new_reads_list = [transcript]
                try:
                    genes_buffer[gene]['spots'][spot_coordinates] = new_reads_list
                except KeyError:
                    
                    try:
                        gene_end_coordinate = gene_end_coordinates[gene]
                    except KeyError:
                        if gene[0:len('__ambiguous[')] == '__ambiguous[':#'ENSG00000127540+ENSG00000267059]'
                            try:
                                ambiguous_genes = gene[len('__ambiguous['):-1].split('+')
                                gene_end_coordinate = max( [ gene_end_coordinates[amb_gene] for amb_gene in ambiguous_genes ] )
                            except KeyError:
                                raise ValueError('ERROR:: gene with id '+gene+' is not found in gtf file\n'+'\n'.join(gene_end_coordinates.keys()))
                        else:
                            raise ValueError('ERROR:: gene with id '+gene+' is not found in gtf file\n'+'\n'.join(gene_end_coordinates.keys()))
                    
                    new_spot_dict = {spot_coordinates:new_reads_list}
                    new_gene_dict = {
                            'spots':new_spot_dict,
                            'gene_end_coordinate':gene_end_coordinate
                        }
                    genes_buffer[gene] = new_gene_dict
            
            #
            # check the gene_buffer to see if we passed the end coordinate of any of the genes
            # if that is the case no more reads will be added to the gene
            # and we can send the completed genes to the parent process
            #
            _tmp = list(genes_buffer.keys()) # to avoid changeing the size of the dict while iterating over it
            deleted_genes=[] # and to ble able to remove the genes from the dict once we sent them to the parent
            for gene in _tmp:
                if rec.reference_start > genes_buffer[gene]['gene_end_coordinate'][1] or chrom != genes_buffer[gene]['gene_end_coordinate'][0]: # true if the current position is past the gene end coordinate
                    if self.verbose: self.print_stat_line( (genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, chrom, start) ) # print stats
                    self.q.put( (gene, genes_buffer[gene]['spots']) ) # HERE we send the gene back to the parent process
                    assert gene not in processed_genes, 'ERROR: the gene '+gene+' cannot be present twice in the genome.\n'
                    processed_genes[gene] = True
                    deleted_genes.append(gene)
            for gene in deleted_genes: # free some memory
                for read_list in genes_buffer[gene]['spots'].values(): tmp_counter_1 += len(read_list)
                genes_buffer[gene]['spots'] = None
                genes_buffer[gene] = {}
                del genes_buffer[gene]
    
            #
            # Update counter and write info if verbose
            #
            tmp_counter_0 += 1
            if self.verbose and tmp_counter_0%100000==0:
                speed_last_100k = round(100000.0/(time.time()-time_last_100k),2)
                time_last_100k = time.time()
                self.print_stat_line( (genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, chrom, start) )
          
        #
        # Close the bam file and yield the last gene(s)
        #
        sam_file.close()
        if self.verbose: self.print_stat_line( (genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, chrom, start) )
        deleted_genes=[]
        for gene in genes_buffer.keys():
            if self.verbose: sys.stderr.write('INFO:: yielding last gene '+gene+' ... \n')
            self.q.put( (gene, genes_buffer[gene]['spots']) )
            assert gene not in processed_genes, 'ERROR: the gene '+gene+' cannot be present twice in the genome.\n'
            processed_genes[gene] = True
            deleted_genes.append(gene)
        for gene in deleted_genes:
            for read_list in genes_buffer[gene]['spots'].values():tmp_counter_1 += len(read_list)
            genes_buffer[gene]['spots'] = None
            genes_buffer[gene] = {}
            del genes_buffer[gene]
        if self.verbose: self.print_stat_line( (genes_buffer, tmp_counter_0, tmp_counter_1, start_time, speed_last_100k, chrom, start) )
        self.q.put( "COMPLETED" ) # send a last keyword to the parent process to define the end of the data

        #
        # cleanup and exit
        #
        import time
        while not self.q.empty():
            if self.verbose: sys.stderr.write('Process='+str(self.p.pid)+' Done processing input wainting for queue to empty.\n')
            time.sleep(1)
        self.q.close()
        self.q.join_thread()
        if self.verbose: sys.stderr.write('Process='+str(self.p.pid)+' returning.\n')
