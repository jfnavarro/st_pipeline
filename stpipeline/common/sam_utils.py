""" 
This module contains some functions and utilities for ST SAM/BAM files
"""

from stpipeline.common.utils import fileOk
from stpipeline.common.stats import qa_stats
import os
import logging 
import pysam
from collections import defaultdict

# TODO this function uses too much memory, optimize it. (Maybe Cython or C++)
def parseUniqueEvents(filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a BAM file where the spot coordinates, 
    gene and UMI are present as extra tags
    The output will be a dictionary 
    [spot][gene] -> (chrom, start, end, clear_name, mapping_quality, strand, umi). 
    :param filename: the input file containing the annotated BAM records
    :return: A dictionary of spots(x,y) to a map of gene names to a list of transcripts 
    (chrom, start, end, clear_name, mapping_quality, strand, umi)
    As map[(x,y)][gene]->list((chrom, start, end, clear_name, mapping_quality, strand, UMI))
    """
    
    logger = logging.getLogger("STPipeline")
    unique_events = defaultdict(lambda : defaultdict(list))
    sam_file = pysam.AlignmentFile(filename, "rb")
    for rec in sam_file.fetch(until_eof=True):
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
        # Get TAGGD tags
        x,y,gene,seq = (None,None,None,None)
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
        if None in [x,y,gene,umi]:
            logger.warning("Warning parsing annotated reads.\n" \
                           "Missing attributes for record {}\n".format(clear_name))
            continue
        
        # Create a new transcript and add it to the dictionary
        transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
        unique_events[(x,y)][gene].append(transcript)

    sam_file.close()
    return unique_events

def parseUniqueEvents_byCoordinate(filename, gff_filename):
    """
    Parses the transcripts present in the filename given as input.
    It expects a coordinate sorted BAM file where the spot coordinates, 
    gene and UMI are present as extra tags
    Will yield a dictionary per gene with a spot coordinate tuple as keys
    foreach gene yield: [spot] -> [(chrom, start, end, clear_name, mapping_quality, strand, umi), ... ] 
    :param filename: the input file containing the annotated BAM records
    :param gff_filename: the gff file containing the gene coordinates
    :return: A generator 
    """
    import sys
    from pympler.asizeof import asizeof
    import time
    sys.stderr.write('INFO:: ENTERING => parseUniqueEvents_byCoordinate\n')
    # gtf format desc http://www.ensembl.org/info/website/upload/gff.html
    sys.stderr.write('INFO:: LOADING gtf file...\n')
    gene_end_coordinates = dict()
    for line in open(gff_filename):
        if line[:4] == 'track' or line[0] == '#': continue
        seqname,source,feature,start,end,score,strand,frame,attributes = line.lstrip().rstrip().split('\t')
        line_dict = {}
        for attribute in attributes.rstrip().split(';'):
            if not attribute.lstrip().rstrip(): continue
            key = attribute.lstrip().split(' ')[0]
            value = attribute.lstrip().split(' ')[1]
            line_dict[key]=value 
        try:
            gene_id = line_dict['gene_id']
            if gene_id[0] == '"' and gene_id[-1] == '"': gene_id=gene_id[1:-1]
        except KeyError:
            raise ValueError('The gene_id attribute is missing in the annotation file ('+ gff_filename+')\nORIGINAL gff line: '+line+'\n')
        try:
            if int(end) > gene_end_coordinates[ gene_id ][1]:
                gene_end_coordinates[ gene_id ] = (seqname,int(end))
        except KeyError:
            gene_end_coordinates[ gene_id ] = (seqname,int(end))

    sys.stderr.write('INFO:: parsing bamfile ... \n')
    genes_buffer = dict()
    processed_genes = dict()
    logger = logging.getLogger("STPipeline")
    sam_file = pysam.AlignmentFile(filename, "rb")
    start_time=time.time()
    tmp_counter_0 = 0
    tmp_counter_1 = 0
    speed_last_100k = 0
    time_last_100k = start_time
    sys.stderr.write('SIZE (MB)\tTOTREADS\tREADSINBUF\tGENESINBUF\tTIME (s)\tAV_SPEED (reads/s)\tCU_SPEED (reads/s)\tPOSITION\n')
    sys.stderr.write(
        str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
        str(tmp_counter_0)+'\t'+\
        str(tmp_counter_0-tmp_counter_1)+'\t'+\
        str(len(genes_buffer))+'\t'+\
        str(round(time.time()-start_time,3))+'\t'+\
        str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
        str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')
    for rec in sam_file.fetch(until_eof=True):
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
        # Get TAGGD tags
        x,y,gene,seq = (None,None,None,None)
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
        if None in [x,y,gene,umi]:
            logger.warning("Warning parsing annotated reads.\n" \
                           "Missing attributes for record {}\n".format(clear_name))
            continue
        
        # Create a new transcript and add it to the dictionary
        transcript = (chrom, start, end, clear_name, mapping_quality, strand, umi)
        try:
            genes_buffer[gene]['spots'][(x,y)].append(transcript)
        except KeyError:
            try:
                genes_buffer[gene] = {
                    'spots':defaultdict(list),
                    'gene_end_coordinate':gene_end_coordinates[gene]
                }
            except KeyError:
                raise ValueError('ERROR:: gene with id '+gene+' is not found in gtf file\n'+'\n'.join(gene_end_coordinates.keys()))
            genes_buffer[gene]['spots'][(x,y)].append(transcript)
        
        _tmp = list(genes_buffer.keys())
        deleted_genes=[]
        for gene in _tmp:
            if rec.reference_start > genes_buffer[gene]['gene_end_coordinate'][1] or chrom != genes_buffer[gene]['gene_end_coordinate'][0]:
                sys.stderr.write(
#                    '                                                                                                  '+'\r'+\
                    str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
                    str(tmp_counter_0)+'\t'+\
                    str(tmp_counter_0-tmp_counter_1)+'\t'+\
                    str(len(genes_buffer))+'\t'+\
                    str(round(time.time()-start_time,3))+'\t'+\
                    str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
                    str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')
                #sys.stderr.write('INFO:: yielding gene '+gene+' ... \n')
                #sys.stderr.write('INFO:: refstart='+str(rec.reference_start)+' gene end='+str(genes_buffer[gene]['gene_end_coordinate'][1])+' chroms'+str(chrom)+' '+str(genes_buffer[gene]['gene_end_coordinate'][0])+' \n')
                yield gene, genes_buffer[gene]['spots']
                #sys.stderr.write('INFO:: back from gene '+gene+'\n')
                assert gene not in processed_genes, 'ERROR: the gene '+gene+' cannot be present twice in the genome.\n'
                processed_genes[gene] = True
                deleted_genes.append(gene)
        for gene in deleted_genes:
            for read_list in genes_buffer[gene]['spots'].values(): tmp_counter_1 += len(read_list)
            genes_buffer[gene]['spots'] = None
            genes_buffer[gene] = {}
            del genes_buffer[gene]
            # sys.stderr.write(
            #     str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
            #     str(tmp_counter_0)+'\t'+\
            #     str(tmp_counter_0-tmp_counter_1)+'\t'+\
            #     str(len(genes_buffer))+'\t'+\
            #     str(round(time.time()-start_time,3))+'\t'+\
            #     str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
            #     str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')

        tmp_counter_0 += 1
        if tmp_counter_0%100000==0:
                speed_last_100k = round(10000/(time.time()-time_last_100k),2)
                time_last_100k = time.time()
                sys.stderr.write(
#                    '                                                                                                  '+'\r'+\
                    str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
                    str(tmp_counter_0)+'\t'+\
                    str(tmp_counter_0-tmp_counter_1)+'\t'+\
                    str(len(genes_buffer))+'\t'+\
                    str(round(time.time()-start_time,3))+'\t'+\
                    str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
                    str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')
      
    sam_file.close()

    sys.stderr.write(
#        '                                                                                                  '+'\r'+\
        str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
        str(tmp_counter_0)+'\t'+\
        str(tmp_counter_0-tmp_counter_1)+'\t'+\
        str(len(genes_buffer))+'\t'+\
        str(round(time.time()-start_time,3))+'\t'+\
        str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
        str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')
    deleted_genes=[]
    for gene in genes_buffer.keys():
        sys.stderr.write('INFO:: yielding last gene '+gene+' ... \n')
        yield gene, genes_buffer[gene]['spots']
        assert gene not in processed_genes, 'ERROR: the gene '+gene+' cannot be present twice in the genome.\n'
        processed_genes[gene] = True
        deleted_genes.append(gene)
    for gene in deleted_genes:
        for read_list in genes_buffer[gene]['spots'].values():tmp_counter_1 += len(read_list)
        genes_buffer[gene]['spots'] = None
        genes_buffer[gene] = {}
        del genes_buffer[gene]
    sys.stderr.write(
#        '                                                                                                  '+'\r'+\
        str(round(asizeof(genes_buffer)/(1000.0*1000.0),2))+'\t'+\
        str(tmp_counter_0)+'\t'+\
        str(tmp_counter_0-tmp_counter_1)+'\t'+\
        str(len(genes_buffer))+'\t'+\
        str(round(time.time()-start_time,3))+'\t'+\
        str(round(tmp_counter_0/(time.time()-start_time),2))+'\t'+\
        str(speed_last_100k)+'\t'+str(chrom)+':'+str(start)+'\n')


def filterMappedReads(mapped_reads,
                      hash_reads,
                      file_output,
                      file_output_discarded=None):
    """ 
    Iterates a BAM file containing mapped reads 
    and discards reads that are not demultiplexed with TaggD
    (for that a dictionary with the read name as key and the X,Y and UMI)
    as values must be given.
    This function will add the X,Y coordinates and UMI as extra tags
    to the output BAM file. 
    It assumes all the reads are aligned (do not contain un-aligned reads),
    filtered for minimum read length and unique (no multimap).
    Demultiplexed reads with the extra tags (x,y and UMI) will be written
    to a file.
    :param mapped_reads: path to a BAM file containing the START alignments
    :param hash_reads: a dictionary of read_names to (x,y,umi) SAM tags
    :param file_output: the path to the file where to write the records
    :param file_output_discarded: the path to the file where to write discarded files
    :type mapped_reads: str
    :type hash_reads: dict
    :type file_output: str
    :type file_output_discarded: str
    :raises: RuntimeError
    """
    logger = logging.getLogger("STPipeline")
    
    if not os.path.isfile(mapped_reads):
        error = "Error, input file not present {}\n".format(mapped_reads)
        logger.error(error)
        raise RuntimeError(error)
    
    # Create output files handlers
    flag_read = "rb"
    flag_write = "wb"
    infile = pysam.AlignmentFile(mapped_reads, flag_read)
    outfile = pysam.AlignmentFile(file_output, flag_write, template=infile)
    if file_output_discarded is not None:
        outfile_discarded = pysam.AlignmentFile(file_output_discarded, 
                                                flag_write, template=infile)
    # Create some counters and loop the records
    dropped_barcode = 0
    present = 0
    for sam_record in infile.fetch(until_eof=True):
        present += 1
        discard_read = False
        # Add the UMI and X,Y coordinates as extra SAM tags
        try:
            # Using as key the read name as it was used to generate the dictionary
            # In order to save memory we truncate the read
            # name to only keep the unique part (lane, tile, x_pos, y_pos)
            # TODO this procedure is specific to only Illumina technology
            key = "".join(sam_record.query_name.split(":")[-4:])
            for tag in hash_reads[key]:
                # TODO add error check here
                tag_tokens = tag.split(":")
                sam_record.set_tag(tag_tokens[0], tag_tokens[2], tag_tokens[1])
            outfile.write(sam_record)
        except KeyError:
            dropped_barcode += 1
            if file_output_discarded is not None:
                outfile_discarded.write(sam_record)
    
    # Close handlers           
    infile.close()
    outfile.close()
    if file_output_discarded is not None:
        outfile_discarded.close()

    if not fileOk(file_output):
        error = "Error filtering mapped reads.\n" \
        "Output file is not present\n {}".format(file_output)
        logger.error(error)
        raise RuntimeError(error)
            
    logger.info("Finish processing aligned reads (R2):" \
                "\nPresent: {0}" \
                "\nDropped - barcode: {1}".format(present,dropped_barcode))
    
    # Update QA object
    qa_stats.reads_after_mapping = present
    qa_stats.reads_after_demultiplexing = (present - dropped_barcode)