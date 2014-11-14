#!/usr/bin/env python
""" This module contains some functions to deal with fastq files
"""

from stpipeline.common.utils import *
import logging 
from itertools import izip
import numpy as np
import scipy.cluster.hierarchy
import collections

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
    nbases = 0
    phred = 64 if qual64 else 33
    
    for qual in qscore[::-1]:
        if (ord(qual) - phred) < min_qual:
            nbases +=1
    
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
                     min_qual=20, min_length=28, qual64=False, outputFolder=None):
    """ 
    Converts reads in rw file appending the first (distance - trim)
    bases of fw and also add FW or RW string to reads names
    It also performs a bwa qualitry trim of the fw and rw reads, when
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
    else:
        logger.error("Error: Input format not recognized " + out_fw + " , " + out_rw)
        raise RuntimeError("Error: Input format not recognized " + out_fw + " , " + out_rw + "\n")

    logger.info("Start Reformatting and Filtering raw reads")
    
    out_fw_handle = safeOpenFile(out_fw, 'w')
    out_fw_writer = writefq(out_fw_handle)
    out_rw_handle = safeOpenFile(out_rw, 'w')
    out_rw_writer = writefq(out_rw_handle)

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

        if line2_trimmed is not None:
            # Add the barcode and polyTs from fw only if rw has not been completely trimmed
            new_seq = line1[1][:trim_fw] + line2_trimmed[1][trim_rw:]
            new_qual = line1[2][:trim_fw] + line2_trimmed[2][trim_rw:]
            record = (line2_trimmed[0], new_seq, new_qual)
            out_rw_writer.send(record)

        else:
            # write fake sequence so bowtie wont fail for having rw and fw with different lenghts
            out_rw_writer.send(getFake(line2))
            dropped_rw += 1  
    
    out_fw_writer.close()
    out_rw_writer.close()
    out_fw_handle.close()
    out_rw_handle.close()
    fw_file.close()
    rw_file.close()
    
    if not fileOk(out_fw) or not fileOk(out_rw):
        logger.error("Error: output file is not present " + out_fw + " , " + out_rw)
        raise RuntimeError("Error: output file is not present " + out_fw + " , " + out_rw + "\n")
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

def hamming_distance(s1, s2):
    """
    Returns the Hamming distance between equal-length sequences.
    """
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def extractMolecularBarcodes(reads, mc_start_position, mc_end_position):
    """ 
    Extracts a list of molecular barcodes from the list of reads given their
    start and end positions
    """
    molecular_barcodes = list()
    for read in reads:
        if mc_end_position > len(read):
            raise ValueError("Molecular barcode could not be found in the read " + read )
        molecular_barcodes.append(read[mc_start_position:mc_end_position])
    return molecular_barcodes

def computeDistanceMatrixFromSequences(reads):
    """
    Computes a distance matrix from a list of reads
    """
    n = len(reads)
    distance_matrix = np.zeros((n,n))

    for i, ele_1 in enumerate(reads):
        for j, ele_2 in enumerate(reads):
            if j >= i:
                break # Since the matrix is symmetrical we don't need to  calculate everything
            difference = hamming_distance(ele_1, ele_2)  
            distance_matrix[i, j] = difference
            distance_matrix[j, i] = difference
    return distance_matrix

def countMolecularBarcodesClustersHierarchical(reads, allowed_missmatches, mc_start_position, 
                                               mc_end_position, min_cluster_size):
    """
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_missmatches)
    It will return a list with the all the clusters and their elements
    It uses a hirarchical clustering approach to then get the flat clusters
    at a certain level (allowed_missmatches)
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    distance_matrix = computeDistanceMatrixFromSequences(molecular_barcodes)
    linkage = scipy.cluster.hierarchy.single(distance_matrix)
    #problem is that linkage will build the tree using relative distances (need to find a way to go from allowed_
    #missmatches to this relative values
    flat_clusters = scipy.cluster.hierarchy.fcluster(linkage, allowed_missmatches, criterion='distance')  
    
    clusters = []
    items = collections.defaultdict(list)
    for i, item in enumerate(flat_clusters):
        items[item].append(i)
    for item, members in items.iteritems():
        if len(members) >= min_cluster_size and len(members) > 1:
            clusters.append([molecular_barcodes[i] for i in members])

    return clusters

def countMolecularBarcodesClustersNaive(reads, allowed_missmatches, 
                                        mc_start_position, mc_end_position, min_cluster_size):
    """
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_missmatches)
    It will return a list with the all the clusters and their elements
    It uses a naive approach to iterate all reads and check for clusters
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    molecular_barcodes.sort()
    clusters_dict = {}
    nclusters = 0
    for i in xrange(0, len(molecular_barcodes)):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcodes[i]]
        else:
            last = clusters_dict[nclusters][-1]
            if hamming_distance(last, molecular_barcodes[i]) <= allowed_missmatches:
                clusters_dict[nclusters].append(molecular_barcodes[i])
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcodes[i]]
    clusters = []
    for item, members in clusters_dict.iteritems():
        if len(members) >= min_cluster_size and len(members) > 1:
            clusters.append(members)
            
    return clusters

def countMolecularBarcodesClustersNaiveFallBack(reads, allowed_missmatches, 
                                                mc_start_position, mc_end_position, min_cluster_size):
    """
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_missmatches)
    It will return a list with the all the clusters and their elements
    It uses a naive approach to iterate all reads and check for clusters
    """
    clusters = []
    centroids = []
    scores = []
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    molecular_barcodes.sort()
    
    for mc in molecular_barcodes:
        matched = False
        
        if len(clusters) == 0:
            clusters.append([mc])
            centroids.append([mc])
            scores.append([])
            continue

        for clustnum in xrange(len(clusters)):
            dist = hamming_distance(mc, centroids[clustnum][0])

            if dist <= allowed_missmatches:
                clusters[clustnum].append(mc)

                if len(scores[clustnum]) == 0:
                    scores[clustnum].append(dist)
                elif dist < scores[clustnum]:
                    scores[clustnum][0] = dist
                    centroids[clustnum][0] = mc

                matched = True
                break

        if not matched:       
            clusters.append([mc])
            centroids.append([mc])
            scores.append([])

    scores_filtered = []
    for cluster in clusters:
        if len(cluster) <= min_cluster_size and len(cluster) > 1:
            scores_filtered.append(cluster)     
    return scores_filtered
