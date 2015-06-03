#!/usr/bin/env python
""" 
This module contains some functions to for clustering reads
by their molecular barcodes using different approaches
"""

import numpy as np
import scipy.cluster.hierarchy
from stpipeline.common.distance import hamming_distance
import collections
import random
from operator import itemgetter

def extractMolecularBarcodes(reads, mc_start_position, mc_end_position):
    """ 
    :param reads a list of tuples (read_name, sequence, quality)
    :param mc_start_position the start position of the molecular barcodes in the reads
    :param mc_end_position the end position of the molecular barcodes in the reads
    Extracts a list of molecular barcodes from the list of reads given their
    start and end positions and returns a list of (molecular_barcode, (read_name, sequence, quality))
    """
    molecular_barcodes = list()
    for read in reads:
        if mc_end_position > len(read[1]) or mc_end_position <= mc_start_position or mc_start_position < 0:
            raise ValueError("Molecular barcode could not be found in the read " + read)
        mc = read[1][mc_start_position:mc_end_position]
        molecular_barcodes.append((mc, read))
    return molecular_barcodes

def computeDistanceMatrixFromSequences(reads):
    """
    :param reads is a list of tuples (read_name, sequence, quality)
    Computes a distance matrix from a list of reads
    """
    n = len(reads)
    distance_matrix = np.zeros((n,n))

    for i, ele_1 in enumerate(reads):
        for j, ele_2 in enumerate(reads):
            if j >= i:
                break # Since the matrix is symmetrical we don't need to  calculate everything
            difference = hamming_distance(ele_1[1], ele_2[1])  
            distance_matrix[i, j] = difference
            distance_matrix[j, i] = difference
    return distance_matrix

def countMolecularBarcodesClustersHierarchical(reads, allowed_mismatches, mc_start_position,
                                               mc_end_position, min_cluster_size):
    """
    :param reads the list of reads to be searched for clusters in the form of tuple (read_name, sequence, quality)
    :param allowed_mismatches how much distance we allow between clusters
    :param mc_start_position start position of the read part that we want to cluster
    :param mc_end_position end position of the read part that we want to cluster
    :param min_cluster_size min number of reads to be count as cluster
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_mismatches). The clusters are built using the molecular
    barcodes present in the reads sequences
    It will return a list with the all the reads, for clusters of reads a random
    read will be chosen. This will quarante that the list of reads returned
    is unique and does not contain biological duplicates
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    distance_matrix = computeDistanceMatrixFromSequences(molecular_barcodes)
    linkage = scipy.cluster.hierarchy.single(distance_matrix)
    #problem is that linkage will build the tree using relative distances (need to find a way to go from allowed_
    #mismatches to this relative values
    flat_clusters = scipy.cluster.hierarchy.fcluster(linkage, allowed_mismatches, criterion='distance')
    
    clusters = []
    
    #TODO finish this
    
    #items = collections.defaultdict(list)
    #for i, item in enumerate(flat_clusters):
    #    items[item].append(i)
    #for item, members in items.iteritems():
    #    if len(members) >= min_cluster_size and len(members) > 1:
    #        clusters.append([molecular_barcodes[i] for i in members])

    return clusters

def countMolecularBarcodesClustersNaive(reads, allowed_mismatches,
                                        mc_start_position, mc_end_position, min_cluster_size):
    """
    :param reads the list of reads to be searched for clusters in the form of tuple (read_name, sequence, quality)
    :param allowed_mismatches how much distance we allow between clusters
    :param mc_start_position start position of the read part that we want to cluster
    :param mc_end_position end position of the read part that we want to cluster
    :param min_cluster_size min number of reads to be count as cluster
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_mismatches). The clusters are built using the molecular
    barcodes present in the reads sequences
    It will return a list with the all the reads, for clusters of reads a random
    read will be chosen. This will quarante that the list of reads returned
    is unique and does not contain biological duplicates
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    molecular_barcodes = sorted(molecular_barcodes, key=itemgetter(0))
    clusters_dict = {}
    nclusters = 0
    for i, molecular_barcode in enumerate(molecular_barcodes):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcode]
        else:
            #compare distant of previous molecular barcodes and new onec
            if hamming_distance(clusters_dict[nclusters][-1][0], molecular_barcode[0]) <= allowed_mismatches:
                clusters_dict[nclusters].append(molecular_barcode)
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcode]
             
    clusters = []
    for item, members in clusters_dict.iteritems():
        if len(members) >= min_cluster_size and len(members) > 1:
            clusters.append(random.choice(members)[1])
        else:
            for member in members:
                clusters.append(member[1])
            
    return clusters
