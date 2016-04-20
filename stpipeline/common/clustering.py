#!/usr/bin/env python
""" 
This module contains some functions to cluster
molecular barcodes (UMIs) by distance
"""

import numpy as np
from scipy.cluster.hierarchy import linkage,fcluster
from collections import defaultdict
from stpipeline.common.distance import hamming_distance
import random

def extractMolecularBarcodes(reads, mc_start_position, mc_end_position):
    """ 
    Extracts a list of molecular barcodes from the list of reads given their
    start and end positions and returns a list of 
    (molecular_barcode, read object, occurrences, number Ns)
    sorted by molecular_barcode, occurrences(reverse) and number of Ns
    :param reads: a list of tuples where the sequence of the read is at the first position
    :param mc_start_position: the start position of the molecular barcodes in the sequence
    :param mc_end_position: the end position of the molecular barcodes in the sequence
    :type mc_start_position: integer
    :type mc_end_position: integer
    :return: a list of tuples
    """
    assert(mc_end_position > mc_start_position and mc_start_position >= 0)
    # Create a list with the molecular barcodes and a hash with the occurrences
    molecular_barcodes = list()
    molecular_barcodes_counts = defaultdict(int)
    for read in reads:
        mc = read[0][mc_start_position:mc_end_position]
        molecular_barcodes.append((mc, read))
        molecular_barcodes_counts[mc] += 1
    # Creates a new list with the molecular barcodes, occurrences and N content
    molecular_barcodes_count_list = list()
    for mc in molecular_barcodes:
        count = molecular_barcodes_counts[mc[0]]
        count_n = mc[0].count("N")
        molecular_barcodes_count_list.append((mc[0],mc[1],count,count_n))
    # Return the list sorted by molecular barcode, occurrence(reverse) and number of Ns
    return sorted(molecular_barcodes_count_list, key=lambda x: (x[0], -x[2], x[3]))

def countMolecularBarcodesClustersHierarchical(molecular_barcodes, 
                                               allowed_mismatches,
                                               min_cluster_size, 
                                               method = "single"):
    """
    Tries to finds clusters of similar molecular barcodes given 
    a minimum cluster size and a minimum distance (allowed_mismatches). 
    It returns a list with the all the non clustered reads, for clusters of 
    multiple reads a random read will be chosen. 
    This will guarantee that the list of reads returned
    is unique and does not contain duplicates
    This approach finds clusters using hierarchical clustering
    :param molecular_barcodes: a list of tuples (molecular_barcode, read)
    :param allowed_mismatches: how much distance we allow between clusters
    :param min_cluster_size: minimum number of reads for a cluster to be
    :param method: the type of distance algorithm when clustering 
                   (single more restrictive or complete less restrictive)
    :type allowed_mismatches: integer
    :type min_cluser_size: integer
    :type method: str 
    :return: a list of unique reads
    """
    # linkage will not work for distance matrices of 1x1 or 2x2 so for these rare cases
    # we use the naive clustering
    if len(molecular_barcodes) <= 2:
        return countMolecularBarcodesClustersNaive(molecular_barcodes, 
                                                   allowed_mismatches,
                                                   min_cluster_size)
    # Distance computation function
    def d(coord):
        i,j = coord
        return hamming_distance(molecular_barcodes[i][0], molecular_barcodes[j][0])
    # Create hierarchical clustering and obtain flat clusters at the distance given
    indices = np.triu_indices(len(molecular_barcodes), 1)
    distance_matrix = np.apply_along_axis(d, 0, indices)
    linkage_cluster = linkage(distance_matrix, method=method)
    flat_clusters = fcluster(linkage_cluster, allowed_mismatches, criterion='distance')
    # Retrieve the original reads from the clusters found and filter them
    clusters = []
    items = defaultdict(list)
    for i, item in enumerate(flat_clusters):
        items[item].append(i)
    for item, members in items.iteritems():
        if len(members) >= min_cluster_size:
            # Cluster so we get a random read
            clusters.append(molecular_barcodes[random.choice(members)][1])
        else:
            # Single cluster so we add all the reads
            clusters += [molecular_barcodes[i][1] for i in members]
    return clusters
    
def countMolecularBarcodesClustersNaive(molecular_barcodes, 
                                        allowed_mismatches, 
                                        min_cluster_size):
    """
    Tries to finds clusters of similar molecular barcodes given 
    a minimum cluster size and a minimum distance (allowed_mismatches). 
    It returns a list with the all the non clustered reads, for clusters of 
    multiple reads a random read will be chosen. 
    This will guarantee that the list of reads returned
    is unique and does not contain duplicates
    This approach is a quick naive approach where the molecular barcodes are sorted
    and then added to clusters until the min distance is above the threshold
    :param molecular_barcodes: a list of tuples (molecular_barcode, read)
    :param allowed_mismatches: how much distance we allow between clusters
    :param min_cluster_size: minimum number of reads for a cluster to be
    :type allowed_mismatches: integer
    :type min_cluser_size: integer
    :return: a list of unique reads
    """
    clusters_dict = {}
    nclusters = 0
    for i, molecular_barcode in enumerate(sorted(molecular_barcodes, key=lambda x: (x[0]))):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcode]
        else:
            # compare distant of previous molecular barcodes and new one
            # if distance is between threshold we add it to the cluster 
            # otherwise we create a new cluster
            if hamming_distance(clusters_dict[nclusters][-1][0], 
                                molecular_barcode[0]) <= allowed_mismatches:
                clusters_dict[nclusters].append(molecular_barcode)
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcode]
    # Filter out computed clusters         
    clusters = []
    for _, members in clusters_dict.iteritems():
        # Check if the cluster is big enough
        # if so we pick a random read
        # otherwise we add all the reads
        if len(members) >= min_cluster_size:
            clusters.append(random.choice(members)[1])
        else:
            clusters += [member[1] for member in members]    
    return clusters
