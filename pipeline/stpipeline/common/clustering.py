#!/usr/bin/env python
""" 
This module contains some functions to for clustering reads
by their molecular barcodes using different approaches
"""

import numpy as np
import scipy.cluster.hierarchy
import counttrie.counttrie as ct
from collections import defaultdict
from stpipeline.common.distance import hamming_distance
import random
from operator import itemgetter

def extractMolecularBarcodes(reads, mc_start_position, mc_end_position):
    """ 
    :param reads a list of tuples (read_name, sequence, quality)
    :param mc_start_position the start position of the molecular barcodes in the reads
    :param mc_end_position the end position of the molecular barcodes in the reads
    Extracts a list of molecular barcodes from the list of reads given their
    start and end positions and returns a list of 
    (molecular_barcode, (read_name, sequence, quality), ocurrences, number Ns)
    sorted by molecular_barcode, ocurrences(reverse) and number of Ns
    """
    assert(mc_end_position > mc_start_position and mc_start_position >= 0)
    # Create a list with the molecular barcodes and a hash with the occurrences
    molecular_barcodes = list()
    molecular_barcodes_counts = defaultdict(int)
    for read in reads:
        if mc_end_position > len(read[1]):
            raise ValueError("UMI could not be found in the read " + read + "\n")
        mc = read[1][mc_start_position:mc_end_position]
        molecular_barcodes.append((mc, read))
        molecular_barcodes_counts[mc] += 1
        
    # Creates a new list with the molecular barcodes, occurences and N content
    molecular_barcodes_count_list = list()
    for mc in molecular_barcodes:
        count = molecular_barcodes_counts[mc[0]]
        count_n = mc[0].count("N")
        molecular_barcodes_count_list.append((mc[0],mc[1],count,count_n))
    
    # Return the list sorted by molecular barcode, occurrence(reverse) and number of Ns
    return sorted(molecular_barcodes_count_list, key=lambda x: (x[0], -x[2], x[3]))

def countMolecularBarcodesClustersHierarchical(reads, allowed_mismatches, mc_start_position,
                                               mc_end_position, min_cluster_size, method = "single"):
    """
    :param reads the list of reads to be searched for clusters in the form of tuple (read_name, sequence, quality)
    :param allowed_mismatches how much distance we allow between clusters
    :param mc_start_position start position of the read part that we want to cluster
    :param mc_end_position end position of the read part that we want to cluster
    :param min_cluster_size min number of reads to be count as cluster
    :param method the type of distance algorithm when clustering (single more restrictive or complete less restrictive)
    This functions tries to finds clusters of similar reads given a min cluster size
    and a minimum distance (allowed_mismatches). The clusters are built using the molecular
    barcodes present in the reads sequences
    It will return a list with the all the reads, for clusters of reads a random
    read will be chosen. This will quarante that the list of reads returned
    is unique and does not contain biological duplicates
    This approach finds clusters using hierarchical clustering
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    def d(coord):
        i,j = coord
        return hamming_distance(molecular_barcodes[i][0], molecular_barcodes[j][0])
    indices = np.triu_indices(len(molecular_barcodes), 1)
    distance_matrix = np.apply_along_axis(d, 0, indices)
    linkage = scipy.cluster.hierarchy.linkage(distance_matrix, method=method)
    flat_clusters = scipy.cluster.hierarchy.fcluster(linkage, allowed_mismatches, criterion='distance')

    # Retrieve the original reads from the clusters found
    clusters = []
    items = defaultdict(list)
    for i, item in enumerate(flat_clusters):
        items[item].append(i)
    for item, members in items.iteritems():
        if len(members) >= min_cluster_size:
            # Cluster so we get a random read
            random_read = molecular_barcodes[random.choice(members)][1]
            clusters.append(random_read)
        else:
            # Single cluster so we add all the reads
            clusters += [molecular_barcodes[i][1] for i in members]

    return clusters

def countMolecularBarcodesPrefixtrie(reads, allowed_mismatches, mc_start_position, 
                                     mc_end_position, min_cluster_size, allow_indels=True):
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
    This approach builds clusters using a prefix trie
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    trie = ct.CountTrie()
    # keep a map of molecular_barcode -> reads to obtain the reads later in the clusters
    reads = defaultdict(list)
    
    # add the molecular barcodes to the prefix tree and also the reads to the map
    for molecular_barcode in molecular_barcodes:
        trie.add(molecular_barcode[0])
        reads[molecular_barcode[0]].append(molecular_barcode[1])
        
    # find clusters and add them to the final list of they are bigger than the minimum size
    clusters = []
    for molecular_barcode in molecular_barcodes:
        cluster = trie.find_equal_length_optimized(molecular_barcode[0], allowed_mismatches, allow_indels)
        size_cluster = 0
        original_reads = []
        for clustered_mc in cluster:
            size_cluster += trie.get_count(clustered_mc)
            # get the original reads
            original_reads += reads[clustered_mc]
            # remove the original reads to not add them twice
            del reads[clustered_mc]
            trie.remove(clustered_mc)
        
        # if it is a cluster we add 1 random read otherwise we add
        # all the reads
        if size_cluster >= min_cluster_size and len(original_reads) > 0:
            clusters.append(random.choice(original_reads))
        else:
            for read in original_reads:
                clusters.append(read)
            
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
    read will be chosen. This will quarantee that the list of reads returned
    is unique and does not contain biological duplicates
    This approach is a quick naive approach where the molecular barcodes are sorted
    and then added to clusters until the min distance is above the threshold
    """
    molecular_barcodes = extractMolecularBarcodes(reads, mc_start_position, mc_end_position)
    clusters_dict = {}
    nclusters = 0
    for i, molecular_barcode in enumerate(molecular_barcodes):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcode]
        else:
            # compare distant of previous molecular barcodes and new one
            # if distance is between threshold we add it to the cluster 
            # otherwise we create a new cluster
            if hamming_distance(clusters_dict[nclusters][-1][0], molecular_barcode[0]) <= allowed_mismatches:
                clusters_dict[nclusters].append(molecular_barcode)
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcode]
             
    clusters = []
    for item, members in clusters_dict.iteritems():
        # Check if the cluster is big enough
        # if so we pick a random read
        # otherwise we add all the reads
        if len(members) >= min_cluster_size:
            clusters.append(random.choice(members)[1])
        else:
            for member in members:
                clusters.append(member[1])
            
    return clusters
