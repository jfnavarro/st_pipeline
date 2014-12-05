#!/usr/bin/env python
""" 
This module contains some functions to for clustering reads using
different approaches
"""

import numpy as np
import scipy.cluster.hierarchy
import collections
from stpipeline.common.distance import hamming_distance

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
