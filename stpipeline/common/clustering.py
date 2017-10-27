""" 
This module contains some functions to cluster
molecular barcodes (UMIs) by hamming distance
"""

import numpy as np
from scipy.cluster.hierarchy import linkage,fcluster
from sklearn.cluster import AffinityPropagation
from collections import defaultdict
from stpipeline.common.cdistance import hamming_distance
import random
from collections import Counter

def countUMIHierarchical(molecular_barcodes, 
                         allowed_mismatches, 
                         method = "single"):
    """
    Tries to finds clusters of similar UMIs using a hierarchical clustering
    with a minimum distance (allowed_mismatches). 
    It returns a list with all the non clustered UMIs, for clusters of 
    multiple UMIs a random one will be selected.
    :param molecular_barcodes: a list of UMIs
    :param allowed_mismatches: how much distance we allow between clusters
    :param method: the type of distance algorithm when clustering 
                   (single more restrictive or complete less restrictive)
    :type allowed_mismatches: integer
    :type method: str 
    :return: a list of unique UMIs
    :rtype: list
    """
    # linkage will not work for distance matrices of 1x1 or 2x2 so for these rare cases
    # we use the naive clustering
    if len(molecular_barcodes) <= 2:
        return countUMINaive(molecular_barcodes, allowed_mismatches)
    # Distance computation function
    def d(coord):
        i,j = coord
        return hamming_distance(molecular_barcodes[i], molecular_barcodes[j])
    # Create hierarchical clustering and obtain flat clusters at the distance given
    indices = np.triu_indices(len(molecular_barcodes), 1)
    distance_matrix = np.apply_along_axis(d, 0, indices)
    linkage_cluster = linkage(distance_matrix, method=method)
    flat_clusters = fcluster(linkage_cluster, allowed_mismatches, criterion='distance')
    # Retrieve the unique clustered UMIs
    items = defaultdict(list)
    for i, item in enumerate(flat_clusters):
        items[item].append(i)
    return [molecular_barcodes[random.choice(members)] for members in items.itervalues()]
  
def countUMINaive(molecular_barcodes, allowed_mismatches):
    """
    Tries to finds clusters of similar UMIs using a naive proximity
    approach where UMIs are sorted and the ones that are consecutive
    and has hamming distance below the given number of miss-matches will
    be clustered together.
    It returns a list with all the non clustered UMIs, for clusters of 
    multiple UMIs a random one will be selected.
    :param molecular_barcodes: a list of UMIs
    :param allowed_mismatches: how much distance we allow between clusters
    :param method: the type of distance algorithm when clustering 
                   (single more restrictive or complete less restrictive)
    :type allowed_mismatches: integer
    :type method: str 
    :return: a list of unique UMIs
    :rtype: list
    """
    clusters_dict = {}
    nclusters = 0
    for i, molecular_barcode in enumerate(sorted(molecular_barcodes)):
        if i == 0:
            clusters_dict[nclusters] = [molecular_barcode]
        else:
            # compare distant of previous molecular barcodes and new one
            # if distance is between threshold we add it to the cluster 
            # otherwise we create a new cluster
            if hamming_distance(clusters_dict[nclusters][-1], molecular_barcode) <= allowed_mismatches:
                clusters_dict[nclusters].append(molecular_barcode)
            else:
                nclusters += 1
                clusters_dict[nclusters] = [molecular_barcode]
    # Return the non clustered UMIs
    return [random.choice(members) for members in clusters_dict.itervalues()]

def breadth_first_search(node, adj_list):
    """ This function has been obtained from 
    https://github.com/CGATOxford/UMI-tools
    The logic behind the algorithm to cluster UMIs using
    an adjacent distance matrix is described in 
    http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
    """
    searched = set()
    found = set()
    queue = set()
    queue.update((node,))
    found.update((node,))
    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.update((node,))
        queue.difference_update(searched)
    return found

def remove_umis(adj_list, cluster, nodes):
    """Removes the specified nodes from the cluster and returns
    the remaining nodes"""
    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node for x in nodes for node in adj_list[x]] + nodes)
    return cluster - nodes_to_remove
    
def dedup_adj(molecular_barcodes, allowed_mismatches):
    """ This function has been obtained from 
    https://github.com/CGATOxford/UMI-tools
    The logic behind the algorithm to cluster UMIs using
    an adjacent distance matrix is described in 
    http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
    """
    c = Counter(molecular_barcodes)
        
    def get_adj_list_adjacency(umis):
        return {umi: [umi2 for umi2 in umis if hamming_distance(umi, umi2) \
                      <= allowed_mismatches] for umi in umis}

    def get_connected_components_adjacency(graph, Counter):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: Counter[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components

    def get_best_adjacency(cluster, adj_list, counts):
        if len(cluster) == 1: return list(cluster)
        sorted_nodes = sorted(cluster, key=lambda x: counts[x], reverse=True)
        for i in range(len(sorted_nodes) - 1):
            if len(remove_umis(adj_list, cluster, sorted_nodes[:i+1])) == 0:
                return sorted_nodes[:i+1]

    def reduce_clusters_adjacency(adj_list, clusters, counts):
        # TS - the "adjacency" variant of this function requires an adjacency
        # list to identify the best umi, whereas the other variants don't
        # As temporary solution, pass adj_list to all variants
        unique_umis = list()
        for cluster in clusters:
            parent_umis = get_best_adjacency(cluster, adj_list, counts)
            unique_umis += parent_umis
        return unique_umis 

    adj_list = get_adj_list_adjacency(c.keys())
    clusters = get_connected_components_adjacency(adj_list, c)
    unique_umis = reduce_clusters_adjacency(adj_list, clusters, c)
    return unique_umis

def dedup_dir_adj(molecular_barcodes, allowed_mismatches):
    """ This function has been obtained from
    https://github.com/CGATOxford/UMI-tools
    The logic behind the algorithm to cluster UMIs using
    an adjacent distance matrix is described in
    http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract
    """
    c = Counter(molecular_barcodes)
    
    def get_adj_list_directional_adjacency(umis, counts):
        return {umi: [umi2 for umi2 in umis if hamming_distance(umi, umi2) <= allowed_mismatches and
                      counts[umi] >= (counts[umi2]*2)-1] for umi in umis}  
    
    def get_connected_components_adjacency(graph, Counter):
        found = list()
        components = list()
        for node in sorted(graph, key=lambda x: Counter[x], reverse=True):
            if node not in found:
                component = breadth_first_search(node, graph)
                found.extend(component)
                components.append(component)
        return components
       
    def reduce_clusters_directional_adjacency(adj_list, clusters, counts):
        return [cluster.pop() for cluster in clusters]
    
    adj_list = get_adj_list_directional_adjacency(c.keys(), c)
    clusters = get_connected_components_adjacency(adj_list, c)
    unique_umis = reduce_clusters_directional_adjacency(adj_list, clusters, c)
    return unique_umis

def affinity_umi_removal(molecular_barcodes_dict, _):
    """
    Tries to finds clusters of similar UMIs using an affinity based approach. 
    It returns a list with all the non clustered UMIs, for clusters of 
    multiple UMIs a random one will be selected.
    :param molecular_barcodes: a list of UMIs
    :return: a list of unique UMIs
    :rtype: list
    """
    molecular_barcodes = molecular_barcodes_dict.keys()
    if len(molecular_barcodes) <= 2:
        return molecular_barcodes
    words = np.asarray(molecular_barcodes)
    lev_similarity = -1 * np.array([[hamming_distance(w1,w2) for w1 in words] for w2 in words])
    affprop = AffinityPropagation(affinity="precomputed", damping=0.8)
    affprop.fit(lev_similarity)
    unique_clusters = list()
    for cluster_id in np.unique(affprop.labels_):
        exemplar = words[affprop.cluster_centers_indices_[cluster_id]]
        cluster = np.unique(words[np.nonzero(affprop.labels_==cluster_id)])
        unique_clusters.append(random.choice(cluster))
    return unique_clusters