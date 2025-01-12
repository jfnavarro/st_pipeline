"""
This module contains some functions to cluster
molecular barcodes (UMIs) sequences by hamming distance
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster  # type: ignore
from collections import defaultdict
from stpipeline.common.cdistance import hamming_distance  # type: ignore
import random
from collections import Counter
from typing import List, Dict, Set, Any


def _breadth_first_search(node: str, adj_list: Dict[str, List[str]]) -> Set[str]:
    """
    Performs a breadth-first search (BFS) to find all connected components starting from a node.
    """
    searched = set()
    found = set(node)
    queue = set(node)
    while len(queue) > 0:
        node = (list(queue))[0]
        found.update(adj_list[node])
        queue.update(adj_list[node])
        searched.add(node)
        queue.difference_update(searched)
    return found


def _remove_umis(adj_list: Dict[str, List[str]], cluster: List[str], nodes: List[str]) -> Set[str]:
    """
    Removes the specified nodes from the cluster and returns
    the remaining nodes
    """
    nodes_to_remove = set([node for x in nodes for node in adj_list[x]] + nodes)
    return set(cluster) - nodes_to_remove


def _get_adj_list_adjacency(umis: List[str], allowed_mismatches: int) -> Dict[str, List[str]]:
    """
    Constructs an adjacency list where each UMI points to all other UMIs within
    the allowed mismatches.
    """
    return {
        umi: [
            umi2 for umi2 in umis if hamming_distance(umi.encode("UTF-8"), umi2.encode("UTF-8")) <= allowed_mismatches
        ]
        for umi in umis
    }


def _get_connected_components_adjacency(adj_list: Dict[str, List[str]], counts: Counter[str]) -> List[List[str]]:
    """
    Traverses the adjacency list to find all connected components.
    """
    found = set()
    components = []
    for node in sorted(adj_list, key=lambda x: counts[x], reverse=True):
        if node not in found:
            nodes = _breadth_first_search(node, adj_list)
            found.update(nodes)
            components.append(list(nodes))
    return components


def _get_best_adjacency(cluster: List[str], adj_list: Dict[str, List[str]], counts: Counter[str]) -> List[str]:
    """
    Identifies the best UMI or set of UMIs from a cluster based on adjacency and counts.
    """
    if len(cluster) == 1:
        return cluster
    sorted_nodes = sorted(cluster, key=lambda x: counts[x], reverse=True)
    for i in range(len(sorted_nodes) - 1):
        if len(_remove_umis(adj_list, cluster, sorted_nodes[: i + 1])) == 0:
            return sorted_nodes[: i + 1]
    return cluster


def _reduce_clusters_adjacency(
    adj_list: Dict[str, List[str]], clusters: List[List[str]], counts: Counter[str]
) -> List[str]:
    """
    Reduces clusters to their best representative UMIs.
    """
    unique_umis = []
    for cluster in clusters:
        unique_umis += _get_best_adjacency(cluster, adj_list, counts)
    return unique_umis


def _get_adj_list_directional_adjacency(
    umis: List[str], counts: Counter[str], allowed_mismatches: int
) -> Dict[str, List[str]]:
    """
    Constructs a directional adjacency list where each UMI points to all other UMIs within
    the allowed mismatches and satisfying the directional count condition.
    """
    return {
        umi: [
            umi2
            for umi2 in umis
            if hamming_distance(umi.encode("UTF-8"), umi2.encode("UTF-8")) <= allowed_mismatches
            and counts[umi] >= (counts[umi2] * 2) - 1
        ]
        for umi in umis
    }


def _reduce_clusters_directional_adjacency(clusters: List[List[str]]) -> List[str]:
    """
    Reduces clusters to their best representative UMIs by selecting one UMI per cluster.
    """
    return [cluster.pop() for cluster in clusters]


def dedup_hierarchical(molecular_barcodes: List[str], allowed_mismatches: int, method: str = "single") -> List[str]:
    """
    Deduplicates molecular barcodes using hierarchical clustering.

    Args:
        molecular_barcodes: A list of UMIs to cluster.
        allowed_mismatches: Maximum allowed mismatches (distance) between UMIs in a cluster.
        method: The linkage method for clustering, "single" for more restrictive or "complete"
                for less restrictive. Defaults to "single".

    Returns:
        A list of unique UMIs after deduplication.

    Raises:
        RuntimeError: If the input list is empty or another error occurs during clustering.
    """
    if len(molecular_barcodes) == 0:
        raise RuntimeError("The input UMIs cannot be empty")

    if len(molecular_barcodes) == 1:
        return molecular_barcodes

    if len(molecular_barcodes) == 2:
        return (
            molecular_barcodes
            if hamming_distance(molecular_barcodes[0].encode("UTF-8"), molecular_barcodes[1].encode("UTF-8"))
            <= allowed_mismatches
            else [random.choice(molecular_barcodes)]
        )

    def d(coord: Any) -> int:
        i, j = coord
        return hamming_distance(molecular_barcodes[i].encode("UTF-8"), molecular_barcodes[j].encode("UTF-8"))  # type: ignore

    indices = np.triu_indices(len(molecular_barcodes), 1)
    distance_matrix = np.apply_along_axis(d, 0, indices)
    linkage_cluster = linkage(distance_matrix, method=method)
    flat_clusters = fcluster(linkage_cluster, allowed_mismatches, criterion="distance")

    items = defaultdict(list)
    for i, item in enumerate(flat_clusters):
        items[item].append(i)

    return [molecular_barcodes[random.choice(members)] for members in list(items.values())]


def dedup_adj(molecular_barcodes: List[str], allowed_mismatches: int) -> List[str]:
    """
    Deduplicates molecular barcodes using an adjacency-based clustering algorithm.

    This function clusters similar UMIs based on an adjacency distance matrix, as described
    in the algorithm from http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract.

    Args:
        molecular_barcodes: A list of UMIs to deduplicate.
        allowed_mismatches: Maximum allowable Hamming distance between UMIs in a cluster.

    Returns:
        A list of unique UMIs after deduplication.

    Raises:
        RuntimeError: If the input list is empty or another error occurs during clustering.
    """
    if len(molecular_barcodes) == 0:
        raise RuntimeError("The input UMIs cannot be empty")
    c = Counter(molecular_barcodes)
    adj_list = _get_adj_list_adjacency(list(c.keys()), allowed_mismatches)
    clusters = _get_connected_components_adjacency(adj_list, c)
    unique_umis = _reduce_clusters_adjacency(adj_list, clusters, c)
    return unique_umis


def dedup_dir_adj(molecular_barcodes: List[str], allowed_mismatches: int) -> List[str]:
    """
    Deduplicates molecular barcodes using a directional adjacency-based clustering algorithm.

    This function clusters similar UMIs based on a directional adjacency distance matrix, as described
    in the algorithm from http://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract.

    Args:
        molecular_barcodes: A list of UMIs to deduplicate.
        allowed_mismatches: Maximum allowable Hamming distance between UMIs in a cluster.

    Returns:
        A list of unique UMIs after deduplication.

    Raises:
        RuntimeError: If the input list is empty or another error occurs during clustering.
    """
    if len(molecular_barcodes) == 0:
        raise RuntimeError("The input UMIs cannot be empty")
    c = Counter(molecular_barcodes)
    adj_list = _get_adj_list_directional_adjacency(list(c.keys()), c, allowed_mismatches)
    clusters = _get_connected_components_adjacency(adj_list, c)
    unique_umis = _reduce_clusters_directional_adjacency(clusters)
    return unique_umis
