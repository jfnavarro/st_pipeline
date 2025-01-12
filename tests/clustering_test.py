#! /usr/bin/env python
"""
Unit-test the package clustering
"""
from collections import Counter
from stpipeline.common.clustering import (
    _breadth_first_search,
    _remove_umis,
    _get_connected_components_adjacency,
    _get_adj_list_adjacency,
    _get_best_adjacency,
    _reduce_clusters_adjacency,
    _get_adj_list_directional_adjacency,
    _reduce_clusters_directional_adjacency,
    dedup_hierarchical,
    dedup_adj,
    dedup_dir_adj,
)


def test_breadth_first_search():
    adj_list = {
        "A": ["B", "C"],
        "B": ["A", "D"],
        "C": ["A"],
        "D": ["B"],
    }
    result = _breadth_first_search("A", adj_list)
    assert result == {"A", "B", "C", "D"}


def test_remove_umis():
    adj_list = {
        "A": ["B"],
        "B": ["A", "C"],
        "C": ["B"],
    }
    cluster = {"A", "B", "C"}
    nodes = ["B"]
    result = _remove_umis(adj_list, cluster, nodes)
    assert result == {"A"}


def test_get_connected_components_adjacency():
    adj_list = {
        "A": ["B"],
        "B": ["A", "C"],
        "C": ["B"],
        "D": [],
    }
    counts = Counter({"A": 3, "B": 2, "C": 1, "D": 4})
    result = _get_connected_components_adjacency(adj_list, counts)
    assert len(result) == 2
    assert {"A", "B", "C"} in result
    assert {"D"} in result


def test_get_adj_list_adjacency():
    umis = ["AAAA", "AAAT", "AATT", "TTTT"]
    allowed_mismatches = 1
    result = _get_adj_list_adjacency(umis, allowed_mismatches)
    assert "AAAA" in result and "AAAT" in result["AAAA"]
    assert "AATT" not in result["AAAA"]


def test_get_best_adjacency():
    adj_list = {
        "A": ["B"],
        "B": ["A", "C"],
        "C": ["B"],
    }
    cluster = ["A", "B", "C"]
    counts = Counter({"A": 3, "B": 2, "C": 1})
    result = _get_best_adjacency(cluster, adj_list, counts)
    assert result == ["A"]


def test_reduce_clusters_adjacency():
    adj_list = {
        "A": ["B"],
        "B": ["A", "C"],
        "C": ["B"],
    }
    clusters = [{"A", "B", "C"}]
    counts = Counter({"A": 3, "B": 2, "C": 1})
    result = _reduce_clusters_adjacency(adj_list, clusters, counts)
    assert result == ["A"]


def test_get_adj_list_directional_adjacency():
    umis = ["AAAA", "AAAT", "AATT", "TTTT"]
    counts = Counter({"AAAA": 4, "AAAT": 3, "AATT": 2, "TTTT": 1})
    allowed_mismatches = 1
    result = _get_adj_list_directional_adjacency(umis, counts, allowed_mismatches)
    assert "AAAA" in result and "AAAT" in result["AAAA"]
    assert "AATT" not in result["AAAA"]


def test_reduce_clusters_directional_adjacency():
    clusters = [{"A", "B", "C"}]
    result = _reduce_clusters_directional_adjacency(clusters)
    assert result == ["A"]


def test_dedup_hierarchical():
    umis = ["AAAA", "AAAT", "AATT", "TTTT"]
    allowed_mismatches = 1
    result = dedup_hierarchical(umis, allowed_mismatches)
    assert len(result) <= len(umis)


def test_dedup_adj():
    umis = ["AAAA", "AAAT", "AATT", "TTTT"]
    allowed_mismatches = 1
    result = dedup_adj(umis, allowed_mismatches)
    assert len(result) <= len(umis)


def test_dedup_dir_adj():
    umis = ["AAAA", "AAAT", "AATT", "TTTT"]
    allowed_mismatches = 1
    result = dedup_dir_adj(umis, allowed_mismatches)
    assert len(result) <= len(umis)
