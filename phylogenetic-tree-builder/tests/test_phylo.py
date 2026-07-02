import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from phylo import distance_matrix, jukes_cantor_distance, neighbor_joining, p_distance, upgma


def test_p_distance_identical():
    assert p_distance("ACGT", "ACGT") == 0.0


def test_p_distance_all_different():
    assert p_distance("ACGT", "TGCA") == 1.0


def test_p_distance_ignores_gaps():
    assert p_distance("AC-T", "AGGT") == 1 / 3


def test_jukes_cantor_greater_than_p_distance():
    seq1 = "ACGTACGTACGT"
    seq2 = "ACGAACGTTCGT"
    p = p_distance(seq1, seq2)
    jc = jukes_cantor_distance(seq1, seq2)
    assert jc >= p


def test_distance_matrix_symmetric():
    seqs = {"a": "ACGTACGT", "b": "ACGTACGA", "c": "TGCATGCA"}
    labels, matrix = distance_matrix(seqs, model="p-distance")
    n = len(labels)
    for i in range(n):
        assert matrix[i][i] == 0.0
        for j in range(n):
            assert abs(matrix[i][j] - matrix[j][i]) < 1e-9


def test_upgma_produces_all_leaves():
    seqs = {"a": "ACGTACGT", "b": "ACGTACGA", "c": "TGCATGCA", "d": "TGCATGCT"}
    labels, matrix = distance_matrix(seqs, model="p-distance")
    tree = upgma(labels, matrix)
    assert sorted(tree.leaf_names()) == sorted(labels)


def test_upgma_clusters_closer_pairs_first():
    # a and b are near-identical; c is far from both
    seqs = {"a": "AAAAAAAA", "b": "AAAAAAAT", "c": "TTTTTTTT"}
    labels, matrix = distance_matrix(seqs, model="p-distance")
    tree = upgma(labels, matrix)
    # the root's first child subtree should contain exactly {a, b}
    subtree_leaf_sets = [set(child.leaf_names()) for child, _ in tree.children]
    assert {"a", "b"} in subtree_leaf_sets


def test_neighbor_joining_produces_all_leaves():
    seqs = {"a": "ACGTACGT", "b": "ACGTACGA", "c": "TGCATGCA", "d": "TGCATGCT", "e": "ACGAACGT"}
    labels, matrix = distance_matrix(seqs, model="p-distance")
    tree = neighbor_joining(labels, matrix)
    assert sorted(tree.leaf_names()) == sorted(labels)


def test_newick_format_is_valid_ish():
    seqs = {"a": "ACGT", "b": "ACGA", "c": "TGCA"}
    labels, matrix = distance_matrix(seqs, model="p-distance")
    tree = upgma(labels, matrix)
    newick = tree.newick()
    assert newick.endswith(";")
    assert newick.count("(") == newick.count(")")
