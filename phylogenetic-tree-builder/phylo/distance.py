"""Pairwise genetic distance calculation from aligned sequences."""
from __future__ import annotations

import math


def p_distance(seq1: str, seq2: str) -> float:
    """Proportion of mismatched sites between two aligned (equal-length) sequences.
    Gap columns ('-') in either sequence are excluded from the comparison."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be aligned (equal length) to compute p-distance")

    compared = 0
    mismatches = 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a != b:
            mismatches += 1
    return mismatches / compared if compared else 0.0


def jukes_cantor_distance(seq1: str, seq2: str) -> float:
    """Jukes-Cantor corrected distance, accounting for multiple substitutions at a site.

    d = -3/4 * ln(1 - 4/3 * p), where p is the raw p-distance. Falls back to a large
    finite distance if p >= 0.75 (saturation, where the correction is undefined).
    """
    p = p_distance(seq1, seq2)
    if p >= 0.75:
        return 10.0  # saturated; cap rather than return infinity
    return -0.75 * math.log(1 - (4 / 3) * p)


def distance_matrix(
    sequences: dict[str, str], model: str = "jukes-cantor"
) -> tuple[list[str], list[list[float]]]:
    """Compute an all-pairs distance matrix. model: 'p-distance' or 'jukes-cantor'."""
    if model == "p-distance":
        dist_fn = p_distance
    elif model == "jukes-cantor":
        dist_fn = jukes_cantor_distance
    else:
        raise ValueError(f"Unknown distance model: {model}")

    labels = list(sequences.keys())
    n = len(labels)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = dist_fn(sequences[labels[i]], sequences[labels[j]])
            matrix[i][j] = d
            matrix[j][i] = d
    return labels, matrix
