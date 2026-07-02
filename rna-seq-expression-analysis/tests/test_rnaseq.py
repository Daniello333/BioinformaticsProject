import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from rnaseq import benjamini_hochberg, cpm_normalize, run_differential_expression, simulate_counts


def test_cpm_normalize_sums_to_million():
    counts = pd.DataFrame({"s1": [10, 20, 70], "s2": [5, 15, 30]})
    cpm = cpm_normalize(counts)
    assert np.allclose(cpm.sum(axis=0), 1e6)


def test_benjamini_hochberg_monotonic_and_bounded():
    pvals = np.array([0.001, 0.02, 0.03, 0.5, 0.9, 0.0001])
    padj = benjamini_hochberg(pvals)
    assert (padj >= pvals - 1e-12).all()
    assert (padj <= 1.0).all()
    assert (padj >= 0.0).all()


def test_benjamini_hochberg_all_significant_when_all_tiny():
    pvals = np.array([1e-10, 1e-9, 1e-8])
    padj = benjamini_hochberg(pvals)
    assert (padj < 0.01).all()


def test_simulate_counts_shape():
    counts, metadata, is_de = simulate_counts(n_genes=100, n_de_genes=10, samples_per_group=3, seed=1)
    assert counts.shape == (100, 6)
    assert metadata.shape[0] == 6
    assert is_de.sum() == 10


def test_run_differential_expression_recovers_de_genes():
    counts, metadata, is_de = simulate_counts(
        n_genes=300, n_de_genes=30, samples_per_group=6, seed=7, fold_change_range=(3.0, 6.0)
    )
    results = run_differential_expression(counts, metadata)
    results = results.loc[counts.index]  # align back to original gene order

    sig_genes = results["padj"] < 0.05
    # Most genes flagged significant should be true positives given the large fold changes.
    true_positive_rate = (sig_genes & is_de).sum() / max(sig_genes.sum(), 1)
    assert true_positive_rate > 0.7
