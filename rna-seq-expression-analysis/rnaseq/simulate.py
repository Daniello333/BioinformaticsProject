"""Simulate an RNA-seq raw count matrix with a known set of differentially
expressed genes, using a negative-binomial model (the standard approximation
for RNA-seq count overdispersion)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_counts(
    n_genes: int = 2000,
    n_de_genes: int = 200,
    samples_per_group: int = 4,
    mean_log_expression: float = 6.0,
    dispersion: float = 0.15,
    fold_change_range: tuple[float, float] = (1.5, 4.0),
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Return (counts, metadata, is_de) where:

    - counts: DataFrame of shape (n_genes, 2 * samples_per_group), integer counts
    - metadata: DataFrame with a 'condition' column ('control' or 'treatment')
    - is_de: boolean array, True for genes simulated as truly differentially expressed
    """
    rng = np.random.default_rng(seed)

    gene_ids = [f"GENE{i:05d}" for i in range(n_genes)]
    sample_ids = [f"control_{i+1}" for i in range(samples_per_group)] + [
        f"treatment_{i+1}" for i in range(samples_per_group)
    ]
    conditions = ["control"] * samples_per_group + ["treatment"] * samples_per_group

    base_mean = rng.lognormal(mean=mean_log_expression, sigma=1.2, size=n_genes)

    is_de = np.zeros(n_genes, dtype=bool)
    de_idx = rng.choice(n_genes, size=n_de_genes, replace=False)
    is_de[de_idx] = True

    directions = rng.choice([-1, 1], size=n_de_genes)
    fold_changes = rng.uniform(*fold_change_range, size=n_de_genes)
    log2fc_true = np.zeros(n_genes)
    log2fc_true[de_idx] = directions * np.log2(fold_changes)

    treatment_mean = base_mean * (2.0**log2fc_true)

    n_samples = len(sample_ids)
    counts = np.zeros((n_genes, n_samples), dtype=int)

    size_param = 1.0 / dispersion
    for g in range(n_genes):
        control_mu = max(base_mean[g], 1e-6)
        treat_mu = max(treatment_mean[g], 1e-6)
        p_control = size_param / (size_param + control_mu)
        p_treat = size_param / (size_param + treat_mu)
        counts[g, :samples_per_group] = rng.negative_binomial(size_param, p_control, size=samples_per_group)
        counts[g, samples_per_group:] = rng.negative_binomial(size_param, p_treat, size=samples_per_group)

    counts_df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    metadata_df = pd.DataFrame({"sample": sample_ids, "condition": conditions}).set_index("sample")

    return counts_df, metadata_df, is_de
