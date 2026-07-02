"""Simulate a SNP genotype case/control dataset for polygenic risk prediction.

Each individual has a genotype (0/1/2 minor allele copies) at n_snps loci.
A small subset of "causal" SNPs contribute additively to disease log-odds
via a logistic model; the rest are pure noise, mimicking a simplified GWAS
polygenic risk scenario.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def simulate_genotypes(
    n_samples: int = 1000,
    n_snps: int = 500,
    n_causal_snps: int = 15,
    maf_range: tuple[float, float] = (0.05, 0.5),
    effect_size_range: tuple[float, float] = (0.5, 1.5),
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    """Return (genotypes, labels, causal_snp_ids).

    genotypes: DataFrame (n_samples, n_snps), values in {0, 1, 2}
    labels: Series of 0 (control) / 1 (case)
    causal_snp_ids: SNP column names that truly drive disease risk
    """
    rng = np.random.default_rng(seed)

    snp_ids = [f"rs{100000 + i}" for i in range(n_snps)]
    maf = rng.uniform(*maf_range, size=n_snps)

    # genotype counts ~ Binomial(2, minor allele frequency), the standard
    # Hardy-Weinberg approximation for unlinked SNPs.
    genotypes = rng.binomial(2, maf, size=(n_samples, n_snps))

    causal_idx = rng.choice(n_snps, size=n_causal_snps, replace=False)
    effect_signs = rng.choice([-1, 1], size=n_causal_snps)
    effect_sizes = effect_signs * rng.uniform(*effect_size_range, size=n_causal_snps)

    log_odds = np.full(n_samples, -1.0)  # intercept, keeps baseline prevalence moderate
    for k, idx in enumerate(causal_idx):
        log_odds += effect_sizes[k] * (genotypes[:, idx] - 2 * maf[idx])  # mean-centered

    prob_case = 1 / (1 + np.exp(-log_odds))
    labels = rng.binomial(1, prob_case)

    genotypes_df = pd.DataFrame(genotypes, columns=snp_ids, index=[f"sample_{i}" for i in range(n_samples)])
    labels_series = pd.Series(labels, index=genotypes_df.index, name="disease_status")
    causal_snp_ids = [snp_ids[i] for i in causal_idx]

    return genotypes_df, labels_series, causal_snp_ids
