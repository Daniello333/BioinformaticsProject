"""Differential expression analysis: CPM normalization, per-gene Welch's t-test,
and Benjamini-Hochberg FDR correction."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def cpm_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts-per-million normalization: scales each sample by its library size."""
    library_sizes = counts.sum(axis=0)
    return counts.div(library_sizes, axis=1) * 1e6


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values (q-values)."""
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]

    adjusted = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity from the largest p-value downward
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result = np.empty(n)
    result[order] = adjusted
    return result


def run_differential_expression(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str = "condition",
    control: str = "control",
    treatment: str = "treatment",
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    """Run per-gene differential expression testing.

    Returns a DataFrame indexed by gene with columns:
    mean_control, mean_treatment, log2fc, pvalue, padj.
    """
    cpm = cpm_normalize(counts)

    control_samples = metadata.index[metadata[group_col] == control]
    treatment_samples = metadata.index[metadata[group_col] == treatment]

    control_vals = cpm[control_samples].to_numpy()
    treatment_vals = cpm[treatment_samples].to_numpy()

    mean_control = control_vals.mean(axis=1)
    mean_treatment = treatment_vals.mean(axis=1)
    log2fc = np.log2(mean_treatment + pseudocount) - np.log2(mean_control + pseudocount)

    _, pvalues = stats.ttest_ind(treatment_vals, control_vals, axis=1, equal_var=False)
    pvalues = np.nan_to_num(pvalues, nan=1.0)
    padj = benjamini_hochberg(pvalues)

    return pd.DataFrame(
        {
            "mean_control": mean_control,
            "mean_treatment": mean_treatment,
            "log2fc": log2fc,
            "pvalue": pvalues,
            "padj": padj,
        },
        index=counts.index,
    ).sort_values("padj")
