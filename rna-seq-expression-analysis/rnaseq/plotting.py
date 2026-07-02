"""Plotting helpers for differential expression results."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def volcano_plot(
    results: pd.DataFrame,
    out_path: str,
    padj_threshold: float = 0.05,
    lfc_threshold: float = 1.0,
    title: str = "Volcano plot",
) -> None:
    """Save a volcano plot (log2FC vs -log10 p-value) highlighting significant genes."""
    neg_log_p = -np.log10(results["pvalue"].clip(lower=1e-300))
    is_sig = (results["padj"] < padj_threshold) & (results["log2fc"].abs() >= lfc_threshold)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        results.loc[~is_sig, "log2fc"],
        neg_log_p[~is_sig],
        s=8,
        c="lightgray",
        label="Not significant",
        alpha=0.7,
    )
    ax.scatter(
        results.loc[is_sig, "log2fc"],
        neg_log_p[is_sig],
        s=8,
        c="crimson",
        label=f"padj < {padj_threshold} & |log2FC| >= {lfc_threshold}",
        alpha=0.8,
    )
    ax.axvline(lfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change (treatment vs control)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
