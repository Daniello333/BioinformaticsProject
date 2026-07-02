"""Plotting helpers for classifier evaluation."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def plot_roc_curve(results_by_model: dict[str, dict], out_path: str, title: str = "ROC Curve") -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    for name, results in results_by_model.items():
        fpr, tpr = results["fpr_tpr"]
        ax.plot(fpr, tpr, label=f"{name} (AUC={results['roc_auc']:.3f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="gray", linewidth=1, label="Chance")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_feature_importance(
    importance: pd.Series, out_path: str, top_n: int = 20, causal_snps: set[str] | None = None
) -> None:
    top = importance.head(top_n)
    colors = ["crimson" if (causal_snps and name in causal_snps) else "steelblue" for name in top.index]

    fig, ax = plt.subplots(figsize=(7, max(4, top_n * 0.3)))
    ax.barh(top.index[::-1], top.values[::-1], color=colors[::-1])
    ax.set_xlabel("Importance")
    ax.set_title(f"Top {top_n} SNPs by importance" + (" (red = truly causal)" if causal_snps else ""))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
