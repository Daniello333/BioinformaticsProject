"""Demo script: load the sample count matrix, run differential expression,
save a volcano plot, and print the top hits."""
from __future__ import annotations

import os

import pandas as pd

from rnaseq import run_differential_expression, volcano_plot

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def main() -> None:
    counts = pd.read_csv(os.path.join(DATA_DIR, "counts.csv"), index_col=0)
    metadata = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"), index_col=0)

    print(f"Loaded {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(metadata)
    print()

    results = run_differential_expression(counts, metadata)

    n_sig = ((results["padj"] < 0.05) & (results["log2fc"].abs() >= 1.0)).sum()
    print(f"Significant genes (padj < 0.05, |log2FC| >= 1.0): {n_sig}\n")

    print("=== Top 10 differentially expressed genes ===")
    print(results.head(10).to_string(float_format=lambda x: f"{x:.4g}"))

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    results_path = os.path.join(OUTPUT_DIR, "diffexp_results.csv")
    results.to_csv(results_path)
    print(f"\nFull results written to {results_path}")

    plot_path = os.path.join(OUTPUT_DIR, "volcano.png")
    volcano_plot(results, plot_path, title="Treatment vs Control")
    print(f"Volcano plot written to {plot_path}")


if __name__ == "__main__":
    main()
