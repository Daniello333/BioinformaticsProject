"""One-off script to (re)generate the committed sample dataset in data/.

Not needed to run the demo (data/ is already populated), but kept so the
dataset is reproducible from the simulator with a fixed seed.
"""
from __future__ import annotations

import os

from rnaseq.simulate import simulate_counts

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def main() -> None:
    counts, metadata, is_de = simulate_counts(
        n_genes=2000,
        n_de_genes=200,
        samples_per_group=6,
        dispersion=0.05,
        fold_change_range=(2.0, 5.0),
        seed=42,
    )
    os.makedirs(DATA_DIR, exist_ok=True)
    counts.to_csv(os.path.join(DATA_DIR, "counts.csv"))
    metadata.to_csv(os.path.join(DATA_DIR, "metadata.csv"))

    import pandas as pd

    pd.Series(is_de, index=counts.index, name="is_true_de").to_csv(
        os.path.join(DATA_DIR, "true_de_labels.csv")
    )
    print(f"Wrote {counts.shape[0]} genes x {counts.shape[1]} samples to {DATA_DIR}")


if __name__ == "__main__":
    main()
