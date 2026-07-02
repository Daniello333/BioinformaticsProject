"""One-off script to (re)generate the committed sample dataset in data/.
Not needed to run the demo (data/ is already populated), but kept so the
dataset is reproducible from the simulator with a fixed seed.
"""
from __future__ import annotations

import os

from mlgenomics.simulate import simulate_genotypes

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def main() -> None:
    genotypes, labels, causal_snps = simulate_genotypes(
        n_samples=1000, n_snps=500, n_causal_snps=15, seed=42
    )
    os.makedirs(DATA_DIR, exist_ok=True)
    genotypes.to_csv(os.path.join(DATA_DIR, "genotypes.csv"))
    labels.to_csv(os.path.join(DATA_DIR, "labels.csv"))

    import pandas as pd

    pd.Series(causal_snps, name="causal_snp_id").to_csv(
        os.path.join(DATA_DIR, "causal_snps.csv"), index=False
    )
    print(
        f"Wrote {genotypes.shape[0]} samples x {genotypes.shape[1]} SNPs to {DATA_DIR} "
        f"({labels.sum()} cases / {(labels == 0).sum()} controls, {len(causal_snps)} causal SNPs)"
    )


if __name__ == "__main__":
    main()
