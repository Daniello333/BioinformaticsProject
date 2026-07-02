"""One-off script to (re)generate data/sequences.fasta: a synthetic aligned
multiple sequence set with a known hierarchical (tree-like) structure, produced
by applying random point mutations along simulated lineages from a common
ancestor. Kept for reproducibility; data/ is already populated.
"""
from __future__ import annotations

import os
import random

BASES = "ACGT"
SEQ_LEN = 300


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if rng.random() < rate:
            seq_list[i] = rng.choice([b for b in BASES if b != seq_list[i]])
    return "".join(seq_list)


def main() -> None:
    rng = random.Random(7)
    ancestor = "".join(rng.choice(BASES) for _ in range(SEQ_LEN))

    # Simulate a rough tree: two deep clades, each with a couple of closely
    # related tips, plus one more distant outgroup.
    clade_a_root = mutate(ancestor, 0.03, rng)
    clade_b_root = mutate(ancestor, 0.03, rng)
    outgroup = mutate(ancestor, 0.15, rng)

    sequences = {
        "Species_A1": mutate(clade_a_root, 0.02, rng),
        "Species_A2": mutate(clade_a_root, 0.02, rng),
        "Species_A3": mutate(mutate(clade_a_root, 0.04, rng), 0.02, rng),
        "Species_B1": mutate(clade_b_root, 0.02, rng),
        "Species_B2": mutate(clade_b_root, 0.02, rng),
        "Outgroup": outgroup,
    }

    out_path = os.path.join(os.path.dirname(__file__), "data", "sequences.fasta")
    with open(out_path, "w") as handle:
        for name, seq in sequences.items():
            handle.write(f">{name}\n{seq}\n")
    print(f"Wrote {len(sequences)} aligned sequences ({SEQ_LEN} bp) to {out_path}")


if __name__ == "__main__":
    main()
