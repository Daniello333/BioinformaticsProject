"""Demo script: load aligned sequences, build a distance matrix, construct
UPGMA and Neighbor-Joining trees, and render dendrograms."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from phylo import distance_matrix, neighbor_joining, plot_dendrogram, upgma

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "sequences.fasta")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")


def read_fasta(path: str) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records[header] = "".join(chunks)
    return records


def main() -> None:
    sequences = read_fasta(DATA_PATH)
    print(f"Loaded {len(sequences)} aligned sequences from {DATA_PATH}\n")

    labels, matrix = distance_matrix(sequences, model="jukes-cantor")
    print("=== Jukes-Cantor distance matrix ===")
    header = "        " + " ".join(f"{l:>10s}" for l in labels)
    print(header)
    for i, label in enumerate(labels):
        row = " ".join(f"{matrix[i][j]:10.4f}" for j in range(len(labels)))
        print(f"{label:8s}{row}")
    print()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    upgma_root = upgma(labels, matrix)
    print("=== UPGMA tree (Newick) ===")
    print(upgma_root.newick())
    plot_dendrogram(upgma_root, os.path.join(OUTPUT_DIR, "upgma_tree.png"), title="UPGMA Tree")
    print(f"Dendrogram written to {OUTPUT_DIR}/upgma_tree.png\n")

    nj_root = neighbor_joining(labels, matrix)
    print("=== Neighbor-Joining tree (Newick) ===")
    print(nj_root.newick())
    plot_dendrogram(nj_root, os.path.join(OUTPUT_DIR, "nj_tree.png"), title="Neighbor-Joining Tree")
    print(f"Dendrogram written to {OUTPUT_DIR}/nj_tree.png")


if __name__ == "__main__":
    main()
