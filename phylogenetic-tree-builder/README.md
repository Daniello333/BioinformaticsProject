# Phylogenetic Tree Builder

Builds phylogenetic trees from aligned sequences using distance-based methods
implemented from scratch: p-distance / Jukes-Cantor corrected distances, UPGMA
(clock-assumed hierarchical clustering), and Neighbor-Joining (no clock
assumption). Trees are output as Newick strings and rendered as dendrograms.

## Contents

- `phylo/distance.py` — `p_distance()`, `jukes_cantor_distance()`, `distance_matrix()`
- `phylo/tree.py` — `TreeNode` (with `.newick()` / `.leaf_names()`), `upgma()`,
  `neighbor_joining()`
- `phylo/plotting.py` — `plot_dendrogram()`
- `data/sequences.fasta` — 6 synthetic aligned sequences (300 bp) simulated
  from a common ancestor with a known clade structure (two 3-species clades
  plus a divergent outgroup) — a natural test case for tree-building
- `generate_sample_data.py` — regenerates `data/sequences.fasta`
- `demo.py` — computes the distance matrix, builds both trees, prints Newick,
  and writes `output/upgma_tree.png` / `output/nj_tree.png`
- `tests/test_phylo.py` — pytest unit tests, including a check that UPGMA
  correctly clusters the most similar pair first

## Usage

```bash
pip install -r requirements.txt
python demo.py
pytest tests/
```

```python
from phylo import distance_matrix, upgma

sequences = {"a": "ACGTACGT...", "b": "ACGTACGA...", "c": "TGCATGCA..."}
labels, matrix = distance_matrix(sequences, model="jukes-cantor")
tree = upgma(labels, matrix)
print(tree.newick())
```

## Notes

- Input sequences must already be aligned (equal length); this project does
  not perform multiple sequence alignment itself. Gap columns ('-') are
  excluded pairwise from p-distance.
- Jukes-Cantor correction accounts for unobserved multiple substitutions at a
  site and is always >= the raw p-distance; it saturates (returns a capped
  distance) when p >= 0.75.
- UPGMA assumes a constant molecular clock and produces an ultrametric,
  rooted tree. Neighbor-Joining does not assume a clock and is generally
  preferred for real data with variable substitution rates across lineages.
