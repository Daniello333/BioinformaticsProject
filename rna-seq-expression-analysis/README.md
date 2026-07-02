# RNA-seq Differential Expression Analysis

A differential expression pipeline for RNA-seq count data: CPM normalization,
per-gene Welch's t-test, Benjamini-Hochberg FDR correction, and volcano plot
visualization. Includes a negative-binomial count simulator so the project is
fully self-contained (no external dataset download required).

## Contents

- `rnaseq/simulate.py` — `simulate_counts()`: generates a synthetic gene x sample
  count matrix (negative-binomial model) with a known set of truly
  differentially expressed genes, for validation
- `rnaseq/diffexp.py` — `cpm_normalize()`, `benjamini_hochberg()`,
  `run_differential_expression()` (Welch's t-test per gene + FDR correction)
- `rnaseq/plotting.py` — `volcano_plot()`
- `data/counts.csv`, `data/metadata.csv` — a fixed, reproducible sample dataset
  (2000 genes x 12 samples, 6 control / 6 treatment, seed=42), plus
  `data/true_de_labels.csv` recording which genes were simulated as truly DE
- `generate_sample_data.py` — regenerates the files in `data/` from the simulator
- `demo.py` — loads `data/`, runs the pipeline, writes `output/diffexp_results.csv`
  and `output/volcano.png`
- `tests/test_rnaseq.py` — pytest unit tests, including a recovery test that
  checks the pipeline correctly flags most simulated DE genes as significant

## Usage

```bash
pip install -r requirements.txt
python demo.py
pytest tests/
```

```python
import pandas as pd
from rnaseq import run_differential_expression, volcano_plot

counts = pd.read_csv("data/counts.csv", index_col=0)
metadata = pd.read_csv("data/metadata.csv", index_col=0)
results = run_differential_expression(counts, metadata)
volcano_plot(results, "output/volcano.png")
```

## Notes

- Counts are simulated with a negative-binomial distribution, the standard
  approximation for RNA-seq overdispersion (mean-variance relationship).
- `run_differential_expression` compares `treatment` vs `control` (as labeled
  in the `condition` column of `metadata`), reporting `log2fc`, raw `pvalue`,
  and BH-adjusted `padj`.
- The `output/` directory is created by `demo.py`; results are not committed
  since they're regenerated deterministically from the seeded sample data.
