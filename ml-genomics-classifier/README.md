# ML Genomics Classifier

A polygenic risk prediction demo: given SNP genotypes (0/1/2 minor allele
copies) across many loci, train classifiers to predict disease case/control
status and check whether the model recovers the handful of loci that
actually drive risk versus the hundreds that are noise — a simplified,
self-contained version of a GWAS-based risk model.

## Contents

- `mlgenomics/simulate.py` — `simulate_genotypes()`: simulates SNP genotypes
  under Hardy-Weinberg (Binomial(2, MAF)) and a logistic disease model driven
  by a small causal subset, returning the truly causal SNP IDs for validation
- `mlgenomics/model.py` — `train_logistic_regression()` (L1-penalized, for
  sparse signal recovery), `train_random_forest()`, `evaluate_model()`
  (accuracy/ROC-AUC/confusion matrix), `feature_importance()`
- `mlgenomics/plotting.py` — `plot_roc_curve()`, `plot_feature_importance()`
  (highlights truly causal SNPs in red)
- `data/genotypes.csv`, `data/labels.csv`, `data/causal_snps.csv` — a fixed,
  reproducible sample dataset (1000 samples x 500 SNPs, 15 causal, seed=42)
- `generate_sample_data.py` — regenerates the files in `data/`
- `demo.py` — trains both models, evaluates them, and reports how many of the
  top-10 most important features are truly causal SNPs
- `tests/test_mlgenomics.py` — pytest unit tests, including a signal-recovery test

## Usage

```bash
pip install -r requirements.txt
python demo.py
pytest tests/
```

```python
import pandas as pd
from sklearn.model_selection import train_test_split
from mlgenomics import train_logistic_regression, evaluate_model, feature_importance

genotypes = pd.read_csv("data/genotypes.csv", index_col=0)
labels = pd.read_csv("data/labels.csv", index_col=0).iloc[:, 0]
X_train, X_test, y_train, y_test = train_test_split(genotypes, labels, test_size=0.25, stratify=labels)

model = train_logistic_regression(X_train, y_train)
print(evaluate_model(model, X_test, y_test))
print(feature_importance(model, list(genotypes.columns)).head(10))
```

## Notes

- Both models use `class_weight="balanced"` since case/control cohorts are
  typically imbalanced; accuracy alone can be misleading here — ROC-AUC is
  the more reliable comparison metric, and the confusion matrix shows both
  error types.
- The L1 logistic regression is a reasonable default for this kind of sparse,
  high-dimensional (SNPs >> causal loci) recovery problem; the random forest
  is included as a nonlinear baseline for comparison.
- This is a simplified, single-population, unlinked-SNP simulation (no
  linkage disequilibrium, population structure, or epistasis) — real GWAS
  data requires substantially more careful handling of confounders.
