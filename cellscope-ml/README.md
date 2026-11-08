# CellScope ML v1.0 — selected source code

**Status: selected source module, not the complete CellScope Studio application.**

This directory publishes a genuine research-code component from the supplied CellScope ML v1.0 archive:

- [`feature_selection.py`](feature_selection.py) — original fold-local variance filtering and ANOVA gene selection, designed to avoid feature-selection leakage in train/test evaluation.
- [`test_feature_selection.py`](test_feature_selection.py) — small additional smoke tests for this published subset.

Run locally:

```bash
pip install numpy scikit-learn pytest
cd cellscope-ml
python -m pytest -q test_feature_selection.py
```

The complete v1.0 project, containing Streamlit, command-line tools, classifiers, donor-aware evaluation, data validation and a larger test suite, **is not yet published in this directory**. Do not interpret this isolated module as a runnable full application or an externally validated classifier.

**Further reading:** [CellScope research case study](https://github.com/Daniello333/Daniello333/blob/main/projects/cellscope-ml.md).

**Author:** Daniel Drabik
