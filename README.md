# Bioinformatics Project

Five self-contained bioinformatics projects, each in its own folder with its
own `requirements.txt`, sample data, a `demo.py` you can run end-to-end, and
pytest unit tests. Every project generates or embeds its own sample data, so
nothing requires downloading external datasets.

## Projects

| Folder | What it does |
|---|---|
| [`sequence-analysis-toolkit/`](sequence-analysis-toolkit/) | FASTA I/O, GC content, transcription/translation, ORF finding, and pairwise alignment (Needleman-Wunsch & Smith-Waterman) implemented from scratch |
| [`variant-calling-analysis/`](variant-calling-analysis/) | VCF parsing, variant classification (SNP/insertion/deletion), quality filtering, and Ti/Tv summary statistics |
| [`rna-seq-expression-analysis/`](rna-seq-expression-analysis/) | RNA-seq differential expression pipeline: CPM normalization, per-gene Welch's t-test, Benjamini-Hochberg FDR correction, volcano plots |
| [`phylogenetic-tree-builder/`](phylogenetic-tree-builder/) | Distance-based phylogenetics: p-distance/Jukes-Cantor, UPGMA and Neighbor-Joining tree construction, Newick output, dendrogram rendering |
| [`ml-genomics-classifier/`](ml-genomics-classifier/) | Polygenic risk prediction from SNP genotypes using logistic regression and random forests, with causal-SNP recovery evaluation |

Each project's README has more detail on its contents, usage, and design notes.

## Running a project

```bash
cd <project-folder>
pip install -r requirements.txt
python demo.py
pytest tests/
```

## Requirements

Python 3.10+. Dependencies vary per project (see each `requirements.txt`);
the sequence and variant-calling toolkits are dependency-free aside from
`pytest`, while the RNA-seq, phylogenetics, and ML projects use
numpy/pandas/scipy/scikit-learn/matplotlib.
