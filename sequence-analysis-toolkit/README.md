# Sequence Analysis Toolkit

A dependency-free toolkit implementing core DNA/protein sequence analysis
algorithms from scratch: FASTA I/O, composition stats, transcription/translation,
ORF finding, and pairwise alignment (Needleman-Wunsch & Smith-Waterman).

## Contents

- `seqtools/fasta.py` — FASTA reading/writing (`read_fasta`, `iter_fasta`, `write_fasta`)
- `seqtools/sequence.py` — GC content, reverse complement, transcription, translation
  (standard codon table), ORF finding on both strands/all frames, Hamming distance
- `seqtools/alignment.py` — Needleman-Wunsch (global) and Smith-Waterman (local)
  dynamic-programming alignment with traceback
- `data/sample.fasta` — three toy sequences for the demo
- `demo.py` — runs every feature against the sample data
- `tests/test_seqtools.py` — pytest unit tests

## Usage

```bash
pip install -r requirements.txt
python demo.py
pytest tests/
```

```python
from seqtools import read_fasta, gc_content, find_orfs, needleman_wunsch

records = read_fasta("data/sample.fasta")
seq = next(iter(records.values()))
print(gc_content(seq))
print(find_orfs(seq, min_protein_length=10))
```

## Notes

- The codon table implements the standard genetic code; `translate()` stops at
  the first stop codon by default (`to_stop=True`).
- `find_orfs` scans all 3 forward frames and all 3 reverse-complement frames,
  returning ORFs sorted by protein length (longest first).
- Alignment gap/match/mismatch scores are configurable keyword arguments.
