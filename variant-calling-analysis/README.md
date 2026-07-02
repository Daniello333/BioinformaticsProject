# Variant Calling / VCF Analysis

A dependency-free VCF (Variant Call Format) parser and analysis toolkit:
parse variant calls, classify them (SNP/insertion/deletion/MNP), filter by
quality/depth/FILTER status, and compute summary statistics including the
transition/transversion (Ti/Tv) ratio.

## Contents

- `varanalysis/vcf_parser.py` — `Variant` dataclass + `parse_vcf()`, a minimal
  VCFv4.x parser (handles `#CHROM` header, `INFO`, `FORMAT`, and sample columns)
- `varanalysis/annotate.py` — `classify_variant()`, `is_transition()`,
  `filter_variants()`, `summarize()`
- `data/sample.vcf` — 12 synthetic variants across 4 chromosomes (SNPs, indels,
  a mix of PASS/low-quality calls)
- `demo.py` — parses, classifies, and filters the sample VCF, printing before/after summaries
- `tests/test_varanalysis.py` — pytest unit tests

## Usage

```bash
pip install -r requirements.txt
python demo.py
pytest tests/
```

```python
from varanalysis import parse_vcf, filter_variants, summarize

variants = parse_vcf("data/sample.vcf")
high_confidence = filter_variants(variants, min_qual=30, min_depth=10, pass_only=True)
print(summarize(high_confidence))
```

## Notes

- Classification uses the first ALT allele; multi-allelic records should be
  split beforehand if you need per-allele classification.
- `Variant.depth()` reads `INFO/DP`, falling back to the first sample's
  `FORMAT/DP` if `INFO/DP` is absent.
- Ti/Tv ratio (`ti_tv_ratio`) is a standard QC metric — genuine germline SNPs
  in humans typically show a ratio around 2.0-2.1 genome-wide (higher in exonic
  regions); a much lower ratio in real data often indicates false positives.
