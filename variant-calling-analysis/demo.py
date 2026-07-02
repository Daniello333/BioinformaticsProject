"""Demo script: parse a VCF, filter low-confidence calls, and summarize variants."""
from __future__ import annotations

import os

from varanalysis import classify_variant, filter_variants, parse_vcf, summarize

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "sample.vcf")


def main() -> None:
    variants = parse_vcf(DATA_PATH)
    print(f"Loaded {len(variants)} variants from {DATA_PATH}\n")

    print("=== Per-variant classification ===")
    for v in variants:
        vtype = classify_variant(v)
        print(f"{v.chrom}:{v.pos} {v.ref}>{','.join(v.alt)} qual={v.qual} depth={v.depth()} type={vtype}")
    print()

    print("=== Raw summary ===")
    for k, v in summarize(variants).items():
        print(f"{k}: {v}")
    print()

    filtered = filter_variants(variants, min_qual=30, min_depth=10, pass_only=True)
    print(f"=== After filtering (QUAL>=30, DP>=10, PASS only): {len(filtered)}/{len(variants)} kept ===")
    for k, v in summarize(filtered).items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
