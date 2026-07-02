"""Variant classification, filtering, and summary statistics."""
from __future__ import annotations

from collections import Counter

from .vcf_parser import Variant

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def classify_variant(variant: Variant) -> str:
    """Classify a variant as snp / insertion / deletion / mnp / complex.

    Uses the first ALT allele. Real VCFs can have multiple ALT alleles per
    record; callers that need per-allele classification should split first.
    """
    ref = variant.ref
    alt = variant.alt[0]

    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(ref) == len(alt):
        return "mnp"
    if len(ref) < len(alt) and alt.startswith(ref):
        return "insertion"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "deletion"
    return "complex"


def is_transition(ref: str, alt: str) -> bool | None:
    """True for transition (purine<->purine or pyrimidine<->pyrimidine), False for
    transversion, None if not a single-base substitution."""
    if len(ref) != 1 or len(alt) != 1:
        return None
    return (ref.upper(), alt.upper()) in TRANSITIONS


def filter_variants(
    variants: list[Variant],
    min_qual: float | None = None,
    min_depth: int | None = None,
    pass_only: bool = False,
) -> list[Variant]:
    """Filter variants by QUAL, depth, and FILTER status."""
    result = []
    for v in variants:
        if pass_only and v.filter not in ("PASS", "."):
            continue
        if min_qual is not None and (v.qual is None or v.qual < min_qual):
            continue
        if min_depth is not None:
            depth = v.depth()
            if depth is None or depth < min_depth:
                continue
        result.append(v)
    return result


def summarize(variants: list[Variant]) -> dict:
    """Produce summary statistics: counts by type, ti/tv ratio, per-chromosome counts."""
    type_counts = Counter(classify_variant(v) for v in variants)

    ti = tv = 0
    for v in variants:
        result = is_transition(v.ref, v.alt[0])
        if result is True:
            ti += 1
        elif result is False:
            tv += 1

    chrom_counts = Counter(v.chrom for v in variants)
    quals = [v.qual for v in variants if v.qual is not None]

    return {
        "total_variants": len(variants),
        "type_counts": dict(type_counts),
        "transitions": ti,
        "transversions": tv,
        "ti_tv_ratio": (ti / tv) if tv else None,
        "chrom_counts": dict(chrom_counts),
        "mean_qual": (sum(quals) / len(quals)) if quals else None,
    }
