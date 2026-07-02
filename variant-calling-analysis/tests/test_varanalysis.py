import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varanalysis import classify_variant, filter_variants, is_transition, parse_vcf, summarize

DATA_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "sample.vcf")


def test_parse_vcf_count():
    variants = parse_vcf(DATA_PATH)
    assert len(variants) == 12


def test_parse_vcf_fields():
    variants = parse_vcf(DATA_PATH)
    v = variants[0]
    assert v.chrom == "chr1"
    assert v.pos == 14653
    assert v.ref == "A"
    assert v.alt == ["G"]
    assert v.qual == 87.5
    assert v.filter == "PASS"
    assert v.info["DP"] == "45"
    assert v.samples["SAMPLE1"]["GT"] == "0/1"


def test_depth_from_info():
    variants = parse_vcf(DATA_PATH)
    assert variants[0].depth() == 45


def test_classify_variant():
    variants = parse_vcf(DATA_PATH)
    by_pos = {(v.chrom, v.pos): v for v in variants}
    assert classify_variant(by_pos[("chr1", 14653)]) == "snp"
    assert classify_variant(by_pos[("chr2", 7712)]) == "insertion"
    assert classify_variant(by_pos[("chr2", 9945)]) == "deletion"
    assert classify_variant(by_pos[("chr3", 20115)]) == "deletion"


def test_is_transition():
    assert is_transition("A", "G") is True
    assert is_transition("C", "T") is True
    assert is_transition("A", "C") is False
    assert is_transition("AT", "G") is None


def test_filter_variants():
    variants = parse_vcf(DATA_PATH)
    filtered = filter_variants(variants, min_qual=30, min_depth=10, pass_only=True)
    assert all(v.qual >= 30 for v in filtered)
    assert all((v.depth() or 0) >= 10 for v in filtered)
    assert all(v.filter == "PASS" for v in filtered)
    assert len(filtered) < len(variants)


def test_summarize():
    variants = parse_vcf(DATA_PATH)
    summary = summarize(variants)
    assert summary["total_variants"] == 12
    assert summary["type_counts"]["snp"] > 0
    assert summary["ti_tv_ratio"] is not None
