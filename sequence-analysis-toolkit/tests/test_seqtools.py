import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from seqtools import (
    find_orfs,
    gc_content,
    needleman_wunsch,
    reverse_complement,
    smith_waterman,
    transcribe,
    translate,
)
from seqtools.sequence import hamming_distance


def test_gc_content():
    assert gc_content("GCGC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("") == 0.0
    assert abs(gc_content("ATGC") - 0.5) < 1e-9


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"


def test_transcribe():
    assert transcribe("ATGC") == "AUGC"


def test_translate():
    assert translate("ATGGCCTAA") == "MA"
    assert translate("ATGTTTTGA") == "MF"
    assert translate("ATG") == "M"


def test_find_orfs_basic():
    # ATG + 3 codons + stop, well past the min length filter
    seq = "ATG" + "GCT" * 3 + "TAA"
    orfs = find_orfs(seq, min_protein_length=3)
    assert len(orfs) == 1
    # protein includes the initiating Met, as with real ORF translation
    assert orfs[0]["protein"] == "MAAA"
    assert orfs[0]["strand"] == "+"


def test_hamming_distance():
    assert hamming_distance("ATGC", "ATGC") == 0
    assert hamming_distance("ATGC", "ATGA") == 1


def test_needleman_wunsch_identical():
    result = needleman_wunsch("ACGT", "ACGT")
    assert result.seq1_aligned == "ACGT"
    assert result.seq2_aligned == "ACGT"
    assert result.score == 4


def test_needleman_wunsch_gap():
    result = needleman_wunsch("ACGT", "ACT")
    assert "-" in result.seq1_aligned or "-" in result.seq2_aligned
    assert len(result.seq1_aligned) == len(result.seq2_aligned)


def test_smith_waterman_local_match():
    # Shared local region "GATTACA" flanked by unrelated sequence
    seq1 = "TTTTGATTACAGGGG"
    seq2 = "CCCCGATTACATTTT"
    result = smith_waterman(seq1, seq2)
    assert "GATTACA" in result.seq1_aligned.replace("-", "")
