"""Core sequence operations: composition, transcription, translation, ORF finding."""
from __future__ import annotations

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def gc_content(seq: str) -> float:
    """Return the fraction of G/C bases in seq (denominator is full sequence length)."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for base in seq if base in "GC")
    return gc / len(seq)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def transcribe(seq: str) -> str:
    """DNA coding strand -> mRNA (T -> U)."""
    return seq.upper().replace("T", "U")


def translate(seq: str, to_stop: bool = True) -> str:
    """Translate a DNA (or RNA) sequence into a protein string using the standard codon table."""
    seq = seq.upper().replace("U", "T")
    protein = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def find_orfs(seq: str, min_protein_length: int = 30) -> list[dict]:
    """Find open reading frames on both strands across all three frames.

    Returns a list of dicts with keys: strand, frame, start, end, protein.
    Coordinates are 0-based, relative to the input (+) strand, half-open [start, end).
    """
    seq = seq.upper()
    orfs: list[dict] = []
    strands = [("+", seq), ("-", reverse_complement(seq))]

    for strand, strand_seq in strands:
        n = len(strand_seq)
        for frame in range(3):
            i = frame
            while i < n - 2:
                codon = strand_seq[i : i + 3]
                if codon == START_CODON:
                    j = i
                    protein_codons = []
                    found_stop = False
                    while j < n - 2:
                        c = strand_seq[j : j + 3]
                        if c in STOP_CODONS:
                            found_stop = True
                            break
                        protein_codons.append(CODON_TABLE.get(c, "X"))
                        j += 3
                    if found_stop and len(protein_codons) >= min_protein_length:
                        if strand == "+":
                            start, end = i, j + 3
                        else:
                            start, end = n - (j + 3), n - i
                        orfs.append(
                            {
                                "strand": strand,
                                "frame": frame,
                                "start": start,
                                "end": end,
                                "protein": "".join(protein_codons),
                            }
                        )
                    i = j + 3 if found_stop else i + 3
                else:
                    i += 3
    return sorted(orfs, key=lambda o: (-len(o["protein"]), o["start"]))


def hamming_distance(seq1: str, seq2: str) -> int:
    """Number of mismatched positions between two equal-length sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length for Hamming distance")
    return sum(a != b for a, b in zip(seq1, seq2))
