"""Demo script exercising the seqtools package end to end."""
from __future__ import annotations

import os

from seqtools import (
    find_orfs,
    gc_content,
    needleman_wunsch,
    read_fasta,
    reverse_complement,
    smith_waterman,
    transcribe,
    translate,
)

DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "sample.fasta")


def main() -> None:
    records = read_fasta(DATA_PATH)
    print(f"Loaded {len(records)} sequences from {DATA_PATH}\n")

    for header, seq in records.items():
        print(f"--- {header} ---")
        print(f"Length: {len(seq)} bp")
        print(f"GC content: {gc_content(seq):.2%}")
        print(f"Reverse complement: {reverse_complement(seq)[:40]}...")
        print(f"mRNA (first 20 nt): {transcribe(seq)[:20]}")
        print(f"Protein (frame 0): {translate(seq)}")
        print()

    seq_list = list(records.items())
    print("=== ORF scan on seq1 ===")
    header, seq = seq_list[0]
    orfs = find_orfs(seq, min_protein_length=3)
    for orf in orfs[:5]:
        print(
            f"strand={orf['strand']} frame={orf['frame']} "
            f"[{orf['start']}:{orf['end']}] protein={orf['protein']}"
        )
    print()

    print("=== Global alignment (Needleman-Wunsch): seq1 vs seq2 ===")
    (h1, s1), (h2, s2) = seq_list[0], seq_list[1]
    result = needleman_wunsch(s1, s2)
    print(result)
    print()

    print("=== Local alignment (Smith-Waterman): seq1 vs seq3 ===")
    (h1, s1), (h3, s3) = seq_list[0], seq_list[2]
    result = smith_waterman(s1, s3)
    print(result)


if __name__ == "__main__":
    main()
