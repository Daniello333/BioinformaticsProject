"""Minimal FASTA reader/writer with no external dependencies."""
from __future__ import annotations

from typing import Iterator


def read_fasta(path: str) -> dict[str, str]:
    """Parse a FASTA file into an ordered dict of {header: sequence}."""
    records: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            records[header] = "".join(chunks)
    return records


def iter_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Stream (header, sequence) pairs without loading the whole file."""
    header = None
    chunks: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks)


def write_fasta(path: str, records: dict[str, str], line_width: int = 70) -> None:
    """Write {header: sequence} records to a FASTA file, wrapping at line_width."""
    with open(path, "w") as handle:
        for header, seq in records.items():
            handle.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                handle.write(seq[i : i + line_width] + "\n")
