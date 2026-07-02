"""Pairwise sequence alignment: Needleman-Wunsch (global) and Smith-Waterman (local)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class AlignmentResult:
    seq1_aligned: str
    seq2_aligned: str
    score: int

    def __str__(self) -> str:
        match_line = "".join(
            "|" if a == b and a != "-" else " " for a, b in zip(self.seq1_aligned, self.seq2_aligned)
        )
        return f"{self.seq1_aligned}\n{match_line}\n{self.seq2_aligned}\nScore: {self.score}"


def _score(a: str, b: str, match: int, mismatch: int) -> int:
    return match if a == b else mismatch


def needleman_wunsch(
    seq1: str, seq2: str, match: int = 1, mismatch: int = -1, gap: int = -2
) -> AlignmentResult:
    """Global alignment via the Needleman-Wunsch dynamic-programming algorithm."""
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap
            dp[i][j] = max(diag, up, left)

    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + _score(
            seq1[i - 1], seq2[j - 1], match, mismatch
        ):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return AlignmentResult("".join(reversed(aligned1)), "".join(reversed(aligned2)), dp[n][m])


def smith_waterman(
    seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -2
) -> AlignmentResult:
    """Local alignment via the Smith-Waterman dynamic-programming algorithm."""
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    best_score = 0
    best_pos = (0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = dp[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch)
            up = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap
            dp[i][j] = max(0, diag, up, left)
            if dp[i][j] > best_score:
                best_score = dp[i][j]
                best_pos = (i, j)

    aligned1, aligned2 = [], []
    i, j = best_pos
    while i > 0 and j > 0 and dp[i][j] > 0:
        if dp[i][j] == dp[i - 1][j - 1] + _score(seq1[i - 1], seq2[j - 1], match, mismatch):
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + gap:
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        else:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    return AlignmentResult("".join(reversed(aligned1)), "".join(reversed(aligned2)), best_score)
