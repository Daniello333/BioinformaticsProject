"""Minimal VCF (Variant Call Format) v4.x parser with no external dependencies."""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Variant:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: list[str]
    qual: float | None
    filter: str
    info: dict[str, str | bool]
    format_fields: list[str] = field(default_factory=list)
    samples: dict[str, dict[str, str]] = field(default_factory=dict)

    def depth(self) -> int | None:
        """Read depth from INFO DP, falling back to the first sample's DP."""
        if "DP" in self.info:
            try:
                return int(self.info["DP"])
            except (TypeError, ValueError):
                pass
        for sample_data in self.samples.values():
            if "DP" in sample_data:
                try:
                    return int(sample_data["DP"])
                except ValueError:
                    pass
        return None


def _parse_info(info_str: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    if info_str == ".":
        return info
    for field_str in info_str.split(";"):
        if "=" in field_str:
            key, value = field_str.split("=", 1)
            info[key] = value
        else:
            info[field_str] = True
    return info


def parse_vcf(path: str) -> list[Variant]:
    """Parse a VCF file into a list of Variant records. Skips header/meta lines."""
    variants: list[Variant] = []
    column_names: list[str] = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                column_names = line[1:].split("\t")
                continue

            fields_ = line.split("\t")
            chrom, pos, vid, ref, alt, qual, flt, info_str = fields_[:8]

            sample_names = column_names[9:] if len(column_names) > 9 else []
            format_fields: list[str] = []
            samples: dict[str, dict[str, str]] = {}
            if len(fields_) > 8:
                format_fields = fields_[8].split(":")
                for name, sample_str in zip(sample_names, fields_[9:]):
                    values = sample_str.split(":")
                    samples[name] = dict(zip(format_fields, values))

            variants.append(
                Variant(
                    chrom=chrom,
                    pos=int(pos),
                    id=vid,
                    ref=ref,
                    alt=alt.split(","),
                    qual=None if qual == "." else float(qual),
                    filter=flt,
                    info=_parse_info(info_str),
                    format_fields=format_fields,
                    samples=samples,
                )
            )
    return variants
