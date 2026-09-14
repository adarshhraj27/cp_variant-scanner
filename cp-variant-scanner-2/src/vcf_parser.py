"""
vcf_parser.py

Parses the simplified tab-delimited VCF-like files used in this project.
Real VCF files follow the same column layout (CHROM, POS, ID, REF, ALT, INFO...),
just with far more metadata packed into INFO. Here we keep INFO as simple
key=value pairs (GENE=, TYPE=, CONSEQUENCE=, INHERITANCE=) to stand in for what
a real annotation tool like ANNOVAR would produce.
"""

import os
from dataclasses import dataclass, field


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = "NONE"
    var_type: str = "SNV"          # SNV / INDEL / CNV / MITO
    consequence: str = "unknown"    # missense / nonsense / frameshift / synonymous / intronic ...
    inheritance: str = "unknown"    # de_novo / inherited / maternal
    info: dict = field(default_factory=dict)

    @property
    def key(self):
        """Unique identifier used to match the same variant across files."""
        return (self.chrom, self.pos, self.ref, self.alt)


def _parse_info(info_str: str) -> dict:
    fields = {}
    for pair in info_str.strip().split(";"):
        if "=" in pair:
            k, v = pair.split("=", 1)
            fields[k] = v
    return fields


class VariantFileError(Exception):
    """Raised for problems reading a variant file, with enough context
    (file, line number, raw content) to fix the input quickly."""


def load_variants(path: str) -> list[Variant]:
    """Load a patient/cohort variant file into a list of Variant objects.

    Malformed lines are skipped with a warning rather than crashing the
    whole run — a single bad row in a multi-thousand-line VCF shouldn't
    take down the pipeline. A missing file raises a clear error instead of
    a raw FileNotFoundError traceback.
    """
    if not os.path.isfile(path):
        raise VariantFileError(f"Variant file not found: {path}")

    variants = []
    with open(path) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                print(f"[WARN] {path}:{line_num}: expected >=5 tab-separated "
                      f"columns, got {len(cols)} — skipping line: {line!r}")
                continue
            try:
                chrom, _id_or_pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
                pos = int(cols[1])
            except (ValueError, IndexError):
                print(f"[WARN] {path}:{line_num}: could not parse POS as an "
                      f"integer — skipping line: {line!r}")
                continue

            info = _parse_info(cols[5]) if len(cols) > 5 and cols[5] else {}
            variants.append(
                Variant(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    gene=info.get("GENE", "NONE") or "NONE",
                    var_type=info.get("TYPE", "SNV") or "SNV",
                    consequence=info.get("CONSEQUENCE", "unknown") or "unknown",
                    inheritance=info.get("INHERITANCE", "unknown") or "unknown",
                    info=info,
                )
            )
    return variants


def load_control_keys(path: str) -> set:
    """Load a control cohort file into a set of (chrom, pos, ref, alt) keys.
    Same skip-and-warn behavior as load_variants for malformed rows."""
    if not os.path.isfile(path):
        raise VariantFileError(f"Control file not found: {path}")

    keys = set()
    with open(path) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                print(f"[WARN] {path}:{line_num}: expected >=4 columns, "
                      f"got {len(cols)} — skipping line: {line!r}")
                continue
            try:
                chrom, pos, ref, alt = cols[0], int(cols[1]), cols[2], cols[3]
            except ValueError:
                print(f"[WARN] {path}:{line_num}: could not parse POS as an "
                      f"integer — skipping line: {line!r}")
                continue
            keys.add((chrom, pos, ref, alt))
    return keys
