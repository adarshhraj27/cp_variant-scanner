"""
gene_filter.py

Restricts a variant list to a candidate gene panel — the same idea as the
paper narrowing its search to genes repeatedly hit across the CP cohort
(e.g., COL4A1). In a real pipeline this step would use genomic coordinates
from a GTF/refGene file; here the panel CSV already carries the coordinate
range per gene so this also works as a positional sanity check.
"""

import csv


def load_gene_panel(path: str) -> dict:
    """Return {gene_symbol: {chrom, start, end, associated_condition}}."""
    panel = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            panel[row["gene"]] = {
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "associated_condition": row["associated_condition"],
            }
    return panel


def _gene_by_coordinates(chrom: str, pos: int, panel: dict):
    """Fallback lookup: if a variant's GENE tag is missing, blank, or wrong
    (e.g. an upstream annotation-tool error), still catch it if its position
    falls inside a panel gene's coordinate range. Returns the gene symbol or
    None."""
    for gene, info in panel.items():
        if info["chrom"] == chrom and info["start"] <= pos <= info["end"]:
            return gene
    return None


def filter_to_panel(variants: list, panel: dict) -> list:
    """Keep only variants that land in the candidate gene panel.

    Mitochondrial variants are always kept (MT genes aren't in the nuclear
    panel). For nuclear variants, matching happens two ways:
      1. By GENE tag, if it's already in the panel (the normal case).
      2. By genomic coordinate, as a fallback — this catches variants whose
         GENE annotation is missing/blank/mislabeled but whose position
         still falls inside a panel gene's range. Each kept variant gets a
         `.match_method` attribute ("name" or "coordinate") so downstream
         code (or a reviewer) can see which path flagged it.
    """
    kept = []
    for v in variants:
        if v.var_type == "MITO":
            v.match_method = "mito"
            kept.append(v)
        elif v.gene in panel:
            v.match_method = "name"
            kept.append(v)
        else:
            fallback_gene = _gene_by_coordinates(v.chrom, v.pos, panel)
            if fallback_gene:
                v.gene = fallback_gene
                v.match_method = "coordinate"
                kept.append(v)
    return kept
