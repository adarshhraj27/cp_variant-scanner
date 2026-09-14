#!/usr/bin/env python3
"""
main.py

CLI entry point for the CP Variant Scanner sprint tool.

Pipeline (mirrors the paper's methodology, simplified):
  1. Load a candidate CP gene panel                (gene_filter.load_gene_panel)
  2. Load one or more patient variant files         (vcf_parser.load_variants)
  3. Keep only variants in the gene panel / MITO     (gene_filter.filter_to_panel)
  4. Remove variants also seen in healthy controls   (control_filter.remove_control_variants)
  5. Classify remaining variants as P/LP, VUS,       (classifier.classify_all)
     or Likely Benign
  6. Report per-patient candidates + cohort-level     (report.py)
     diagnostic yield and a variant-breakdown chart

Usage:
    python main.py --patients data/sample_patient_variants.vcf \
                    --controls data/control_cohort_variants.vcf \
                    --panel data/cp_gene_panel.csv \
                    --out output/
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from vcf_parser import load_variants, load_control_keys, VariantFileError
from gene_filter import load_gene_panel, filter_to_panel
from control_filter import remove_control_variants
from classifier import classify_all, CLASS_P_LP
from report import (
    print_candidate_table,
    plot_variant_breakdown,
    compute_diagnostic_yield,
)


def _dedupe(variants: list) -> list:
    """Drop exact duplicate variant rows (same chrom/pos/ref/alt), keeping
    the first occurrence. Guards against double-counting if the same
    variant appears twice in one file (e.g. a re-run appended to a log, or
    a multi-sample VCF with a repeated line)."""
    seen = set()
    unique = []
    for v in variants:
        if v.key not in seen:
            seen.add(v.key)
            unique.append(v)
    return unique


def run_pipeline(patient_paths, control_paths, panel_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.isfile(panel_path):
        raise VariantFileError(f"Gene panel file not found: {panel_path}")
    panel = load_gene_panel(panel_path)
    if not panel:
        raise VariantFileError(f"Gene panel file is empty: {panel_path}")

    control_keys = set()
    for cpath in control_paths:
        control_keys |= load_control_keys(cpath)

    per_patient_classifications = {}
    all_analyzed_variants = []

    for ppath in patient_paths:
        patient_id = os.path.splitext(os.path.basename(ppath))[0]
        variants = _dedupe(load_variants(ppath))

        if not variants:
            print(f"\n=== Patient: {patient_id} ===")
            print("  No usable variant rows found in this file — skipping.")
            per_patient_classifications[patient_id] = []
            continue

        panel_variants = filter_to_panel(variants, panel)
        novel, in_controls = remove_control_variants(panel_variants, control_keys)
        classify_all(novel)

        print(f"\n=== Patient: {patient_id} ===")
        print(f"  Variants in gene panel / mito: {len(panel_variants)}")
        print(f"  Removed as seen-in-controls:   {len(in_controls)}")
        print(f"  Remaining for classification:  {len(novel)}\n")
        print_candidate_table(novel)

        per_patient_classifications[patient_id] = [v.classification for v in novel]
        all_analyzed_variants.extend(novel)

    yield_pct = compute_diagnostic_yield(per_patient_classifications)
    print(f"\n=== Cohort Summary ===")
    print(f"Patients analyzed: {len(patient_paths)}")
    print(f"Diagnostic yield (>=1 P/LP variant): {yield_pct:.1f}%")

    chart_path = os.path.join(out_dir, "variant_breakdown.png")
    plot_variant_breakdown(all_analyzed_variants, chart_path)
    print(f"Saved variant breakdown chart to: {chart_path}")


def main():
    parser = argparse.ArgumentParser(description="CP Variant Scanner - Sprint 1 tool")
    parser.add_argument(
        "--patients", nargs="+", default=["data/sample_patient_variants.vcf"],
        help="One or more patient variant files"
    )
    parser.add_argument(
        "--controls", nargs="+", default=["data/control_cohort_variants.vcf"],
        help="One or more control cohort variant files"
    )
    parser.add_argument(
        "--panel", default="data/cp_gene_panel.csv",
        help="CSV of candidate CP genes"
    )
    parser.add_argument(
        "--out", default="output/",
        help="Directory to write reports/charts to"
    )
    args = parser.parse_args()
    try:
        run_pipeline(args.patients, args.controls, args.panel, args.out)
    except VariantFileError as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
