"""
report.py

Summarizes classified variants and reproduces the kind of breakdown shown
in the paper's Figure (share of the cohort with P/LP variants, VUS, or
none, split by variant type). Also computes "diagnostic yield" the same
way the paper defines it: the percentage of patients with at least one
Pathogenic/Likely Pathogenic finding.
"""

from collections import Counter
import matplotlib.pyplot as plt

from classifier import CLASS_P_LP, CLASS_VUS, CLASS_BENIGN


def summarize_by_type(variants: list) -> Counter:
    """Count P/LP variants by variant type (SNV, INDEL, MITO, ...) — mirrors
    the paper's pie chart of which variant classes explained CP cases."""
    counts = Counter()
    for v in variants:
        if v.classification == CLASS_P_LP:
            counts[v.var_type] += 1
    return counts


def summarize_by_classification(variants: list) -> Counter:
    counts = Counter()
    for v in variants:
        counts[v.classification] += 1
    return counts


def print_candidate_table(variants: list) -> None:
    candidates = [v for v in variants if v.classification in (CLASS_P_LP, CLASS_VUS)]
    if not candidates:
        print("No candidate (P/LP or VUS) variants found.")
        return

    print(f"{'GENE':<10}{'CHROM':<7}{'POS':<12}{'TYPE':<8}{'CONSEQUENCE':<14}{'CLASS'}")
    print("-" * 70)
    for v in sorted(candidates, key=lambda x: x.classification):
        print(f"{v.gene:<10}{v.chrom:<7}{v.pos:<12}{v.var_type:<8}{v.consequence:<14}{v.classification}")


def plot_variant_breakdown(all_panel_variants: list, out_path: str) -> None:
    """Pie chart of classification breakdown across all panel/MITO variants
    considered — the same style of figure used in the presentation."""
    counts = summarize_by_classification(all_panel_variants)
    labels = list(counts.keys())
    sizes = list(counts.values())

    fig, ax = plt.subplots(figsize=(8, 6))
    wedges, _, autotexts = ax.pie(sizes, autopct="%1.1f%%", startangle=90)
    ax.set_title("Variant Classification Breakdown")
    ax.legend(
        wedges,
        [f"{l} ({c})" for l, c in zip(labels, sizes)],
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def compute_diagnostic_yield(per_patient_classifications: dict) -> float:
    """per_patient_classifications: {patient_id: [classification, ...]}
    Returns % of patients with at least one P/LP variant."""
    if not per_patient_classifications:
        return 0.0
    n_positive = sum(
        1 for classes in per_patient_classifications.values() if CLASS_P_LP in classes
    )
    return 100.0 * n_positive / len(per_patient_classifications)
