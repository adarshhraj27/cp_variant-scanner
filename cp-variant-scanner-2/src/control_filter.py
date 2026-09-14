"""
control_filter.py

Mirrors the paper's use of healthy control cohorts (CHILD, Inova) as a
baseline: any variant that also shows up in unaffected controls is treated
as common/benign background variation rather than a CP candidate, and is
removed from further consideration.
"""


def remove_control_variants(variants: list, control_keys: set) -> tuple[list, list]:
    """Split variants into (novel, seen_in_controls)."""
    novel, in_controls = [], []
    for v in variants:
        if v.key in control_keys:
            in_controls.append(v)
        else:
            novel.append(v)
    return novel, in_controls
