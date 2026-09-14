"""
classifier.py

Stands in for the paper's pathogenicity call, which required agreement
across 5 of 7 prediction tools (SIFT, PolyPhen, CADD, etc.) plus manual
curation. We don't have access to those tools here, so this uses a
transparent rule-based heuristic on consequence + inheritance instead.
Swap `classify_variant` out for a real annotation-tool call (e.g., a
wrapper around ANNOVAR/CADD output) without touching the rest of the
pipeline — that's the point of keeping this in its own module.
"""

# Known pathogenic mitochondrial positions (illustrative — e.g. m.3243A>G
# is the classic MELAS mutation) used only to demo how prior knowledge
# can be layered onto the heuristic.
KNOWN_PATHOGENIC_MITO_POS = {3243}

CLASS_P_LP = "Pathogenic/Likely Pathogenic"
CLASS_VUS = "Variant of Uncertain Significance"
CLASS_BENIGN = "Likely Benign"


def classify_variant(v) -> str:
    if v.var_type == "MITO":
        return CLASS_P_LP if v.pos in KNOWN_PATHOGENIC_MITO_POS else CLASS_VUS

    if v.consequence in ("nonsense", "frameshift"):
        return CLASS_P_LP

    if v.consequence == "missense":
        return CLASS_P_LP if v.inheritance == "de_novo" else CLASS_VUS

    if v.consequence in ("synonymous", "intronic"):
        return CLASS_BENIGN

    return CLASS_VUS


def classify_all(variants: list) -> list:
    """Attach a `.classification` attribute to each variant and return the list."""
    for v in variants:
        v.classification = classify_variant(v)
    return variants
