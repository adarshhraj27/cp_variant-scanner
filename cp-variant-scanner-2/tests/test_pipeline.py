"""
Basic unit tests for the CP Variant Scanner pipeline.
Run with:  python -m pytest tests/  (from the project root)
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from vcf_parser import load_variants, load_control_keys, Variant, VariantFileError
from gene_filter import load_gene_panel, filter_to_panel
from control_filter import remove_control_variants
from classifier import classify_variant, classify_all, CLASS_P_LP, CLASS_VUS, CLASS_BENIGN

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def test_load_variants():
    variants = load_variants(os.path.join(DATA_DIR, "sample_patient_variants.vcf"))
    assert len(variants) == 12
    assert variants[0].gene == "COL4A1"
    assert variants[0].chrom == "13"


def test_load_control_keys():
    keys = load_control_keys(os.path.join(DATA_DIR, "control_cohort_variants.vcf"))
    assert ("13", 110890456, "C", "T") in keys
    assert len(keys) == 5


def test_gene_panel_filter_keeps_panel_genes_and_mito():
    panel = load_gene_panel(os.path.join(DATA_DIR, "cp_gene_panel.csv"))
    variants = load_variants(os.path.join(DATA_DIR, "sample_patient_variants.vcf"))
    filtered = filter_to_panel(variants, panel)
    genes_kept = {v.gene for v in filtered}
    assert "NONE" not in genes_kept          # off-panel genes dropped
    assert any(v.var_type == "MITO" for v in filtered)  # mito always kept


def test_control_filter_splits_correctly():
    panel = load_gene_panel(os.path.join(DATA_DIR, "cp_gene_panel.csv"))
    variants = load_variants(os.path.join(DATA_DIR, "sample_patient_variants.vcf"))
    control_keys = load_control_keys(os.path.join(DATA_DIR, "control_cohort_variants.vcf"))
    panel_variants = filter_to_panel(variants, panel)
    novel, in_controls = remove_control_variants(panel_variants, control_keys)
    assert len(novel) + len(in_controls) == len(panel_variants)
    assert all(v.key not in control_keys for v in novel)


def test_classifier_rules():
    de_novo_missense = Variant("1", 1, "A", "G", gene="X", consequence="missense", inheritance="de_novo")
    inherited_missense = Variant("1", 2, "A", "G", gene="X", consequence="missense", inheritance="inherited")
    nonsense = Variant("1", 3, "A", "G", gene="X", consequence="nonsense", inheritance="inherited")
    synonymous = Variant("1", 4, "A", "G", gene="X", consequence="synonymous", inheritance="inherited")
    known_mito = Variant("MT", 3243, "A", "G", var_type="MITO")
    unknown_mito = Variant("MT", 9999, "A", "G", var_type="MITO")

    assert classify_variant(de_novo_missense) == CLASS_P_LP
    assert classify_variant(inherited_missense) == CLASS_VUS
    assert classify_variant(nonsense) == CLASS_P_LP
    assert classify_variant(synonymous) == CLASS_BENIGN
    assert classify_variant(known_mito) == CLASS_P_LP
    assert classify_variant(unknown_mito) == CLASS_VUS


def test_classify_all_attaches_attribute():
    variants = [Variant("1", 1, "A", "G", consequence="synonymous")]
    classify_all(variants)
    assert hasattr(variants[0], "classification")
    assert variants[0].classification == CLASS_BENIGN


# --- Edge cases / stress-test style scenarios -------------------------------

def test_missing_patient_file_raises_clear_error():
    try:
        load_variants("data/this_file_does_not_exist.vcf")
        assert False, "expected VariantFileError"
    except VariantFileError as e:
        assert "not found" in str(e)


def test_malformed_lines_are_skipped_not_fatal():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\tINFO\n")
        f.write("13\t110850123\t.\tG\tA\tGENE=COL4A1;TYPE=SNV;CONSEQUENCE=missense\n")
        f.write("this line is garbage and too short\n")           # too few columns
        f.write("13\tNOT_A_NUMBER\t.\tG\tA\tGENE=COL4A1\n")        # bad POS
        f.write("13\t110850124\t.\tG\tA\n")                        # no INFO field at all
        path = f.name
    try:
        variants = load_variants(path)
        # 2 good rows survive (one with INFO, one without); 2 bad rows are skipped
        assert len(variants) == 2
        assert variants[1].gene == "NONE"  # missing INFO defaults safely
    finally:
        os.unlink(path)


def test_empty_gene_panel_raises_clear_error():
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write("gene,chrom,start,end,associated_condition\n")  # header only, no rows
        path = f.name
    try:
        panel = load_gene_panel(path)
        assert panel == {}
    finally:
        os.unlink(path)


def test_coordinate_fallback_catches_missing_gene_tag():
    panel = load_gene_panel(os.path.join(DATA_DIR, "cp_gene_panel.csv"))
    # Variant with a blank/wrong GENE tag but a position inside COL4A1's range
    v = Variant("13", 110850123, "G", "A", gene="NONE")
    filtered = filter_to_panel([v], panel)
    assert len(filtered) == 1
    assert filtered[0].gene == "COL4A1"
    assert filtered[0].match_method == "coordinate"


def test_variant_outside_any_panel_range_is_dropped():
    panel = load_gene_panel(os.path.join(DATA_DIR, "cp_gene_panel.csv"))
    v = Variant("7", 1000000, "G", "A", gene="NONE")
    filtered = filter_to_panel([v], panel)
    assert filtered == []


def test_dedupe_removes_exact_duplicate_rows():
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from main import _dedupe
    v1 = Variant("13", 110850123, "G", "A", gene="COL4A1")
    v2 = Variant("13", 110850123, "G", "A", gene="COL4A1")  # exact duplicate
    v3 = Variant("13", 110850124, "G", "A", gene="COL4A1")  # different position
    result = _dedupe([v1, v2, v3])
    assert len(result) == 2
