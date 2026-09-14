# CP Variant Scanner — Sprint 1 Tool

BMI 461 — Sprint project inspired by:
> Fehlings, D. L., Zarrei, M., et al. (2024). *Comprehensive whole-genome sequence
> analyses provide insights into the genomic architecture of cerebral palsy.*
> Nature Genetics.

## What this is

The paper sequenced 327 CP trios and screened for pathogenic variants across
a target gene panel, filtering out variation that also appeared in healthy
control cohorts (CHILD, Inova), then classified what was left as
Pathogenic/Likely Pathogenic (P/LP), a Variant of Uncertain Significance
(VUS), or benign.

This tool reproduces that same **workflow shape** — gene-panel filter →
control-cohort filter → pathogenicity classification → cohort report — as a
small, runnable pipeline. It's meant as a *starting point* for the sprint
("Adaptation for Sprint" in the presentation): a way to quickly flag
candidate variants in new or siloed genetic data using the gene panel and
logic established by the paper, rather than a clinical-grade replication of
their actual pipeline (which used BWA, GATK, CNVnator/ERDS, Manta/DELLY,
ExpansionHunter, and 5-of-7 pathogenicity predictor consensus on real,
access-controlled patient data we don't have).

## How it maps to the paper

| Paper's step | This tool's equivalent |
|---|---|
| WGS variant calling (GATK, Manta, DELLY, ExpansionHunter, MitoMaster) | `vcf_parser.py` reads pre-called variants (SNV/INDEL/CNV/MITO tags in the INFO field stand in for those tools' output) |
| Candidate gene panel (COL4A1 and other repeatedly-hit genes) | `data/cp_gene_panel.csv` + `gene_filter.py` |
| Comparison against CHILD/Inova healthy control cohorts | `data/control_cohort_variants.vcf` + `control_filter.py` |
| 5-of-7 pathogenicity predictor consensus (SIFT, PolyPhen, CADD, etc.) | `classifier.py` — a transparent, swappable rule-based heuristic (consequence + inheritance pattern) |
| Diagnostic yield (% of cohort explained by a P/LP variant) | `report.compute_diagnostic_yield()` |
| Variant-type breakdown figure | `report.plot_variant_breakdown()` |

Because the paper's actual patient data (EGA, Brain-CODE, Inova) requires a
formal Data Access Committee request, everything here runs on small
**synthetic sample data** in `data/` so the pipeline is runnable and
testable out of the box. Swap in real (access-approved) VCFs and a real
annotation tool's output and the pipeline logic doesn't need to change.

## Project structure

```
cp-variant-scanner/
├── main.py                          # CLI entry point
├── requirements.txt
├── data/
│   ├── cp_gene_panel.csv            # candidate CP-associated genes
│   ├── sample_patient_variants.vcf  # synthetic patient trio variants
│   └── control_cohort_variants.vcf  # synthetic healthy-control variants
├── src/
│   ├── vcf_parser.py                # parses variant + control files
│   ├── gene_filter.py               # restricts to the candidate gene panel
│   ├── control_filter.py            # removes variants seen in controls
│   ├── classifier.py                # P/LP vs VUS vs Benign heuristic
│   └── report.py                    # tables, diagnostic yield, pie chart
├── tests/
│   └── test_pipeline.py             # unit tests for every module
└── output/
    └── variant_breakdown.png        # generated on each run
```

## Running it

```bash
pip install -r requirements.txt

python main.py \
  --patients data/sample_patient_variants.vcf \
  --controls data/control_cohort_variants.vcf \
  --panel data/cp_gene_panel.csv \
  --out output/
```

To run it on more than one patient (this is where "diagnostic yield" becomes
meaningful — the paper's 11.3% figure is a cohort statistic, not a
single-patient one):

```bash
python main.py --patients data/patient_001.vcf data/patient_002.vcf ...
```

### Tests

```bash
python -m pytest tests/
```

## Extending this for the sprint

A few directions the "Focusing on the Sprint" slide pointed at, and where
they'd plug in:

- **Swap the classifier for a real annotation tool.** `classifier.py` is
  isolated on purpose — replace `classify_variant()` with a call to a real
  ANNOVAR/CADD output parser without touching the rest of the pipeline.
- **Scale the control comparison.** Right now `control_filter.py` does exact
  (chrom, pos, ref, alt) matching against one small file. For real data,
  this is where you'd instead query a population frequency database (e.g.
  gnomAD) and filter by allele frequency instead of exact match.
- **Automate new-patient screening.** As the presentation notes, one gap in
  this research area is limited data volume. `main.py` already accepts
  multiple `--patients` files, so pointing it at a folder of new,
  unlabeled samples and checking which land on the gene panel is a direct
  way to help find more candidate patients/carriers for follow-up.
