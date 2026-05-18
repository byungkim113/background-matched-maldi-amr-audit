# Overleaf Copy/Paste Guide

Use this folder as the Overleaf project root.

## Recommended: master paper draft

Use `overleaf_master_paper.tex` as the main Overleaf document for the cleaned
paper draft. It is self-contained except for the figure PDFs.

Upload these figure files in a folder named `figures/`:

```text
figures/figure_1_framework.pdf
figures/figure_2_primary_background_audit.pdf
figures/figure_4_cross_resistance_network.pdf
figures/figure_5_public_wgs_proteomic_support.pdf
```

## Option A: Upload the manuscript folder

Upload the full `manuscript/` folder with this structure:

```text
main.tex
references.bib
figures/
tables/
source_data/
supplementary/
```

Then set `main.tex` as the main document.

## Option B: Copy/paste only the manuscript

1. Create a new blank Overleaf project.
2. Paste `main.tex` into the Overleaf `main.tex`.
3. Upload all PDFs from `figures/`.
4. Upload all `.tex` files from `tables/`.
5. Keep the same folder names:

```text
figures/figure_1_framework.pdf
figures/figure_2_primary_background_audit.pdf
figures/figure_3_model_family_replication.pdf
figures/figure_4_cross_resistance_network.pdf
figures/figure_5_public_wgs_proteomic_support.pdf
tables/table_1_primary_audit.tex
tables/table_2_model_replication.tex
tables/table_3_cross_resistance_edges.tex
tables/table_4_public_wgs_auc.tex
tables/table_5_biomarker_enrichment.tex
tables/table_6_upec_clone_control.tex
```

## Supplementary information

To compile the supplementary file, upload:

```text
supplementary/supplementary_information.tex
tables/table_6_claim_guardrails.tex
source_data/
```

Set `supplementary/supplementary_information.tex` as the main file when compiling the supplement.

## Regenerating figures locally

From the repository root:

```bash
python scripts/make_ncomms_figures.py
```

This regenerates:

- `manuscript/figures/*.pdf`
- `manuscript/tables/*.tex`
- `manuscript/source_data/*.csv`
