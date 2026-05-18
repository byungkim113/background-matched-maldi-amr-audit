# Weis LR E. coli Six-Drug Background Audit

This directory records the background-matched audit summary from the multi-drug
E. coli Weis-LR prediction export.

## Run Scope

- Input prediction table: `weis_predictions_long.csv`
- Input table SHA256: `5affe3d28d49d64daaaf37a63a1aa1559851e5c313cb2499618c10224a869e29`
- Rows read by audit: 3,184
- Model label: `Weis-lr`
- Organism: Escherichia coli
- Drugs: Amoxicillin-Clavulanic acid, Cefepime, Ceftazidime, Ceftriaxone,
  Ciprofloxacin, Norfloxacin
- Test sites: DRIAMS-B, DRIAMS-C, DRIAMS-D
- Bootstrap iterations: 500
- Permutation iterations: 500
- Kaggle output directory: `/kaggle/working/weis_lr_ecoli6_background_audit`
- Summary CSV SHA256: `7794e192913b38357738292a60dc68d79b92e078e7fc0b3bb8779c2ea1dc5b3f`

The isolate-level prediction table is not committed here; this directory keeps
the aggregate audit summary used for interpretation.

## Summary

- Audit summary rows: 17
- Fully interpretable rows: 3
- Caution rows: 6
- Not interpretable rows: 8

Fully interpretable rows:

| Site | Drug | Raw AUC | Background-Centered AUC | Matched Retention | N Matched | Interpretation |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| DRIAMS-C | Amoxicillin-Clavulanic acid | 0.626248 | 0.578182 | 0.745856 | 135 | partially_retained_or_uncertain |
| DRIAMS-D | Amoxicillin-Clavulanic acid | 0.652425 | 0.524015 | 0.924812 | 369 | background_driven_collapse |
| DRIAMS-D | Ciprofloxacin | 0.667111 | 0.619607 | 0.914948 | 355 | focal_signal_retained |

## Interpretation

The multi-drug Weis-LR audit shows mixed behavior after matching on
co-resistance background. DRIAMS-D Ciprofloxacin retains a focal signal after
background matching, while DRIAMS-D Amoxicillin-Clavulanic acid collapses toward
chance after adjustment. Many other site-drug rows lack enough within-background
resistant and susceptible isolates to support a strong matched interpretation.

This result should be described as an audit of a Weis-style LR multi-drug export,
not as an additional exact Weis official-panel parity result.
