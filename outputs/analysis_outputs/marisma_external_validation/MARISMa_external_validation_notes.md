# MARISMa External Validation Snapshot

This folder contains the locked MARISMa prediction/audit snapshot generated from
the Mega/CNN E. coli mechanism6 checkpoint on the MARISMa Bruker spectra and AMR
labels.

## What Was Run

- Input spectra: MARISMa Bruker raw spectra, vectorized to the same 6,000-bin
  2,000-20,000 Da format used by the DRIAMS Mega/CNN workflow.
- Input labels: `AMR.csv` from the MARISMa Kaggle dataset snapshot.
- Model: locked Mega/CNN E. coli mechanism6 run,
  `exp_ecoli_mechanism6_drugid_mae30`.
- Prediction output:
  `marisma_mega_predictions_long.csv`.
- Audit output:
  `marisma_background_audit/background_matched_audit_summary.csv`.

## Important Dataset Note

The AMR label file contains records from 2018-2024 overall, but the three target
organisms used here (*Escherichia coli*, *Staphylococcus aureus*, and
*Staphylococcus epidermidis*) appear only in 2024 in this AMR snapshot. The
current Mega/CNN external check therefore uses the 2024 MARISMa target-organism
labels, not a multi-year MARISMa time series.

## Result Summary

For all six E. coli targets in the mechanism6 panel, raw external AUC was near
chance. Background matching was interpretable for all rows, but it did not rescue
the model because the starting raw signal was weak.

| Drug | Raw AUC | Background-centered AUC | Matched retention | Audit category |
| --- | ---: | ---: | ---: | --- |
| Ciprofloxacin | 0.519 | 0.507 | 0.862 | weak raw signal |
| Norfloxacin | 0.541 | 0.475 | 0.503 | weak raw signal |
| Amoxicillin-Clavulanic acid | 0.519 | 0.538 | 0.964 | weak raw signal |
| Ceftriaxone/Cefotaxime analog | 0.478 | 0.536 | 0.144 | weak raw signal |
| Ceftazidime | 0.479 | 0.562 | 0.144 | weak raw signal |
| Cefepime | 0.473 | 0.416 | 0.205 | weak raw signal |

## Interpretation

This is not a positive deployment validation. It is an external stress test.
The audit catches that the DRIAMS-trained Mega/CNN model should not be trusted
directly on this MARISMa snapshot: raw discrimination is close to chance and the
background-centered results remain weak or unstable.

For the paper, the safe claim is:

> In an independent MARISMa Bruker snapshot, the DRIAMS-trained Mega/CNN model did
> not transfer directly. The background-matched audit therefore flags the model as
> unsuitable for deployment without MARISMa-specific recalibration, retraining, or
> additional validation.

Do not use this result to claim MARISMa confirms DRIAMS transferability. Its value
is that the audit produces a clear "do not deploy as-is" decision when a model is
externally brittle.
