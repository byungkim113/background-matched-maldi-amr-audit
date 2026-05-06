# Background-Matched Transfer Audit

This is a model-agnostic audit for MALDI-TOF AMR prediction. It does not train a
model and does not depend on Mega_Model. It only needs a prediction table from
any model, including Weis-style logistic regression, random forest, SVM,
LightGBM, CNNs, or future hospital models.

## Purpose

Raw external AUC can be high because a model learned the resistant population
background rather than focal-drug resistance. This audit asks a stricter
question:

> After matching isolates by co-resistance background, does prediction for the
> focal drug still work?

If the signal survives, the model has evidence of focal-drug signal. If it
collapses, raw performance was likely background/ecology driven.

## Required Input

The easiest input is a long CSV with one row per isolate/drug prediction:

```text
isolate_id,site,year,organism,drug,label,prob
A1,DRIAMS-A,2018,Escherichia coli,Ciprofloxacin,S,0.14
A1,DRIAMS-A,2018,Escherichia coli,Ceftriaxone,R,0.81
A2,DRIAMS-B,2018,Escherichia coli,Ciprofloxacin,R,0.77
```

Required meanings:

- `isolate_id`: stable isolate/sample identifier
- `site`: hospital, lab, country, or dataset split
- `year`: collection year; use a constant like `unknown` if unavailable
- `organism`: organism/species name
- `drug`: antibiotic name
- `label`: true AST result; accepts `S/R`, `0/1`, `susceptible/resistant`
- `prob`: model score where higher means more resistant

Column names are configurable, so you do not have to rename files by hand.

## Weis-Style Usage

The Weis repository is at:

<https://github.com/BorgwardtLab/maldi_amr>

The companion exporter in `scripts/export_weis_predictions_for_audit.py` can rerun the upstream Weis/Borgwardt model code and write one row per isolate/drug prediction. For a broad Weis-code rerun, use `--panel weis-core --external-row-policy all`. For the E. coli six-drug panel used in this project, use `--panel custom` with the relevant `--species` and `--drugs`.

After running a Weis-style or Weis-code model, export one row per isolate/drug prediction.
If the output columns are:

```text
sample_id,hospital,collection_year,species,antibiotic,phenotype,score
```

run:

```bash
python run_background_audit_framework.py \
  --predictions-csv weis_predictions.csv \
  --id-col sample_id \
  --site-col hospital \
  --year-col collection_year \
  --organism-col species \
  --drug-col antibiotic \
  --label-col phenotype \
  --prob-col score \
  --output-dir weis_background_audit
```

For stricter matching within the same collection year:

```bash
python run_background_audit_framework.py \
  --predictions-csv weis_predictions.csv \
  --id-col sample_id \
  --site-col hospital \
  --year-col collection_year \
  --organism-col species \
  --drug-col antibiotic \
  --label-col phenotype \
  --prob-col score \
  --match-year \
  --output-dir weis_background_audit_match_year
```

## Outputs

The script writes:

```text
normalized_predictions.csv
background_matched_audit_summary.csv
background_matched_audit_summary.md
background_matched_retained_rows.csv
background_matched_sensitivity.csv
cross_resistance_edges.csv
cross_resistance_prevalence.csv
background_audit_with_resistance_ecology.csv
background_audit_report.md
cross_resistance_network.svg
```

Key columns in `background_matched_audit_summary.csv`:

- `raw_auc`: ordinary AUC for the model
- `matched_auc`: AUC after keeping only background strata containing both R and S
- `stratum_centered_auc`: AUC after subtracting each background stratum's mean score
- `matched_retention`: fraction of rows retained after matching
- `adequacy_label`: whether the matched subset is interpretable
- `interpretation_category`: plain-language category

## Interpretation

| Pattern | Interpretation |
| --- | --- |
| raw AUC high, centered AUC high | focal-drug signal survives background matching |
| raw AUC high, centered AUC low | raw signal likely background/ecology driven |
| raw AUC low, centered AUC low | little useful signal |
| low retention / low n | do not overinterpret |

## Recommended Reporting

For a paper or model card, report:

1. Raw external AUC
2. Background-matched AUC
3. Stratum-centered AUC
4. Matched retention
5. Adequacy label
6. Cross-resistance network edges for the panel

This makes clear whether a MALDI-AMR model is predicting the focal resistance
phenotype or exploiting the population structure surrounding resistance.
