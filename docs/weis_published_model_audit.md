# Weis/Borgwardt Published-Model Audit Protocol

This protocol separates three claims that should not be mixed.

1. **Exact Weis et al. replication** means reproducing the paper's reported raw
   metrics under the same upstream code, DRIAMS preprocessing, split policy,
   model family, seed, cross-validation search and result files.
2. **Weis-code rerun** means training and scoring through the official
   `BorgwardtLab/maldi_amr` implementation, but not yet proving exact equality
   to the paper's stored results.
3. **Background-matched published-model audit** means taking isolate-level
   predictions from that Weis-code rerun and applying this repository's audit.
   For this purpose, scoring all eligible external isolates is preferable to
   scoring only the paper's small stratified external subset.

The current manuscript should only call the existing committed output a
Weis/Borgwardt-style compatibility audit. Promote it to a clean published-model
audit only after the steps below produce a traced prediction export.

## Required Inputs

- Official Weis/Borgwardt repository:
  `https://github.com/BorgwardtLab/maldi_amr`
- Raw DRIAMS Dryad extraction containing `DRIAMS-A`, `DRIAMS-B`, `DRIAMS-C` and
  `DRIAMS-D`.
- A Python environment that can import:
  `maldi_learn`, `libTLDA`, `lightgbm`, `sklearn`, `numpy` and `pandas`.

The exporter records the upstream git commit when `--weis-repo` is a git
checkout. Keep that report with the generated prediction CSV.

## Step 1: Paper-Parity Rerun

Use the original Weis/Borgwardt organism-drug panel and the stratified external
row policy. This run is for checking whether the local environment reproduces
the upstream raw AUROC/AUPRC values closely enough to cite as a replication.

```bash
python scripts/export_weis_predictions_for_audit.py \
  --weis-repo /path/to/BorgwardtLab/maldi_amr \
  --driams-root /path/to/driams \
  --audit-script run_background_audit_framework.py \
  --panel weis-core \
  --model lr \
  --external-row-policy stratified \
  --seed 35 \
  --n-folds 5 \
  --bootstrap-n 200 \
  --permutation-n 200 \
  --output-dir outputs/weis_lr_paper_parity_rerun
```

Repeat with `--model lightgbm` if the paper claim being checked is the
LightGBM table/figure rather than the logistic-regression figure.

The key files are:

- `weis_raw_results.json`
- `weis_predictions_long.csv`
- `weis_reproduction_report.json`
- `audit/background_matched_audit_summary.csv`

Do not use this run alone for the main background-matched conclusion if the
stratified subset leaves too few valid co-resistance strata.

Then compare the raw metrics with the upstream stored result JSONs:

```bash
python scripts/compare_weis_raw_metrics.py \
  --weis-raw-results outputs/weis_lr_paper_parity_rerun/weis_raw_results.json \
  --reference-results-root /path/to/BorgwardtLab/maldi_amr/results/validation_per_species_and_antibiotic/lr \
  --output-csv outputs/weis_lr_paper_parity_rerun/weis_metric_parity.csv \
  --summary-md outputs/weis_lr_paper_parity_rerun/weis_metric_parity.md \
  --tolerance 1e-6
```

Use the matching reference subdirectory for the selected model family, for
example `validation_per_species_and_antibiotic/lightgbm` for LightGBM.

## Step 2: Published-Code Audit Run

Use the same upstream code, model family, seed and panel, but score all eligible
external isolates so the audit has enough within-background support.

```bash
python scripts/export_weis_predictions_for_audit.py \
  --weis-repo /path/to/BorgwardtLab/maldi_amr \
  --driams-root /path/to/driams \
  --audit-script run_background_audit_framework.py \
  --panel weis-core \
  --model lr \
  --external-row-policy all \
  --seed 35 \
  --n-folds 5 \
  --bootstrap-n 500 \
  --permutation-n 500 \
  --output-dir outputs/weis_lr_full_external_audit
```

This is the cleanest statement for the paper if parity checks pass:

> We reran the official Weis/Borgwardt implementation to obtain isolate-level
> predictions and applied the background-matched audit to those predictions.

If parity checks have not passed, use this weaker statement:

> We applied the audit to a Weis/Borgwardt-code rerun; these rows demonstrate
> audit compatibility with a published workflow but are not claimed as an exact
> replication of Weis et al.

## Step 3: Promotion Criteria

Promote the result from compatibility audit to published-model audit only when:

- `weis_reproduction_report.json` records an official `BorgwardtLab/maldi_amr`
  commit.
- The chosen model family matches the claim (`lr` or `lightgbm`).
- Raw metrics from the parity rerun are close to the upstream stored/public
  metrics for the same model and organism-drug panel.
- The full-external audit run uses the same trained-code path and clearly states
  that `external_row_policy=all` was chosen for the audit, not for exact paper
  row-set replication.
- The manuscript says exactly which level is being claimed.

## Manuscript Language

Use this if parity is not yet established:

> As a compatibility analysis, we reran the official Weis/Borgwardt modeling
> code to produce isolate-level predictions and applied the same background
> audit. This analysis demonstrates that the audit can consume predictions from
> a published MALDI-AMR workflow, but it is not treated as an exact replication
> of the original benchmark.

Use this only after parity is established:

> We reproduced the Weis/Borgwardt shallow-model workflow sufficiently to match
> the reported raw performance for the audited model family, exported
> isolate-level predictions, and then applied the background-matched audit.
