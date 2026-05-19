# Reproducibility Guide

This guide separates three levels of reproduction:

1. Run the model-agnostic audit on an existing prediction CSV.
2. Regenerate the current paper-facing tables and figures from committed outputs.
3. Reproduce training/evaluation from raw DRIAMS data.

Most reviewers should start with levels 1 and 2. Level 3 requires local access to raw DRIAMS data.

## 1. Environment Setup

Pip setup (recommended for most users):

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

RDKit is only required for Morgan-fingerprint drug conditioning. Install it separately via conda-forge (`conda install -c conda-forge rdkit`) if you need that path, as it is not included in `requirements.txt`.

## 2. Verify The Repository

```bash
python -m pytest tests
```

Syntax-check the main scripts:

```bash
PYTHONPYCACHEPREFIX=/tmp/pycache python -m py_compile \
  run_background_audit_framework.py \
  scripts/run_background_audit.py \
  scripts/export_mega_predictions_for_audit.py \
  scripts/export_weis_predictions_for_audit.py \
  scripts/compare_weis_raw_metrics.py \
  scripts/upec_wgs_validation_analysis.py \
  scripts/updated_proteomic_overlap_analysis.py \
  scripts/run_public_upec_analysis.py \
  scripts/marisma_end_to_end_kaggle.py \
  scripts/make_final_framework_tables_figures.py \
  scripts/make_paper_figures.py
```

## 3. Run The Audit On Any Model

Prepare a long prediction CSV with one row per isolate/drug prediction:

```text
isolate_id,site,year,organism,drug,label,prob
```

Then run:

```bash
python scripts/run_background_audit.py \
  --predictions-csv /path/to/predictions_long.csv \
  --output-dir outputs/background_audit_custom \
  --bootstrap-n 500 \
  --permutation-n 500
```

If your CSV uses different column names, call the core script directly:

```bash
python run_background_audit_framework.py \
  --predictions-csv /path/to/predictions_long.csv \
  --output-dir outputs/background_audit_custom \
  --id-col sample_id \
  --site-col hospital \
  --year-col collection_year \
  --organism-col species \
  --drug-col antibiotic \
  --label-col phenotype \
  --prob-col score
```

Main outputs:

- `background_matched_audit_summary.csv`
- `background_matched_retained_rows.csv`
- `background_matched_sensitivity.csv`
- `background_audit_with_resistance_ecology.csv`
- `cross_resistance_edges.csv`
- `background_audit_report.md`

## 4. Export Mega/CNN Predictions From A Locked Run

Use this when you have a completed Mega run directory with `config.json` and model checkpoints:

```bash
python scripts/export_mega_predictions_for_audit.py \
  --run-dir runs/exp_ecoli_mechanism6_drugid_mae30 \
  --data-root /path/to/driams \
  --output-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --model-name Mega-CNN
```

Then audit:

```bash
python scripts/run_background_audit.py \
  --predictions-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --output-dir outputs/mega_cnn_background_audit
```

## 5. Train The Clinical4 Mega/CNN Experiment

```bash
python scripts/run_training_clinical4.py \
  --data-root /path/to/driams \
  --output-dir runs \
  --experiment exp_clinical4_mae30 \
  --mae-epochs 30 \
  --early-stop macro \
  --seed-policy all
```

This runs the four-pair clinical panel:

- *E. coli* / Ciprofloxacin
- *E. coli* / Amoxicillin-Clavulanic acid
- *S. aureus* / Oxacillin
- *S. epidermidis* / Erythromycin

## 6. Train The E. coli Mechanism6 Mega/CNN Experiment

```bash
python scripts/run_training_ecoli6.py \
  --data-root /path/to/driams \
  --output-dir runs \
  --experiment exp_ecoli_mechanism6_drugid_mae30 \
  --mae-epochs 30 \
  --early-stop primary \
  --seed-policy all \
  --prevalence-shift none
```

This runs the six-pair E. coli panel:

- Ciprofloxacin
- Norfloxacin
- Amoxicillin-Clavulanic acid
- Ceftriaxone
- Ceftazidime
- Cefepime

## 7. Run LightGBM Baselines

```bash
python scripts/run_lgbm_baselines.py \
  --data-root /path/to/driams \
  --pair-profile ecoli_mechanism6 \
  --output-dir runs \
  --experiment lgbm_ecoli_mechanism6 \
  --with-random-cv
```

This produces `runs/lgbm_ecoli_mechanism6/metrics/lgbm_results.csv` and, if requested, the temporal-versus-random-CV diagnostic.

## 8. Run Public UPEC WGS/Proteomic Support Analyses

Using the processed Bruker median-peak feature table committed in this repository:

```bash
python scripts/run_public_upec_analysis.py \
  --metadata data_manifests/upec_master_metadata.tsv \
  --median-peaks data_manifests/Bruker_csv_medianpeaks_df.csv \
  --wgs-output-dir outputs/analysis_outputs/upec_wgs_validation_outputs \
  --proteomic-output-dir outputs/analysis_outputs/updated_proteomic_overlap_outputs \
  --folds 5 \
  --permutations 10000
```

Key outputs:

- `outputs/analysis_outputs/upec_wgs_validation_outputs/centroid_binary_cv_results.csv`
- `outputs/analysis_outputs/upec_wgs_validation_outputs/st131_resistance_associations.csv`
- `outputs/analysis_outputs/updated_proteomic_overlap_outputs/updated_proteomic_overlap_permutation_enrichment.csv`

## 9. Regenerate Final Tables And Figures

```bash
python scripts/make_paper_figures.py --output-dir /tmp/maldi_amr_final_framework_outputs
```

This writes reproduced artifacts to a scratch folder so a clean clone is not
modified. The committed paper-facing artifact folder is:

```text
outputs/final_framework_outputs/
```

It contains the current paper-facing tables, figures, and argument draft.

## 10. Weis/Borgwardt Published-Code Audit

The notebook:

```text
notebooks/weis_lightgbm_background_audit_kaggle.ipynb
```

is a Kaggle-oriented workflow. It clones the Borgwardt/Weis code, exports
predictions with isolate IDs, then runs the same model-agnostic background
audit. See `docs/weis_published_model_audit.md` for the exact claim hierarchy:
exact replication, Weis-code rerun, and background-matched published-model
audit.

For a paper-parity check, use the upstream repository, the original
organism/drug panel and the stratified external subset:

```bash
git clone https://github.com/BorgwardtLab/maldi_amr.git /kaggle/working/maldi_amr

python scripts/export_weis_predictions_for_audit.py \
  --weis-repo /kaggle/working/maldi_amr \
  --driams-root /kaggle/input/datasets/drscarlat/driams \
  --audit-script scripts/run_background_audit.py \
  --panel weis-core \
  --model lr \
  --external-row-policy stratified \
  --seed 35 \
  --n-folds 5 \
  --bootstrap-n 200 \
  --permutation-n 200 \
  --output-dir outputs/weis_lr_paper_parity_rerun
```

For the main background-matched audit, use the same official model code but
score all eligible external isolates:

```bash
python scripts/export_weis_predictions_for_audit.py \
  --weis-repo /kaggle/working/maldi_amr \
  --driams-root /kaggle/input/datasets/drscarlat/driams \
  --audit-script scripts/run_background_audit.py \
  --panel weis-core \
  --model lr \
  --external-row-policy all \
  --seed 35 \
  --n-folds 5 \
  --bootstrap-n 500 \
  --permutation-n 500 \
  --output-dir outputs/weis_lr_full_external_audit
```

Check paper-parity metrics against the official stored JSONs before using
"replication" language:

```bash
python scripts/compare_weis_raw_metrics.py \
  --weis-raw-results outputs/weis_lr_paper_parity_rerun/weis_raw_results.json \
  --reference-results-root /kaggle/working/maldi_amr/results/validation_per_species_and_antibiotic/lr \
  --output-csv outputs/weis_lr_paper_parity_rerun/weis_metric_parity.csv \
  --summary-md outputs/weis_lr_paper_parity_rerun/weis_metric_parity.md
```

Use `--external-row-policy all` for paper-facing background matching. The legacy `--external-row-policy stratified` option intentionally scores only the stratified external subset and can leave too few rows for matched strata, especially at DRIAMS-B.

Repeat with `--model lightgbm` if the audited published-model claim is about
the Weis/Borgwardt LightGBM family rather than logistic regression.

Use the custom E. coli six-drug panel only for comparability with this
repository's primary E. coli mechanism panel:

```bash
python scripts/export_weis_predictions_for_audit.py \
  --weis-repo /kaggle/working/maldi_amr \
  --driams-root /kaggle/input/datasets/drscarlat/driams \
  --audit-script scripts/run_background_audit.py \
  --panel custom \
  --species "Escherichia coli" \
  --drugs "Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,Ceftriaxone,Ceftazidime,Cefepime" \
  --model lightgbm \
  --external-row-policy all \
  --seed 35 \
  --n-folds 2 \
  --output-dir outputs/weis_lightgbm_ecoli6_external_audit
```

Treat these as Weis-code reruns. Exact publication parity still requires comparing the raw metrics against the upstream Weis result JSONs and confirming the same preprocessing, splits, and hyperparameter search space.

To audit multiple upstream model families, rerun the same command with different `--model` values such as `lr`, `rf`, `svm-linear`, `svm-rbf`, `lightgbm`, and `mlp`, writing each model to its own output directory.

## 11. MARISMa External Stress Test

The script:

```text
scripts/marisma_end_to_end_kaggle.py
```

is a Kaggle-oriented external validation workflow for the MARISMa Bruker data.
It has three stages:

1. `preprocess`: read MARISMa `AMR.csv`, locate raw Bruker folders, and convert
   spectra to 6,000-bin vectors over 2,000-20,000 Da.
2. `predict`: apply a locked Mega/CNN checkpoint to the MARISMa vectors and write
   a long prediction CSV.
3. `audit`: run the model-agnostic background-matched audit on those MARISMa
   predictions.

Smoke-test preprocessing:

```bash
python scripts/marisma_end_to_end_kaggle.py \
  --stage preprocess \
  --amr-csv /kaggle/input/datasets/bfdf121/marisma2/AMR.csv \
  --marisma-root /kaggle/input/datasets/bfdf121/marisma/MARISMa \
  --output-dir /kaggle/working/marisma_preprocessed_smoke \
  --max-spectra 50
```

Full preprocessing:

```bash
python scripts/marisma_end_to_end_kaggle.py \
  --stage preprocess \
  --amr-csv /kaggle/input/datasets/bfdf121/marisma2/AMR.csv \
  --marisma-root /kaggle/input/datasets/bfdf121/marisma/MARISMa \
  --output-dir /kaggle/working/marisma_preprocessed_full
```

Prediction from a locked Mega/CNN run:

```bash
python scripts/marisma_end_to_end_kaggle.py \
  --stage predict \
  --mega-model-path /kaggle/working/Mega_Model.py \
  --run-dir /kaggle/input/datasets/bfdf121/newruns/runs/exp_ecoli_mechanism6_drugid_mae30 \
  --vectors-npy /kaggle/input/datasets/bfdf121/marisma-vectors/results-5/marisma_preprocessed_full/marisma_vectors_6000.npy \
  --manifest-csv /kaggle/input/datasets/bfdf121/marisma-vectors/results-5/marisma_preprocessed_full/marisma_prediction_manifest.csv \
  --prediction-csv /kaggle/working/marisma_mega_predictions_long.csv \
  --tta-passes 1
```

Audit:

```bash
python scripts/marisma_end_to_end_kaggle.py \
  --stage audit \
  --prediction-csv /kaggle/working/marisma_mega_predictions_long.csv \
  --output-dir /kaggle/working/marisma_external_validation \
  --audit-output-dir /kaggle/working/marisma_external_validation/marisma_isolate_background_audit
```

The audit stage aggregates spot-level MARISMa rows to isolate/drug rows before
running the framework, excludes isolate/drug groups with conflicting labels, and
writes `marisma_duplicate_handling_report.json`.

The current locked snapshot is committed under:

```text
outputs/analysis_outputs/marisma_external_validation/
```

In that snapshot, raw MARISMa AUCs for the six E. coli mechanism6 targets are
near chance. The audit therefore labels this as an external stress-test failure,
not a positive deployment validation.
