# Prediction Exporters

This repository separates model training from model auditing. The background-matched audit does not need model checkpoints directly; it needs locked isolate-level prediction CSVs.

## Mega/CNN Locked-Prediction Exporter

Use:

```text
scripts/export_mega_predictions_for_audit.py
```

Label:

```text
Mega/CNN locked-prediction exporter
```

Purpose:

- Loads a completed `Mega_Model.py` run.
- Reads the run `config.json`.
- Loads the selected `maldi_amr_seed*.pt` checkpoints.
- Runs forward passes on A-2018 and external DRIAMS sites.
- Writes one long CSV row per isolate/drug prediction.

Expected output:

```text
runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv
```

Required audit columns:

```text
isolate_id,site,year,organism,drug,label,prob
```

Typical command:

```bash
python scripts/export_mega_predictions_for_audit.py \
  --run-dir runs/exp_ecoli_mechanism6_drugid_mae30 \
  --data-root /path/to/driams \
  --output-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --model-name Mega-CNN
```

Then run:

```bash
python scripts/run_background_audit.py \
  --predictions-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --output-dir outputs/mega_cnn_background_audit
```

## Portable/External Model Exporters

For non-Mega models, any script is acceptable if it writes the same long prediction schema. The audit is intentionally model-agnostic.

The Weis/Borgwardt published-code exporter is:

```text
scripts/export_weis_predictions_for_audit.py
```

That script is supplementary and Kaggle-oriented because the upstream
Weis/Borgwardt workflow has its own dependencies and runtime assumptions. It
reruns the original BorgwardtLab `maldi_amr` model code and writes isolate-level
predictions that the background audit can consume. See
`docs/weis_published_model_audit.md` before describing the output as an exact
replication.

Recommended full external-row audit command:

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
  --output-dir outputs/weis_lr_full_external_audit
```

Important options:

- `--panel weis-core` uses the organism/drug panel from the original Weis/Borgwardt repository.
- `--panel custom` uses the supplied `--species` and `--drugs` arguments, which is useful for the E. coli six-drug panel used in this project.
- `--external-row-policy all` scores all eligible external rows and is the recommended setting for background matching.
- `--external-row-policy stratified` reproduces the older small external subset diagnostic and should not be used for primary audit figures.
- `--model` is passed directly into upstream `amr_maldi_ml.models.run_experiment`; original supported values include `lr`, `svm-rbf`, `svm-linear`, `rf`, `lightgbm`, and `mlp`. Run separate exports for separate model families.

The output includes `weis_reproduction_report.md` and
`weis_reproduction_report.json`, including the upstream git commit when
available. Treat these results as a Weis-code rerun until the raw AUCs are
checked against the upstream Weis result JSONs for exact publication parity.
