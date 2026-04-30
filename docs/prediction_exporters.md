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

The Weis/Borgwardt compatibility script is:

```text
scripts/export_weis_predictions_for_audit.py
```

That script is supplementary and Kaggle-oriented because the upstream Weis/Borgwardt workflow has its own dependencies and runtime assumptions.
