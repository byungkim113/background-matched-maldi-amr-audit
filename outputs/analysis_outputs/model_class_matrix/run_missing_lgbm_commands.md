# Missing LGBM Model-Class Cells

Run this on Kaggle or another machine with the DRIAMS data mounted. The pipeline
exports isolate-level LGBM predictions, runs the model-agnostic background audit,
and rebuilds the matrix.

```bash
python scripts/run_model_class_matrix_pipeline.py \
  --data-root /kaggle/input/datasets/drscarlat/driams \
  --ecoli-run-dir /kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30 \
  --saureus-run-dir /kaggle/working/runs/exp_saureus_panel_oxa_background_mae30
```

Use `--dry-run` first to print every command without running it.

The equivalent manual commands are below.

## E. coli LGBM single-task

```bash
python scripts/export_lgbm_predictions_for_audit.py \
  --data-root /kaggle/input/datasets/drscarlat/driams \
  --pair-profile ecoli_mechanism6 \
  --run-dir /kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30 \
  --variants single \
  --output-dir /kaggle/working/lgbm_exports/ecoli

python run_background_audit_framework.py \
  --predictions-csv /kaggle/working/lgbm_exports/ecoli/lgbm_single_predictions_long.csv \
  --background-signature-col background_signature \
  --model-name LGBM-single-ecoli6 \
  --output-dir /kaggle/working/ecoli_lgbm_single_background_audit
```

## S. aureus/Oxacillin LGBM single-task and multi-task

```bash
python scripts/export_lgbm_predictions_for_audit.py \
  --data-root /kaggle/input/datasets/drscarlat/driams \
  --pair-profile saureus_panel \
  --run-dir /kaggle/working/runs/exp_saureus_panel_oxa_background_mae30 \
  --variants single,multi \
  --train-if-missing \
  --output-dir /kaggle/working/lgbm_exports/saureus_oxa

python run_background_audit_framework.py \
  --predictions-csv /kaggle/working/lgbm_exports/saureus_oxa/lgbm_multi_predictions_long.csv \
  --background-signature-col background_signature \
  --model-name LGBM-multi-saureus-oxa \
  --output-dir /kaggle/working/saureus_lgbm_multi_oxa_background_audit

python run_background_audit_framework.py \
  --predictions-csv /kaggle/working/lgbm_exports/saureus_oxa/lgbm_single_predictions_long.csv \
  --background-signature-col background_signature \
  --model-name LGBM-single-saureus-oxa \
  --output-dir /kaggle/working/saureus_lgbm_single_oxa_background_audit
```
