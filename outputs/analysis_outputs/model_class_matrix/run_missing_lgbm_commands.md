# LGBM Model-Class Cells Completed

All planned LGBM model-class audit cells have been exported, audited, and
included in `model_class_matrix.csv`.

Completed additions:

- E. coli expanded-panel LGBM single-task background audit.
- S. aureus/Oxacillin LGBM multi-task background audit.
- S. aureus/Oxacillin LGBM single-task background audit.

To reproduce these outputs on Kaggle or another machine with DRIAMS mounted,
run:

```bash
python scripts/run_model_class_matrix_pipeline.py --dry-run
python scripts/run_model_class_matrix_pipeline.py
```
