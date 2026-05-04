# Background-Matched MALDI-AMR Audit

This repository contains a working snapshot of a MALDI-TOF antimicrobial
resistance evaluation framework. The central idea is that MALDI-AMR models
should not be judged only by raw AUC: apparent resistance prediction can reflect
focal-drug signal, co-resistance structure, lineage/population structure, and
hospital-specific resistant ecology.

The main method in this repository is a model-agnostic background-matched audit:
given isolate-level model predictions and AST labels for multiple antibiotics,
it asks whether focal-drug prediction survives after matching isolates by their
co-resistance background.

## Repository Layout

```text
Mega_Model.py
  Main experimental model engine for DRIAMS MALDI-AMR training/evaluation.

run_background_audit_framework.py
  Model-agnostic background-matched audit. This is the core framework script.

scripts/
  run_training_clinical4.py
  run_training_ecoli6.py
  run_lgbm_baselines.py
  export_mega_predictions_for_audit.py
  run_background_audit.py
  run_public_upec_analysis.py
  make_paper_figures.py
  make_final_framework_tables_figures.py
  export_weis_predictions_for_audit.py
  marisma_end_to_end_kaggle.py
  background_matched_contrastive_kaggle.py
  build_cross_resistance_network.py

tests/
  Regression tests for the audit framework and Mega_Model helper behavior.

outputs/
  Current tables, figures, and intermediate analysis artifacts.

data_manifests/
  Public UPEC/WGS bridge metadata and file manifests. Raw DRIAMS spectra and
  raw UPEC spectrum archives are not included.

docs/
  Reproducibility guide, data/source notes, and longer audit notes.

manuscript/
  Overleaf-ready Nature Communications-style manuscript draft, vector PDF
  figures, LaTeX tables, references, and literature-positioning notes.
```

## Key Script Labels

- `scripts/export_mega_predictions_for_audit.py` — **Mega/CNN locked-prediction exporter**. Use this after a completed Mega run to create `mega_predictions_long.csv`, the isolate-level prediction table required by the background-matched audit.
- `run_background_audit_framework.py` — **model-agnostic audit engine**. Use this on any long prediction CSV from Mega/CNN, LGBM, Weis-style models, or external models.
- `scripts/run_background_audit.py` — **thin audit wrapper** for the default prediction CSV format.
- `scripts/marisma_end_to_end_kaggle.py` — **MARISMa external stress-test workflow**. Vectorizes MARISMa Bruker spectra, exports Mega/CNN predictions, and runs the same audit on the external MARISMa snapshot.

## Core Claim

MALDI-TOF AMR models are background-sensitive predictors. Their apparent AMR
performance can come from a mixture of focal-drug resistance signal,
co-resistance structure, lineage/population structure, and hospital-specific
resistant ecology.

The audit asks:

> If two isolates have the same co-resistance background, can the model still
> distinguish resistant versus susceptible isolates for the focal drug?

If the answer is yes, the model retains focal/background-independent signal. If
the answer is no, raw AUC was likely inflated by resistant-population background.

## Minimal Workflow

For the fuller command-by-command version, see
[`docs/reproduce.md`](docs/reproduce.md). Source citations and redistribution
notes are in [`docs/sources.md`](docs/sources.md) and
[`docs/data_availability.md`](docs/data_availability.md).

Train or evaluate a Mega model:

```bash
python scripts/run_training_ecoli6.py --data-root /path/to/driams
```

Export isolate-level predictions from a completed Mega run:

```bash
python scripts/export_mega_predictions_for_audit.py \
  --run-dir runs/exp_ecoli_mechanism6_drugid_mae30
```

Run the background-matched audit:

```bash
python scripts/run_background_audit.py \
  --predictions-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --output-dir outputs/background_audit
```

Run the public WGS-linked UPEC support analyses:

```bash
python scripts/run_public_upec_analysis.py \
  --median-peaks /path/to/Bruker_csv_medianpeaks_df.csv
```

Build the current paper tables/figures:

```bash
python scripts/make_paper_figures.py
```

Build the Nature Communications-style manuscript figures and LaTeX tables:

```bash
python scripts/make_ncomms_figures.py
```

## Required Prediction Format

The model-agnostic audit expects a long CSV with one row per isolate/drug
prediction:

```text
isolate_id,site,year,organism,drug,label,prob
```

Optional:

```text
model_name
```

The framework then builds co-resistance background signatures internally using
the other drug labels available for the same isolate.

## Current Evidence Snapshot

The included outputs summarize:

- Raw external AUC versus background-centered AUC for CNN and LGBM variants.
- Cross-resistance network structure in the E. coli expanded panel.
- Public WGS-linked UPEC evidence that MALDI spectra encode ST131 lineage.
- Published ST131 biomarker overlap/enrichment for discriminative MALDI peaks.
- MARISMa external stress-test outputs showing that the DRIAMS-trained Mega/CNN
  model does not transfer as-is to the current MARISMa E. coli snapshot; the
  audit flags all six E. coli targets as weak raw signal.

These outputs are analysis artifacts, not raw clinical data.

## Data Availability Notes

This repository does not include raw DRIAMS spectra, raw AST exports, raw WGS
FASTQ files, or raw Bruker archives. Reproducing the full pipeline requires
access to the corresponding public or controlled datasets and local paths passed
through the command-line wrappers.

The repository does include derived analysis outputs, public UPEC manifest
tables, processed Bruker median-peak features for the public support analysis,
and a locked Mega CNN checkpoint archive for the current snapshot.

## Status

This is a research snapshot prepared for organization and sharing. Before formal
submission, the code should be further cleaned into installable modules and the
main claims should be tied to final locked result files.
