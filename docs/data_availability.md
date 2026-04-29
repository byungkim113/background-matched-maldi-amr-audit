# Data Availability And Redistribution

This repository is organized so reviewers can inspect the framework code, run the model-agnostic audit on prediction CSVs, and regenerate the included paper-facing tables and figures. It intentionally does not redistribute most raw clinical or sequencing data.

## Included In This Repository

- Source code for the background-matched audit framework.
- Mega/CNN training and export wrappers.
- LightGBM baseline wrapper.
- Public UPEC/Basel manifest tables and processed Bruker median-peak features used by the WGS-linked support analysis.
- Derived audit outputs, final tables, and final figures under `outputs/`.
- A locked Mega CNN checkpoint archive under `model_checkpoints/mega_cnn_archive_2026-04-22/`.

## Not Included

- Raw DRIAMS spectra, raw DRIAMS AST exports, or patient-level clinical records.
- Raw ENA FASTQ files from PRJEB55855.
- Raw Bruker FID archives or raw OSF spectrum archives.
- Any private hospital metadata beyond the derived, de-identified analysis tables already committed.

## External Data Required For Full Reproduction

To reproduce the full DRIAMS experiments from raw spectra, download the DRIAMS data independently and pass the local root with `--data-root`. The expected structure is a directory containing `DRIAMS-A`, `DRIAMS-B`, `DRIAMS-C`, and `DRIAMS-D`.

To reproduce the public UPEC support analysis from upstream files, download the Cuénod/Basel UPEC metadata and Bruker median-peak data from the upstream sources listed in [sources.md](sources.md). The repository includes the processed files used for the current result snapshot, but not the raw archives.

## Prediction CSVs Versus Checkpoints

The audit itself does not require model checkpoints. It only needs a long prediction CSV with:

```text
isolate_id,site,year,organism,drug,label,prob
```

For fastest review, provide or use a locked `mega_predictions_long.csv` generated from the final model. For full model reruns, use the training wrappers and the checkpoints or retrain from raw DRIAMS data.

## Privacy And Redistribution Notes

No raw patient-level data are intentionally included. All raw clinical datasets should be obtained from their original repositories or data providers under their own terms of use. Before making this repository fully public, re-check upstream licenses for every data-derived table and decide whether checkpoint hosting should remain in GitHub or move to Zenodo/OSF with a DOI.
