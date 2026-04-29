#!/usr/bin/env python3
"""Export Weis/Borgwardt MALDI-AMR predictions for background-matched auditing.

The Weis repository's stored JSON result files contain scores and labels, but
not isolate identifiers. Background-matched auditing needs isolate identifiers
and co-resistance background labels, so this script reruns the Weis-style
external validation workflow and writes an ID-preserving long prediction table.

Typical Kaggle usage:

    python export_weis_predictions_for_audit.py \
      --weis-repo /kaggle/working/maldi_amr \
      --driams-root /kaggle/input/datasets/drscarlat/driams \
      --species "Escherichia coli" \
      --drugs "Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,Ceftriaxone,Ceftazidime,Cefepime" \
      --model lr \
      --seed 35 \
      --output-dir /kaggle/working/weis_background_audit

The script writes:

    weis_predictions_long.csv
    weis_raw_results.json
    audit/background_matched_audit_summary.csv
    audit/background_audit_report.md
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import pathlib
import subprocess
import sys
from datetime import datetime
from typing import Iterable

import numpy as np
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score


LABEL_CHAR = {0: "S", 1: "R"}
UNKNOWN_CHAR = "U"


def import_from_path(module_name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {module_name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def split_csv_arg(text: str) -> list[str]:
    return [part.strip() for part in str(text or "").split(",") if part.strip()]


def load_weis_modules(weis_repo: pathlib.Path):
    ml_dir = weis_repo / "amr_maldi_ml"
    if not ml_dir.exists():
        raise FileNotFoundError(f"Weis amr_maldi_ml directory not found: {ml_dir}")
    sys.path.insert(0, str(ml_dir))
    print(f"[setup] Loading Weis model code from {ml_dir}", flush=True)
    models = import_from_path("weis_models_for_audit", ml_dir / "models.py")
    return models


def normalize_ast_value(value) -> int | None:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    upper = text.upper()
    if upper in {"R", "RESISTANT", "TRUE", "1"}:
        return 1
    if upper in {"S", "SUSCEPTIBLE", "FALSE", "0"}:
        return 0
    try:
        number = float(text)
    except ValueError:
        return None
    if number == 0:
        return 0
    if number == 1:
        return 1
    return None


def code_series_from_dataset(dataset) -> list[str]:
    meta = dataset.y
    for col in ("code", "Code", "sample_id", "id", "uuid"):
        if col in meta.columns:
            return [str(value) for value in meta[col].values]
    return [str(idx) for idx in meta.index]


def load_dataset_for_focal(
    *,
    driams_root: str,
    site: str,
    years: str | Iterable[str],
    species: str,
    antibiotic: str,
):
    from maldi_learn.driams import load_driams_dataset
    from maldi_learn.filters import DRIAMSBooleanExpressionFilter
    from maldi_learn.utilities import case_based_stratification, stratify_by_species_and_label

    print(
        f"[data] Loading focal dataset site={site} years={years} "
        f"species={species} drug={antibiotic}",
        flush=True,
    )

    extra_filters = []
    if site == "DRIAMS-A":
        extra_filters.append(DRIAMSBooleanExpressionFilter("workstation != HospitalHygiene"))
        id_suffix = "strat"
        strat_fn = case_based_stratification
    else:
        id_suffix = "clean"
        strat_fn = stratify_by_species_and_label

    dataset = load_driams_dataset(
        driams_root,
        site,
        years=years,
        species=species,
        antibiotics=antibiotic,
        handle_missing_resistance_measurements="remove_if_all_missing",
        spectra_type="binned_6000",
        on_error="warn",
        id_suffix=id_suffix,
        extra_filters=extra_filters,
    )
    x = np.asarray([spectrum.intensities for spectrum in dataset.X])
    y = dataset.to_numpy(antibiotic)
    codes = code_series_from_dataset(dataset)
    print(
        f"[data] Loaded focal dataset site={site} drug={antibiotic}: "
        f"n={len(y)} pos={int(np.sum(y == 1))} neg={int(np.sum(y == 0))} "
        f"x_shape={x.shape}",
        flush=True,
    )
    return dataset, x, y, codes, strat_fn, dataset.y


def load_background_label_map(
    *,
    driams_root: str,
    site: str,
    years: str | Iterable[str],
    species: str,
    drugs: list[str],
) -> dict[str, dict[str, int | None]]:
    from maldi_learn.driams import load_driams_dataset
    from maldi_learn.filters import DRIAMSBooleanExpressionFilter

    print(
        f"[data] Loading background label map site={site} years={years} "
        f"species={species} drugs={','.join(drugs)}",
        flush=True,
    )

    extra_filters = []
    id_suffix = "clean"
    if site == "DRIAMS-A":
        extra_filters.append(DRIAMSBooleanExpressionFilter("workstation != HospitalHygiene"))
        id_suffix = "strat"

    dataset = load_driams_dataset(
        driams_root,
        site,
        years=years,
        species=species,
        antibiotics=drugs,
        handle_missing_resistance_measurements="remove_if_all_missing",
        spectra_type="binned_6000",
        on_error="warn",
        id_suffix=id_suffix,
        extra_filters=extra_filters,
    )
    codes = code_series_from_dataset(dataset)
    by_code: dict[str, dict[str, int | None]] = {}
    for idx, code in enumerate(codes):
        labels: dict[str, int | None] = {}
        for drug in drugs:
            if drug in dataset.y.columns:
                labels[drug] = normalize_ast_value(dataset.y.iloc[idx][drug])
            else:
                try:
                    labels[drug] = int(dataset.to_numpy(drug)[idx])
                except Exception:
                    labels[drug] = None
        by_code[str(code)] = labels
    print(f"[data] Loaded background label map site={site}: n_codes={len(by_code)}", flush=True)
    return by_code


def background_signature(labels: dict[str, int | None], focal_drug: str, drugs: list[str]) -> str:
    parts = []
    for drug in drugs:
        if drug == focal_drug:
            continue
        label = labels.get(drug)
        parts.append(f"{drug}={LABEL_CHAR.get(label, UNKNOWN_CHAR)}")
    return "|".join(parts) if parts else "NO_BACKGROUND_DRUGS"


def class_one_probability(estimator, x_test) -> np.ndarray:
    proba = estimator.predict_proba(x_test)
    classes = list(getattr(estimator, "classes_", []))
    if 1 in classes:
        return proba[:, classes.index(1)]
    if proba.shape[1] == 2:
        return proba[:, 1]
    raise ValueError(f"Cannot find resistant class probability from classes={classes}")


def metrics_from_probability(y_true, prob: np.ndarray) -> dict:
    pred = (prob >= 0.5).astype(int)
    result = {
        "accuracy": float(accuracy_score(y_true, pred)),
        "auprc": float(average_precision_score(y_true, prob)),
    }
    try:
        result["auroc"] = float(roc_auc_score(y_true, prob))
    except ValueError:
        result["auroc"] = None
    return result


def to_positional_indices(meta, indexer, *, name: str) -> np.ndarray:
    """Convert maldi-learn split output to positional indices for numpy arrays."""
    arr = np.asarray(indexer)
    if arr.dtype == bool:
        if arr.shape[0] != len(meta):
            raise ValueError(f"{name} boolean mask length {arr.shape[0]} != metadata length {len(meta)}")
        return np.flatnonzero(arr)
    if np.issubdtype(arr.dtype, np.integer):
        return arr.astype(int)

    positions = meta.index.get_indexer(indexer)
    if np.any(positions < 0):
        meta_index_as_text = meta.index.astype(str)
        lookup = {value: pos for pos, value in enumerate(meta_index_as_text)}
        positions = np.asarray([lookup.get(str(value), -1) for value in indexer], dtype=int)
    if np.any(positions < 0):
        missing = [str(value) for value, pos in zip(indexer, positions) if pos < 0][:5]
        raise ValueError(f"Could not map {name} split labels to positional indices; examples={missing}")
    return positions.astype(int)


def write_csv(path: pathlib.Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


PREDICTION_FIELDS = [
    "isolate_id",
    "site",
    "year",
    "organism",
    "drug",
    "label",
    "prob",
    "background_signature",
    "model_name",
]


def append_progress(output_dir: pathlib.Path, message: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().isoformat(timespec="seconds")
    line = f"{stamp} {message}"
    print(line, flush=True)
    with (output_dir / "progress.log").open("a") as handle:
        handle.write(line + "\n")


def write_partial_outputs(output_dir: pathlib.Path, prediction_rows: list[dict], raw_results: list[dict]) -> None:
    write_csv(output_dir / "partial_weis_predictions_long.csv", prediction_rows, PREDICTION_FIELDS)
    (output_dir / "partial_weis_raw_results.json").write_text(json.dumps(raw_results, indent=2) + "\n")


def run(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "RUN_STARTED.txt").write_text(
        f"Started {datetime.now().isoformat(timespec='seconds')}\n"
        f"model={args.model}\n"
        f"train_site={args.train_site}\n"
        f"test_sites={args.test_sites}\n"
        f"species={args.species}\n"
        f"drugs={args.drugs}\n"
        f"seed={args.seed}\n"
        f"n_folds={args.n_folds}\n"
    )
    append_progress(output_dir, "[start] Weis prediction export for background audit")
    append_progress(
        output_dir,
        f"[start] model={args.model} train_site={args.train_site} "
        f"test_sites={args.test_sites} species={args.species} "
        f"seed={args.seed} n_folds={args.n_folds}",
    )
    models = load_weis_modules(args.weis_repo)
    os.environ["DRIAMS_ROOT"] = str(args.driams_root)

    from maldi_learn.driams import DRIAMSDatasetExplorer

    drugs = split_csv_arg(args.drugs)
    test_sites = split_csv_arg(args.test_sites)

    explorer = DRIAMSDatasetExplorer(str(args.driams_root))
    train_years = explorer.available_years(args.train_site)
    append_progress(output_dir, f"[setup] train_years for {args.train_site}: {train_years}")
    raw_results = []
    prediction_rows: list[dict] = []
    trained_by_drug: dict[str, dict] = {}

    for test_site in test_sites:
        append_progress(output_dir, f"[site] Starting test_site={test_site}")
        test_years = explorer.available_years(test_site)
        append_progress(output_dir, f"[site] test_years for {test_site}: {test_years}")
        background_map = load_background_label_map(
            driams_root=str(args.driams_root),
            site=test_site,
            years=test_years,
            species=args.species,
            drugs=drugs,
        )

        for drug in drugs:
            append_progress(output_dir, f"[pair] Starting {args.species} / {drug} -> {test_site}")
            _, x_train, y_train, train_codes, train_strat_fn, train_meta = load_dataset_for_focal(
                driams_root=str(args.driams_root),
                site=args.train_site,
                years=train_years,
                species=args.species,
                antibiotic=drug,
            )
            _, x_test_full, y_test_full, test_codes, test_strat_fn, test_meta = load_dataset_for_focal(
                driams_root=str(args.driams_root),
                site=test_site,
                years=test_years,
                species=args.species,
                antibiotic=drug,
            )

            train_idx, _ = train_strat_fn(train_meta, antibiotic=drug, random_state=args.seed)
            _, test_idx = test_strat_fn(test_meta, antibiotic=drug, random_state=args.seed)
            train_idx = to_positional_indices(train_meta, train_idx, name=f"{args.train_site}/{drug}/train")
            test_idx = to_positional_indices(test_meta, test_idx, name=f"{test_site}/{drug}/test")
            append_progress(
                output_dir,
                f"[split] {drug} -> {test_site}: train_n={len(train_idx)} test_n={len(test_idx)}",
            )

            x_train_used = x_train[train_idx]
            y_train_used = y_train[train_idx]
            x_test = x_test_full[test_idx]
            y_test = y_test_full[test_idx]
            codes_test = [test_codes[idx] for idx in test_idx]

            if drug not in trained_by_drug:
                print(
                    f"Training Weis-style {args.model} once for {args.species} / {drug} "
                    f"on {args.train_site}; n_train={len(y_train_used)} n_folds={args.n_folds}",
                    flush=True,
                )
                results, estimator = models.run_experiment(
                    x_train_used,
                    y_train_used,
                    x_test,
                    y_test,
                    args.model,
                    args.n_folds,
                    random_state=args.seed,
                    verbose=True,
                    return_best_estimator=True,
                )
                trained_by_drug[drug] = {
                    "estimator": estimator,
                    "best_params": results.get("best_params"),
                    "n_train": int(len(y_train_used)),
                }
                append_progress(
                    output_dir,
                    f"[train] Finished {drug}; best_params={trained_by_drug[drug]['best_params']}",
                )
            else:
                append_progress(output_dir, f"[train] Reusing trained {args.model} model for {drug}")
                estimator = trained_by_drug[drug]["estimator"]
                results = metrics_from_probability(y_test, class_one_probability(estimator, x_test))

            probs = class_one_probability(estimator, x_test)
            append_progress(
                output_dir,
                f"[eval] {drug} -> {test_site}: auroc={results.get('auroc')} "
                f"auprc={results.get('auprc')} accuracy={results.get('accuracy')}",
            )

            raw_results.append(
                {
                    "train_site": args.train_site,
                    "test_site": test_site,
                    "train_years": train_years,
                    "test_years": test_years,
                    "species": args.species,
                    "drug": drug,
                    "model": args.model,
                    "seed": args.seed,
                    "n_train": int(trained_by_drug[drug]["n_train"]),
                    "n_test": int(len(y_test)),
                    "auroc": results.get("auroc"),
                    "auprc": results.get("auprc"),
                    "accuracy": results.get("accuracy"),
                    "best_params": trained_by_drug[drug].get("best_params"),
                }
            )

            year = "_".join(str(y) for y in test_years)
            before_rows = len(prediction_rows)
            for code, label, prob in zip(codes_test, y_test, probs):
                labels = background_map.get(str(code), {})
                prediction_rows.append(
                    {
                        "isolate_id": str(code),
                        "site": test_site,
                        "year": year,
                        "organism": args.species,
                        "drug": drug,
                        "label": int(label),
                        "prob": float(prob),
                        "background_signature": background_signature(labels, drug, drugs),
                        "model_name": f"Weis-{args.model}",
                    }
                )
            append_progress(
                output_dir,
                f"[rows] Added {len(prediction_rows) - before_rows} prediction rows for "
                f"{drug} -> {test_site}; total_rows={len(prediction_rows)}",
            )
            write_partial_outputs(output_dir, prediction_rows, raw_results)
            append_progress(output_dir, f"[write] Updated partial outputs in {output_dir}")

    predictions_csv = output_dir / "weis_predictions_long.csv"
    write_csv(
        predictions_csv,
        prediction_rows,
        PREDICTION_FIELDS,
    )
    (output_dir / "weis_raw_results.json").write_text(json.dumps(raw_results, indent=2) + "\n")
    print(f"[write] Wrote prediction table to {predictions_csv} rows={len(prediction_rows)}", flush=True)
    print(f"[write] Wrote raw Weis metrics to {output_dir / 'weis_raw_results.json'}", flush=True)

    audit_dir = output_dir / "audit"
    cmd = [
        sys.executable,
        str(args.audit_script),
        "--predictions-csv",
        str(predictions_csv),
        "--background-signature-col",
        "background_signature",
        "--model-name",
        f"Weis-{args.model}",
        "--output-dir",
        str(audit_dir),
        "--bootstrap-n",
        str(args.bootstrap_n),
        "--permutation-n",
        str(args.permutation_n),
    ]
    if args.match_year:
        cmd.append("--match-year")
    print("[audit] Running background audit command:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    print(f"[done] Wrote Weis prediction table to {predictions_csv}", flush=True)
    print(f"[done] Wrote Weis background audit to {audit_dir}", flush=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weis-repo", type=pathlib.Path, required=True)
    parser.add_argument("--driams-root", type=pathlib.Path, required=True)
    parser.add_argument("--audit-script", type=pathlib.Path, default=pathlib.Path("run_background_audit_framework.py"))
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--train-site", default="DRIAMS-A")
    parser.add_argument("--test-sites", default="DRIAMS-B,DRIAMS-C,DRIAMS-D")
    parser.add_argument("--species", default="Escherichia coli")
    parser.add_argument(
        "--drugs",
        default="Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,Ceftriaxone,Ceftazidime,Cefepime",
    )
    parser.add_argument("--model", default="lr")
    parser.add_argument("--seed", type=int, default=35)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--bootstrap-n", type=int, default=500)
    parser.add_argument("--permutation-n", type=int, default=500)
    parser.add_argument("--match-year", action="store_true")
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
