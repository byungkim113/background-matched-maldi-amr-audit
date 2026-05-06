#!/usr/bin/env python3
"""Export Weis/Borgwardt MALDI-AMR predictions for background-matched auditing.

The Weis repository's stored JSON result files contain scores and labels, but
not isolate identifiers. Background-matched auditing needs isolate identifiers
and co-resistance background labels, so this script reruns the upstream Weis
model code and writes an ID-preserving long prediction table.

Typical Kaggle usage:

    python export_weis_predictions_for_audit.py \
      --weis-repo /kaggle/working/maldi_amr \
      --driams-root /kaggle/input/datasets/drscarlat/driams \
      --panel weis-core \
      --model lightgbm \
      --external-row-policy all \
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
import time
from datetime import datetime
from typing import Iterable

import numpy as np
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score


LABEL_CHAR = {0: "S", 1: "R"}
UNKNOWN_CHAR = "U"
WEIS_SOURCE_URL = "https://github.com/BorgwardtLab/maldi_amr"
WEIS_CORE_PANELS = {
    "Escherichia coli": [
        "Ciprofloxacin",
        "Ceftriaxone",
        "Cefepime",
        "Piperacillin-Tazobactam",
        "Tobramycin",
    ],
    "Staphylococcus aureus": [
        "Oxacillin",
        "Penicillin",
        "Ciprofloxacin",
        "Fusidic acid",
    ],
    "Klebsiella pneumoniae": [
        "Ciprofloxacin",
        "Ceftriaxone",
        "Cefepime",
        "Amoxicillin-Clavulanic acid",
        "Meropenem",
        "Tobramycin",
        "Piperacillin-Tazobactam",
    ],
}


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


def resolve_species_drug_panels(args: argparse.Namespace) -> list[tuple[str, list[str]]]:
    if args.panel == "custom":
        return [(args.species, split_csv_arg(args.drugs))]
    if args.panel == "weis-core":
        return [(species, list(drugs)) for species, drugs in WEIS_CORE_PANELS.items()]
    raise ValueError(f"Unknown panel: {args.panel}")


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


def select_external_indices(
    meta,
    split_fn,
    *,
    drug: str,
    seed: int,
    policy: str,
    name: str,
) -> np.ndarray:
    """Select external rows for scoring.

    ``all`` is the audit default because background matching needs isolate-level
    coverage across the full external site. ``stratified`` preserves the older
    Weis-style subset behavior for diagnostics and metric comparisons.
    """
    if policy == "all":
        return np.arange(len(meta), dtype=int)
    if policy == "stratified":
        _, test_idx = split_fn(meta, antibiotic=drug, random_state=seed)
        return to_positional_indices(meta, test_idx, name=name)
    raise ValueError(f"Unknown external row policy: {policy}")


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


def planned_pair_count(test_sites: Iterable[str], species_panels: list[tuple[str, list[str]]]) -> int:
    return len(list(test_sites)) * sum(len(drugs) for _, drugs in species_panels)


def pair_progress_percent(pair_index: int, total_pairs: int) -> float | None:
    if total_pairs <= 0:
        return None
    return round((pair_index / total_pairs) * 100.0, 1)


def format_pair_progress(
    pair_index: int,
    total_pairs: int,
    test_site: str,
    species: str,
    drug: str,
) -> str:
    percent = pair_progress_percent(pair_index, total_pairs)
    if percent is None:
        return f"[pair {pair_index}/?] {test_site} | {species} / {drug}"
    return f"[pair {pair_index}/{total_pairs} {percent:.1f}%] {test_site} | {species} / {drug}"


def format_seconds(seconds: float) -> str:
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {secs}s"
    if minutes:
        return f"{minutes}m {secs}s"
    return f"{secs}s"


def write_status(
    output_dir: pathlib.Path,
    *,
    stage: str,
    pair_index: int,
    total_pairs: int,
    test_site: str | None = None,
    species: str | None = None,
    drug: str | None = None,
    rows_written: int = 0,
    elapsed_seconds: float | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    percent = pair_progress_percent(pair_index, total_pairs)
    status = {
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "stage": stage,
        "pair_index": pair_index,
        "total_pairs": total_pairs,
        "percent_pairs_complete": percent,
        "test_site": test_site,
        "species": species,
        "drug": drug,
        "rows_written": rows_written,
    }
    if elapsed_seconds is not None:
        status["elapsed"] = format_seconds(elapsed_seconds)
        status["elapsed_seconds"] = round(float(elapsed_seconds), 1)
    (output_dir / "STATUS.json").write_text(json.dumps(status, indent=2) + "\n")

    percent_text = "unknown" if percent is None else f"{percent:.1f}%"
    current = (
        f"{status['updated_at']} | stage={stage} | pair={pair_index}/{total_pairs} "
        f"({percent_text}) | site={test_site or '-'} | species={species or '-'} | "
        f"drug={drug or '-'} | rows={rows_written}"
    )
    if elapsed_seconds is not None:
        current += f" | elapsed={format_seconds(elapsed_seconds)}"
    (output_dir / "CURRENT_STAGE.txt").write_text(current + "\n")


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


def write_reproduction_report(
    output_dir: pathlib.Path,
    args: argparse.Namespace,
    species_panels: list[tuple[str, list[str]]],
    raw_results: list[dict],
    predictions_csv: pathlib.Path,
) -> None:
    report_json = {
        "source_url": WEIS_SOURCE_URL,
        "source_repo": str(args.weis_repo),
        "model_code": str(args.weis_repo / "amr_maldi_ml" / "models.py"),
        "panel": args.panel,
        "external_row_policy": args.external_row_policy,
        "model": args.model,
        "seed": args.seed,
        "n_folds": args.n_folds,
        "train_site": args.train_site,
        "test_sites": split_csv_arg(args.test_sites),
        "species_panels": [{"species": species, "drugs": drugs} for species, drugs in species_panels],
        "prediction_csv": str(predictions_csv),
        "n_prediction_rows": sum(int(row.get("n_test", 0)) for row in raw_results),
        "n_raw_result_rows": len(raw_results),
        "note": (
            "This is an ID-preserving rerun through the original BorgwardtLab/maldi_amr "
            "model code. Exact paper parity should additionally be checked by comparing "
            "raw metrics against upstream Weis result JSONs for the same model, seed, "
            "train/test split, and preprocessing."
        ),
    }
    (output_dir / "weis_reproduction_report.json").write_text(json.dumps(report_json, indent=2) + "\n")

    panel_lines = "\n".join(
        f"- {species}: {', '.join(drugs)}"
        for species, drugs in species_panels
    )
    markdown = f"""# Weis-Code Background Audit Export

Source repository: [{WEIS_SOURCE_URL}]({WEIS_SOURCE_URL})

Loaded model code: `{args.weis_repo / "amr_maldi_ml" / "models.py"}`

This export reruns the original Weis/Borgwardt model implementation and writes
isolate-level prediction rows for the background-matched audit.

## Configuration

- Model: `{args.model}`
- Panel: `{args.panel}`
- External row policy: `{args.external_row_policy}`
- Train site: `{args.train_site}`
- Test sites: `{args.test_sites}`
- Seed: `{args.seed}`
- Folds: `{args.n_folds}`
- Prediction CSV: `{predictions_csv}`

## Species and Drugs

{panel_lines}

## Parity Note

This should be described as a Weis-code rerun with ID-preserving all-row external
scoring. Calling it an exact reproduction of the published paper additionally
requires matching the upstream stored raw metrics under the same model, seed,
split logic, and preprocessing.
"""
    (output_dir / "weis_reproduction_report.md").write_text(markdown)


def run(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    started_at = time.monotonic()
    (output_dir / "RUN_STARTED.txt").write_text(
        f"Started {datetime.now().isoformat(timespec='seconds')}\n"
        f"model={args.model}\n"
        f"panel={args.panel}\n"
        f"external_row_policy={args.external_row_policy}\n"
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
        f"test_sites={args.test_sites} panel={args.panel} "
        f"external_row_policy={args.external_row_policy} "
        f"seed={args.seed} n_folds={args.n_folds}",
    )
    models = load_weis_modules(args.weis_repo)
    os.environ["DRIAMS_ROOT"] = str(args.driams_root)

    from maldi_learn.driams import DRIAMSDatasetExplorer

    species_panels = resolve_species_drug_panels(args)
    test_sites = split_csv_arg(args.test_sites)
    total_pairs = planned_pair_count(test_sites, species_panels)
    pair_index = 0
    write_status(
        output_dir,
        stage="setup",
        pair_index=pair_index,
        total_pairs=total_pairs,
        rows_written=0,
        elapsed_seconds=time.monotonic() - started_at,
    )
    append_progress(
        output_dir,
        f"[plan] total_pair_site_jobs={total_pairs} "
        f"test_sites={len(test_sites)} organism_panels={len(species_panels)}",
    )

    explorer = DRIAMSDatasetExplorer(str(args.driams_root))
    train_years = explorer.available_years(args.train_site)
    append_progress(output_dir, f"[setup] train_years for {args.train_site}: {train_years}")
    raw_results = []
    prediction_rows: list[dict] = []
    trained_by_pair: dict[tuple[str, str], dict] = {}

    for test_site in test_sites:
        append_progress(output_dir, f"[site] Starting test_site={test_site}")
        test_years = explorer.available_years(test_site)
        append_progress(output_dir, f"[site] test_years for {test_site}: {test_years}")

        for species, drugs in species_panels:
            append_progress(output_dir, f"[species] Starting {species} on {test_site}")
            background_map = load_background_label_map(
                driams_root=str(args.driams_root),
                site=test_site,
                years=test_years,
                species=species,
                drugs=drugs,
            )

            for drug in drugs:
                pair_index += 1
                pair_label = format_pair_progress(pair_index, total_pairs, test_site, species, drug)
                write_status(
                    output_dir,
                    stage="loading_pair",
                    pair_index=pair_index,
                    total_pairs=total_pairs,
                    test_site=test_site,
                    species=species,
                    drug=drug,
                    rows_written=len(prediction_rows),
                    elapsed_seconds=time.monotonic() - started_at,
                )
                append_progress(output_dir, f"{pair_label} [stage loading_pair] starting")
                pair_key = (species, drug)
                _, x_train, y_train, train_codes, train_strat_fn, train_meta = load_dataset_for_focal(
                    driams_root=str(args.driams_root),
                    site=args.train_site,
                    years=train_years,
                    species=species,
                    antibiotic=drug,
                )
                _, x_test_full, y_test_full, test_codes, test_strat_fn, test_meta = load_dataset_for_focal(
                    driams_root=str(args.driams_root),
                    site=test_site,
                    years=test_years,
                    species=species,
                    antibiotic=drug,
                )

                train_idx, _ = train_strat_fn(train_meta, antibiotic=drug, random_state=args.seed)
                train_idx = to_positional_indices(
                    train_meta,
                    train_idx,
                    name=f"{args.train_site}/{species}/{drug}/train",
                )
                test_idx = select_external_indices(
                    test_meta,
                    test_strat_fn,
                    drug=drug,
                    seed=args.seed,
                    policy=args.external_row_policy,
                    name=f"{test_site}/{species}/{drug}/test",
                )
                append_progress(
                    output_dir,
                    f"{pair_label} [stage split] "
                    f"train_n={len(train_idx)} test_n={len(test_idx)} "
                    f"external_row_policy={args.external_row_policy}",
                )
                write_status(
                    output_dir,
                    stage="split_done",
                    pair_index=pair_index,
                    total_pairs=total_pairs,
                    test_site=test_site,
                    species=species,
                    drug=drug,
                    rows_written=len(prediction_rows),
                    elapsed_seconds=time.monotonic() - started_at,
                )

                x_train_used = x_train[train_idx]
                y_train_used = y_train[train_idx]
                x_test = x_test_full[test_idx]
                y_test = y_test_full[test_idx]
                codes_test = [test_codes[idx] for idx in test_idx]

                if pair_key not in trained_by_pair:
                    write_status(
                        output_dir,
                        stage="training_model",
                        pair_index=pair_index,
                        total_pairs=total_pairs,
                        test_site=test_site,
                        species=species,
                        drug=drug,
                        rows_written=len(prediction_rows),
                        elapsed_seconds=time.monotonic() - started_at,
                    )
                    append_progress(
                        output_dir,
                        f"{pair_label} [stage training_model] "
                        f"model={args.model} train_n={len(y_train_used)} n_folds={args.n_folds}",
                    )
                    print(
                        f"Training Weis-code {args.model} once for {species} / {drug} "
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
                    trained_by_pair[pair_key] = {
                        "estimator": estimator,
                        "best_params": results.get("best_params"),
                        "n_train": int(len(y_train_used)),
                    }
                    append_progress(
                        output_dir,
                        f"{pair_label} [stage training_done] "
                        f"best_params={trained_by_pair[pair_key]['best_params']}",
                    )
                else:
                    append_progress(output_dir, f"{pair_label} [stage training_reuse] Reusing trained {args.model} model")
                    estimator = trained_by_pair[pair_key]["estimator"]
                    results = metrics_from_probability(y_test, class_one_probability(estimator, x_test))

                write_status(
                    output_dir,
                    stage="scoring",
                    pair_index=pair_index,
                    total_pairs=total_pairs,
                    test_site=test_site,
                    species=species,
                    drug=drug,
                    rows_written=len(prediction_rows),
                    elapsed_seconds=time.monotonic() - started_at,
                )
                probs = class_one_probability(estimator, x_test)
                append_progress(
                    output_dir,
                    f"{pair_label} [stage eval] auroc={results.get('auroc')} "
                    f"auprc={results.get('auprc')} accuracy={results.get('accuracy')}",
                )

                raw_results.append(
                    {
                        "source_url": WEIS_SOURCE_URL,
                        "panel": args.panel,
                        "external_row_policy": args.external_row_policy,
                        "train_site": args.train_site,
                        "test_site": test_site,
                        "train_years": train_years,
                        "test_years": test_years,
                        "species": species,
                        "drug": drug,
                        "model": args.model,
                        "seed": args.seed,
                        "n_train": int(trained_by_pair[pair_key]["n_train"]),
                        "n_test": int(len(y_test)),
                        "auroc": results.get("auroc"),
                        "auprc": results.get("auprc"),
                        "accuracy": results.get("accuracy"),
                        "best_params": trained_by_pair[pair_key].get("best_params"),
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
                            "organism": species,
                            "drug": drug,
                            "label": int(label),
                            "prob": float(prob),
                            "background_signature": background_signature(labels, drug, drugs),
                            "model_name": f"Weis-{args.model}",
                        }
                    )
                append_progress(
                    output_dir,
                    f"{pair_label} [stage rows] Added {len(prediction_rows) - before_rows} "
                    f"prediction rows; total_rows={len(prediction_rows)}",
                )
                write_partial_outputs(output_dir, prediction_rows, raw_results)
                write_status(
                    output_dir,
                    stage="partial_outputs_written",
                    pair_index=pair_index,
                    total_pairs=total_pairs,
                    test_site=test_site,
                    species=species,
                    drug=drug,
                    rows_written=len(prediction_rows),
                    elapsed_seconds=time.monotonic() - started_at,
                )
                append_progress(output_dir, f"{pair_label} [stage write] Updated partial outputs in {output_dir}")

    predictions_csv = output_dir / "weis_predictions_long.csv"
    write_csv(
        predictions_csv,
        prediction_rows,
        PREDICTION_FIELDS,
    )
    (output_dir / "weis_raw_results.json").write_text(json.dumps(raw_results, indent=2) + "\n")
    write_reproduction_report(output_dir, args, species_panels, raw_results, predictions_csv)
    write_status(
        output_dir,
        stage="running_background_audit",
        pair_index=pair_index,
        total_pairs=total_pairs,
        rows_written=len(prediction_rows),
        elapsed_seconds=time.monotonic() - started_at,
    )
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
    append_progress(output_dir, "[audit] Running background audit command: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    write_status(
        output_dir,
        stage="done",
        pair_index=pair_index,
        total_pairs=total_pairs,
        rows_written=len(prediction_rows),
        elapsed_seconds=time.monotonic() - started_at,
    )
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
    parser.add_argument(
        "--panel",
        choices=["custom", "weis-core"],
        default="custom",
        help=(
            "custom uses --species/--drugs. weis-core uses the organism-drug "
            "pair list from BorgwardtLab/maldi_amr."
        ),
    )
    parser.add_argument("--species", default="Escherichia coli")
    parser.add_argument(
        "--drugs",
        default="Ciprofloxacin,Norfloxacin,Amoxicillin-Clavulanic acid,Ceftriaxone,Ceftazidime,Cefepime",
    )
    parser.add_argument(
        "--external-row-policy",
        choices=["all", "stratified"],
        default="all",
        help=(
            "all scores every eligible external isolate, which is preferred for "
            "background matching. stratified preserves the older split subset."
        ),
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
