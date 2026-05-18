#!/usr/bin/env python3
"""Export LightGBM single-task and multi-task predictions for background audit."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MEGA = ROOT / "Mega_Model.py"

FIELDS = [
    "model_name",
    "site",
    "year",
    "isolate_id",
    "organism",
    "drug",
    "label",
    "prob",
    "background_signature",
    "background_known_count",
    "background_resistant_count",
]

LABEL_CHAR = {0: "S", 1: "R"}
UNKNOWN_CHAR = "U"


def load_mega(path: Path):
    spec = importlib.util.spec_from_file_location("Mega_Model_lgbm_export", path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load Mega_Model module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["Mega_Model_lgbm_export"] = module
    spec.loader.exec_module(module)
    return module


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Export LGBM predictions in long audit format.")
    parser.add_argument("--mega-model", type=Path, default=DEFAULT_MEGA)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--pair-profile", default="ecoli_mechanism6")
    parser.add_argument("--run-dir", type=Path, default=None, help="Optional run dir with config.json/models.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--variants", default="single,multi", help="Comma-separated: single,multi")
    parser.add_argument("--train-if-missing", action="store_true")
    return parser


def sample_to_prediction_row(
    sample: tuple[str, float, int],
    *,
    prob: float,
    active_lookup: dict[int, tuple[str, str]],
    site: str,
    model_name: str,
) -> dict:
    path_text, label, task_id = sample
    path = Path(path_text)
    organism, drug = active_lookup[int(task_id)]
    return {
        "model_name": model_name,
        "site": site,
        "year": path.parent.name,
        "isolate_id": path.stem,
        "organism": organism,
        "drug": drug,
        "label": int(float(label)),
        "prob": float(prob),
    }


def isolate_key(row: dict) -> tuple[str, str, str, str]:
    return (str(row["site"]), str(row["year"]), str(row["organism"]), str(row["isolate_id"]))


def add_background_signatures(rows: Sequence[dict], active_lookup: dict[int, tuple[str, str]]) -> list[dict]:
    drugs_by_organism: dict[str, list[str]] = defaultdict(list)
    for _, (organism, drug) in sorted(active_lookup.items()):
        if drug not in drugs_by_organism[organism]:
            drugs_by_organism[organism].append(drug)

    labels: dict[tuple[str, str, str, str], dict[str, int]] = defaultdict(dict)
    for row in rows:
        labels[isolate_key(row)][str(row["drug"])] = int(row["label"])

    enriched = []
    for row in rows:
        organism = str(row["organism"])
        focal = str(row["drug"])
        parts = []
        known = 0
        resistant = 0
        isolate_labels = labels[isolate_key(row)]
        for drug in drugs_by_organism.get(organism, []):
            if drug == focal:
                continue
            label = isolate_labels.get(drug)
            char = LABEL_CHAR.get(label, UNKNOWN_CHAR)
            parts.append(f"{drug}={char}")
            if label in (0, 1):
                known += 1
                resistant += int(label == 1)
        new_row = dict(row)
        new_row["background_signature"] = "|".join(parts) if parts else "NO_BACKGROUND_DRUGS"
        new_row["background_known_count"] = known
        new_row["background_resistant_count"] = resistant
        enriched.append(new_row)
    return enriched


def load_valid_spectrum_matrix(mega, samples: Sequence[tuple[str, float, int]]):
    import numpy as np

    valid_samples = []
    spectra = []
    labels = []
    task_ids = []
    for sample in samples:
        path, label, task_id = sample
        try:
            spectra.append(mega.load_spectrum(path))
        except Exception as exc:
            print(f"WARNING: skipping {path}: {exc}", flush=True)
            continue
        valid_samples.append(sample)
        labels.append(float(label))
        task_ids.append(int(task_id))
    if not valid_samples:
        return [], None, None, None
    return (
        valid_samples,
        np.stack(spectra),
        np.array(labels, dtype=np.float32),
        np.array(task_ids, dtype=np.int32),
    )


def predict_rows_for_site(
    mega,
    *,
    samples: Sequence[tuple[str, float, int]],
    site: str,
    active_lookup: dict[int, tuple[str, str]],
    variant: str,
    model,
) -> list[dict]:
    import numpy as np

    valid_samples, X, _, task_ids = load_valid_spectrum_matrix(mega, samples)
    if not valid_samples:
        return []

    probs = np.full(len(valid_samples), np.nan, dtype=np.float64)
    model_name = "LGBM-single" if variant == "single" else "LGBM-multi"
    if variant == "single":
        for task_id, task_model in model.items():
            mask = task_ids == int(task_id)
            if mask.any():
                probs[mask] = task_model.predict(X[mask])
    elif variant == "multi":
        X_aug = np.hstack([X, task_ids.reshape(-1, 1).astype(np.float32)])
        probs[:] = model.predict(X_aug)
    else:
        raise ValueError(f"Unknown variant: {variant}")

    rows = []
    for sample, prob in zip(valid_samples, probs):
        if np.isfinite(prob):
            rows.append(
                sample_to_prediction_row(
                    sample,
                    prob=float(prob),
                    active_lookup=active_lookup,
                    site=site,
                    model_name=model_name,
                )
            )
    return rows


def load_config(run_dir: Path | None) -> dict:
    if not run_dir:
        return {}
    config_path = run_dir / "config.json"
    return json.loads(config_path.read_text()) if config_path.exists() else {}


def load_lgbm_models(models_dir: Path, variants: set[str]):
    models = {}
    if "single" in variants:
        single = {}
        for path in sorted(models_dir.glob("lgbm_single_org*.pkl")):
            task_text = path.stem.replace("lgbm_single_org", "")
            if task_text.isdigit():
                with path.open("rb") as f:
                    single[int(task_text)] = pickle.load(f)
        if single:
            models["single"] = single
    if "multi" in variants:
        path = models_dir / "lgbm_multi.pkl"
        if path.exists():
            with path.open("rb") as f:
                models["multi"] = pickle.load(f)
    return models


def train_lgbm_models(mega, data_root: Path, active_pairs, variants: set[str]):
    all_source = mega.load_all_organisms(str(data_root), mega.TRAIN_SITE, active_pairs)
    train_s, val_s, _ = mega.make_split(all_source, active_pairs)
    models = {}
    if "single" in variants:
        models["single"] = mega.train_lgbm_singletask(train_s, val_s, active_pairs)
    if "multi" in variants:
        models["multi"] = mega.train_lgbm_multitask(train_s, val_s)
    return models


def build_eval_sites(mega, data_root: Path, active_pairs) -> dict[str, list[tuple[str, float, int]]]:
    all_source = mega.load_all_organisms(str(data_root), mega.TRAIN_SITE, active_pairs)
    _, _, test_s = mega.make_split(all_source, active_pairs)
    eval_sites = {}
    if test_s:
        eval_sites[f"A-{mega.TEST_YEAR}"] = test_s
    for site in mega.TEST_SITES:
        if not (data_root / site).exists():
            continue
        site_rows = mega.load_all_organisms(str(data_root), site, active_pairs)
        samples = [(p, l, o) for p, l, o, _ in site_rows]
        if samples:
            eval_sites[site] = samples
    return eval_sites


def write_rows(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})


def main() -> None:
    args = build_parser().parse_args()
    variants = {item.strip() for item in args.variants.split(",") if item.strip()}
    unknown = variants - {"single", "multi"}
    if unknown:
        raise ValueError(f"Unknown variants: {', '.join(sorted(unknown))}")

    mega = load_mega(args.mega_model)
    config = load_config(args.run_dir)
    pair_profile = config.get("pair_profile", args.pair_profile)
    mega.init_config(pair_profile)
    active_pairs = [tuple(row) for row in config.get("active_pairs", [])] or mega.screen_active_pairs(
        str(args.data_root), pair_profile
    )
    active_lookup = {int(task_id): (organism, drug) for task_id, organism, drug in active_pairs}

    models = {}
    if args.run_dir:
        models = load_lgbm_models(args.run_dir / "models", variants)
    missing = variants - set(models)
    if missing:
        if not args.train_if_missing:
            raise FileNotFoundError(
                "Missing LGBM model files for variants "
                f"{', '.join(sorted(missing))}. Re-run with --train-if-missing."
            )
        models.update(train_lgbm_models(mega, args.data_root, active_pairs, missing))

    eval_sites = build_eval_sites(mega, args.data_root, active_pairs)
    for variant in sorted(variants):
        if not models.get(variant):
            print(f"Skipping {variant}: no model available", flush=True)
            continue
        rows = []
        for site, samples in eval_sites.items():
            print(f"Predicting {variant} {site}: n={len(samples)}", flush=True)
            rows.extend(
                predict_rows_for_site(
                    mega,
                    samples=samples,
                    site=site,
                    active_lookup=active_lookup,
                    variant=variant,
                    model=models[variant],
                )
            )
        rows = add_background_signatures(rows, active_lookup)
        output = args.output_dir / f"lgbm_{variant}_predictions_long.csv"
        write_rows(output, rows)
        print(f"Wrote {output} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
