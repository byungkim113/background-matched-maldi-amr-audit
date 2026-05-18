#!/usr/bin/env python3
"""MARISMa v2 preprocessing pipeline for background-matched MALDI-AMR audits.

This script is intentionally staged. The ``preprocess`` stage reads MARISMa
Bruker spectra, converts them into the 6000-bin representation used by the
DRIAMS/Mega experiments, and writes a long manifest of organism-drug labels.
The ``predict`` stage loads a completed Mega_Model run and writes isolate-level
MARISMa predictions for overlapping organism-drug pairs.

Example Kaggle use:

    python /kaggle/working/marisma_end_to_end_kaggle.py \\
      --stage preprocess \\
      --amr-csv /kaggle/input/datasets/bfdf121/marisma/AMR.csv \\
      --marisma-root /kaggle/input/datasets/bfdf121/marisma/MARISMa \\
      --output-dir /kaggle/working/marisma_preprocessed
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


TARGET_MZ_MIN = 2000.0
TARGET_MZ_MAX = 20000.0
N_BINS = 6000
TARGET_MZ = np.linspace(TARGET_MZ_MIN, TARGET_MZ_MAX, N_BINS, dtype=np.float64)


@dataclass(frozen=True)
class TargetPair:
    organism: str
    paper_drug: str
    marisma_drug: str
    relationship: str
    ecology_block: str


TARGET_PAIRS = [
    TargetPair("Escherichia coli", "Ciprofloxacin", "Ciprofloxacin", "exact", "fluoroquinolone"),
    TargetPair("Escherichia coli", "Norfloxacin", "Norfloxacin", "exact", "fluoroquinolone"),
    TargetPair("Escherichia coli", "Levofloxacin", "Levofloxacin", "external_extension", "fluoroquinolone"),
    TargetPair(
        "Escherichia coli",
        "Amoxicillin-Clavulanic acid",
        "Amoxicillin/Clavulanic acid",
        "spelling_alias",
        "beta-lactam/inhibitor",
    ),
    TargetPair(
        "Escherichia coli",
        "Ceftriaxone",
        "Cefotaxime",
        "third_generation_cephalosporin_analog",
        "third-generation cephalosporin",
    ),
    TargetPair("Escherichia coli", "Ceftazidime", "Ceftazidime", "exact", "third-generation cephalosporin"),
    TargetPair("Escherichia coli", "Cefepime", "Cefepime", "exact", "fourth-generation cephalosporin"),
    TargetPair(
        "Escherichia coli",
        "Cotrimoxazole",
        "Trimethoprim/Sulfamethoxazole",
        "name_alias",
        "folate-pathway/MDR block",
    ),
    TargetPair("Staphylococcus aureus", "Oxacillin", "Oxacillin", "exact", "MRSA/beta-lactam"),
    TargetPair("Staphylococcus epidermidis", "Erythromycin", "Erythromycin", "exact", "macrolide"),
]


BACKGROUND_PANELS = {
    "Escherichia coli": [
        "Ciprofloxacin",
        "Norfloxacin",
        "Levofloxacin",
        "Amoxicillin/Clavulanic acid",
        "Cefotaxime",
        "Ceftazidime",
        "Cefepime",
        "Trimethoprim/Sulfamethoxazole",
    ],
    "Staphylococcus aureus": [
        "Oxacillin",
        "Erythromycin",
        "Clindamycin",
        "Ciprofloxacin",
        "Levofloxacin",
        "Trimethoprim/Sulfamethoxazole",
        "Penicillin",
    ],
    "Staphylococcus epidermidis": [
        "Oxacillin",
        "Erythromycin",
        "Clindamycin",
        "Ciprofloxacin",
        "Levofloxacin",
        "Trimethoprim/Sulfamethoxazole",
    ],
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Preprocess MARISMa v2 Bruker spectra and export Mega_Model "
            "prediction rows for background-matched MALDI-AMR audits."
        )
    )
    p.add_argument(
        "--stage",
        choices=["preprocess", "predict", "audit", "all"],
        default="preprocess",
        help=(
            "preprocess builds vectors/manifest; predict scores an existing Mega_Model run; "
            "audit aggregates spot-level predictions to isolate/drug rows and runs the background audit; "
            "all runs preprocess then predict."
        ),
    )
    p.add_argument("--amr-csv", type=Path, default=None, help="Path to MARISMa AMR.csv.")
    p.add_argument("--marisma-root", type=Path, default=None, help="Path to directory containing year folders.")
    p.add_argument("--output-dir", type=Path, default=Path("/kaggle/working/marisma_preprocessed"))
    p.add_argument("--min-s", type=int, default=50)
    p.add_argument("--min-r", type=int, default=50)
    p.add_argument("--max-spectra", type=int, default=None, help="Optional smoke-test cap on unique spectra.")
    p.add_argument("--progress-every", type=int, default=250)
    p.add_argument("--include-target", action="append", default=None, help="Optional paper drug names to include.")
    p.add_argument("--skip-existing", action="store_true", help="Skip vectorization if all preprocess outputs exist.")
    p.add_argument("--mega-model-path", type=Path, default=Path("Mega_Model.py"), help="Path to Mega_Model.py.")
    p.add_argument("--run-dir", type=Path, default=None, help="Mega run directory containing config.json and models/.")
    p.add_argument("--vectors-npy", type=Path, default=None, help="Override path to marisma_vectors_6000.npy.")
    p.add_argument("--manifest-csv", type=Path, default=None, help="Override path to marisma_prediction_manifest.csv.")
    p.add_argument("--prediction-csv", type=Path, default=None, help="Output path for marisma_mega_predictions_long.csv.")
    p.add_argument("--model-name", default="Mega-CNN-MARISMa")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--tta-passes", type=int, default=None, help="Override config TTA passes. Use 1 for smoke tests.")
    p.add_argument("--max-prediction-rows", type=int, default=None, help="Optional smoke-test cap after pair filtering.")
    p.add_argument(
        "--audit-script",
        type=Path,
        default=Path("run_background_audit_framework.py"),
        help="Path to the model-agnostic audit script for --stage audit.",
    )
    p.add_argument(
        "--audit-output-dir",
        type=Path,
        default=None,
        help="Audit output directory. Defaults to <output-dir>/marisma_isolate_background_audit.",
    )
    p.add_argument("--bootstrap-n", type=int, default=500, help="Bootstrap replicates for --stage audit.")
    p.add_argument("--permutation-n", type=int, default=500, help="Permutation replicates for --stage audit.")
    return p


def log(message: str) -> None:
    print(message, flush=True)


def resolve_existing_path(user_path: Path | None, candidates: Iterable[Path], label: str) -> Path:
    checked: list[Path] = []
    if user_path is not None:
        checked.append(user_path)
        if user_path.exists():
            return user_path
    for path in candidates:
        checked.append(path)
        if path.exists():
            return path
    checked_text = "\n  ".join(str(p) for p in checked)
    raise FileNotFoundError(f"{label} not found. Checked:\n  {checked_text}")


def default_amr_candidates() -> list[Path]:
    return [
        Path("/kaggle/input/datasets/bfdf121/marisma/AMR.csv"),
        Path("/kaggle/input/marisma/AMR.csv"),
        Path("/kaggle/input/datasets/AMR.csv"),
        Path("AMR.csv"),
    ]


def default_marisma_root_candidates() -> list[Path]:
    return [
        Path("/kaggle/input/datasets/bfdf121/marisma/MARISMa"),
        Path("/kaggle/input/marisma/MARISMa"),
        Path("/kaggle/input/datasets/MARISMa"),
        Path("MARISMa"),
    ]


def normalize_label(value) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text or text in {"NAN", "NA"}:
        return None
    return text


def binary_label_value(value) -> int | None:
    text = normalize_label(value)
    if text == "S":
        return 0
    if text in {"R", "R*", "ESBL"}:
        return 1
    return None


def parse_jcamp_params(text: str, required: Iterable[str]) -> dict[str, float]:
    params: dict[str, float] = {}
    for key in required:
        match = re.search(rf"##\${re.escape(key)}=\s*([\-0-9.eE]+)", text)
        if match is None:
            raise ValueError(f"Missing Bruker parameter {key}")
        params[key] = float(match.group(1))
    return params


def bruker_mz_axis(params: dict[str, float]) -> np.ndarray:
    """Reconstruct the Bruker MALDI m/z axis from flex/XMASS metadata.

    Formula follows the quadratic Bruker TOF calibration used by
    readBrukerFlexData-style readers. It was validated against MARISMa CalStar
    reference peaks with sub-Da agreement in the user's Kaggle smoke test.
    """
    delay = float(params["DELAY"])
    dw = float(params["DW"])
    ml1 = float(params["ML1"])
    ml2 = float(params["ML2"])
    ml3 = float(params["ML3"])
    td = int(round(float(params["TD"])))

    if td <= 0:
        raise ValueError(f"Invalid Bruker TD/SI length: {td}")
    if dw <= 0:
        raise ValueError(f"Invalid Bruker DW value: {dw}")

    times = delay + np.arange(td, dtype=np.float64) * dw
    b = math.sqrt(1e12 / ml1)
    a = ml3
    c = ml2 - times

    if abs(a) < 1e-12:
        root = -c / b
    else:
        discriminant = b * b - 4.0 * a * c
        if np.any(discriminant < 0):
            raise ValueError("Negative Bruker calibration discriminant")
        sqrt_disc = np.sqrt(discriminant)
        r1 = (-b + sqrt_disc) / (2.0 * a)
        r2 = (-b - sqrt_disc) / (2.0 * a)
        root = np.where((r1 > 0) & ((r1 < r2) | (r2 <= 0)), r1, r2)

    mz = root * root
    if not np.isfinite(mz).all():
        raise ValueError("Non-finite values in reconstructed m/z axis")
    if np.any(np.diff(mz) <= 0):
        raise ValueError("Reconstructed m/z axis is not strictly increasing")
    return mz


def analysis_dir_for_spot(spot_dir: Path) -> Path:
    spot_dir = Path(spot_dir)
    if (spot_dir / "pdata" / "1" / "1r").exists() and (spot_dir / "acqu").exists():
        return spot_dir
    analysis = spot_dir / "1" / "1SLin"
    if (analysis / "pdata" / "1" / "1r").exists() and (analysis / "acqu").exists():
        return analysis
    raise FileNotFoundError(f"Could not find Bruker 1/1SLin analysis under {spot_dir}")


def read_bruker_spot(spot_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    analysis = analysis_dir_for_spot(spot_dir)
    acqu_text = (analysis / "acqu").read_text(errors="ignore")
    params = parse_jcamp_params(acqu_text, ["DELAY", "DW", "ML1", "ML2", "ML3", "TD"])
    mz = bruker_mz_axis(params)
    intensity = np.fromfile(analysis / "pdata" / "1" / "1r", dtype="<i4").astype(np.float32)
    if intensity.shape[0] != mz.shape[0]:
        raise ValueError(f"m/z and intensity length mismatch: {mz.shape[0]} vs {intensity.shape[0]}")
    intensity[intensity < 0] = 0.0
    return mz, intensity


def vectorize_spectrum(mz: np.ndarray, intensity: np.ndarray) -> np.ndarray:
    if mz.shape[0] != intensity.shape[0]:
        raise ValueError("m/z and intensity arrays must have the same length")
    if mz.shape[0] < 2:
        raise ValueError("Spectrum must contain at least two points")
    y = np.asarray(intensity, dtype=np.float32).copy()
    y[y < 0] = 0.0
    vector = np.interp(TARGET_MZ, mz, y, left=0.0, right=0.0).astype(np.float32)
    vector = np.log1p(vector)
    std = float(vector.std())
    if std > 0:
        vector = (vector - float(vector.mean())) / std
    else:
        vector = vector - float(vector.mean())
    if not np.isfinite(vector).all():
        raise ValueError("Non-finite values in vectorized spectrum")
    return vector.astype(np.float32)


def resolve_spot_path(marisma_root: Path, path_value: str) -> Path:
    text = str(path_value).strip()
    parts = [part for part in text.split("/") if part]
    if parts and parts[0].lower() == "marisma":
        parts = parts[1:]
    return Path(marisma_root).joinpath(*parts)


def selected_targets(include_target: list[str] | None) -> list[TargetPair]:
    if not include_target:
        return TARGET_PAIRS
    wanted = {target.strip().lower() for target in include_target}
    targets = [target for target in TARGET_PAIRS if target.paper_drug.lower() in wanted or target.marisma_drug.lower() in wanted]
    if not targets:
        raise ValueError(f"No target pairs matched --include-target values: {include_target}")
    return targets


def count_sr_labels(df: pd.DataFrame, target: TargetPair) -> tuple[int, int]:
    labels = df.loc[df["Species"].eq(target.organism), target.marisma_drug].map(binary_label_value)
    return int((labels == 0).sum()), int((labels == 1).sum())


def build_label_manifest(df: pd.DataFrame, targets: list[TargetPair], min_s: int, min_r: int) -> pd.DataFrame:
    required = {"Identifier", "target_position", "Year", "Path", "Species", "Microorganism"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"AMR.csv is missing required columns: {missing}")

    rows: list[dict] = []
    for target in targets:
        if target.marisma_drug not in df.columns:
            log(f"[label] skip {target.organism} / {target.paper_drug}: missing column {target.marisma_drug}")
            continue
        sub = df[df["Species"].eq(target.organism)].copy()
        n_s, n_r = count_sr_labels(df, target)
        if n_s < min_s or n_r < min_r:
            log(f"[label] skip {target.organism} / {target.paper_drug}: n_s={n_s} n_r={n_r}")
            continue
        labels = sub[target.marisma_drug].map(binary_label_value)
        keep = labels.notna()
        if not keep.any():
            continue
        for idx, rec in sub.loc[keep].iterrows():
            row = {
                "source_row_index": int(idx),
                "isolate_id": rec["Identifier"],
                "spot_id": rec["target_position"],
                "site": "MARISMa",
                "year": int(rec["Year"]) if not pd.isna(rec["Year"]) else "",
                "organism": target.organism,
                "species": rec["Species"],
                "microorganism": rec["Microorganism"],
                "paper_drug": target.paper_drug,
                "marisma_drug": target.marisma_drug,
                "drug": target.paper_drug,
                "drug_relationship": target.relationship,
                "ecology_block": target.ecology_block,
                "label": int(labels.loc[idx]),
                "path": rec["Path"],
            }
            panel = BACKGROUND_PANELS.get(target.organism, [])
            for drug in panel:
                if drug in df.columns:
                    bg = binary_label_value(rec.get(drug))
                    row[f"background__{drug}"] = "" if bg is None else int(bg)
            rows.append(row)

    if not rows:
        return pd.DataFrame()
    manifest = pd.DataFrame(rows).drop_duplicates()
    return manifest.sort_values(["organism", "drug", "year", "isolate_id", "spot_id"]).reset_index(drop=True)


def markdown_table(df: pd.DataFrame, columns: list[str] | None = None, max_rows: int | None = None) -> str:
    """Render a small markdown table without pandas' optional tabulate dependency."""
    if df.empty:
        return "_No rows._"
    out = df.copy()
    if columns is not None:
        out = out[columns]
    if max_rows is not None:
        out = out.head(max_rows)
    header = "| " + " | ".join(out.columns.astype(str)) + " |"
    divider = "| " + " | ".join(["---"] * len(out.columns)) + " |"
    rows = []
    for rec in out.to_dict("records"):
        vals = []
        for col in out.columns:
            val = rec[col]
            if isinstance(val, float):
                val = f"{val:.4g}" if math.isfinite(val) else ""
            vals.append(str(val))
        rows.append("| " + " | ".join(vals) + " |")
    return "\n".join([header, divider, *rows])


def write_markdown_report(output_dir: Path, summary: dict, label_manifest: pd.DataFrame, failures: pd.DataFrame) -> None:
    lines = [
        "# MARISMa v2 Preprocessing Report",
        "",
        "## Summary",
        "",
        f"- Label rows requested: **{summary['n_label_rows_initial']:,}**",
        f"- Label rows with processed vectors: **{summary['n_label_rows_with_vectors']:,}**",
        f"- Unique spectra attempted: **{summary['n_unique_spectra_attempted']:,}**",
        f"- Unique spectra vectorized: **{summary['n_unique_spectra_vectorized']:,}**",
        f"- Vector shape: **{summary['vector_shape']}**",
        f"- Failures: **{summary['n_failures']:,}**",
        "",
        "## Outputs",
        "",
        "- `marisma_vectors_6000.npy`: model-ready 6000-bin spectra.",
        "- `marisma_prediction_manifest.csv`: long label manifest with `vector_index`.",
        "- `marisma_unique_spectra_manifest.csv`: one row per vector.",
        "- `marisma_preprocess_failures.csv`: spectra that could not be read.",
        "",
        "## Drug Mapping Note",
        "",
        "`E. coli / Ceftriaxone` is mapped to MARISMa `E. coli / Cefotaxime` as a third-generation cephalosporin analog. Keep both `paper_drug` and `marisma_drug` columns in downstream reports.",
        "",
    ]
    if not label_manifest.empty:
        counts = label_manifest.groupby(["organism", "paper_drug", "marisma_drug", "label"]).size().unstack(fill_value=0)
        counts = counts.reset_index()
        lines.extend(["## Label Counts After Vectorization", "", markdown_table(counts), ""])
    if not failures.empty:
        lines.extend(["## First Failures", "", markdown_table(failures, max_rows=20), ""])
    (output_dir / "marisma_preprocess_report.md").write_text("\n".join(lines) + "\n")


def preprocess(args: argparse.Namespace) -> None:
    amr_csv = resolve_existing_path(args.amr_csv, default_amr_candidates(), "AMR.csv")
    marisma_root = resolve_existing_path(args.marisma_root, default_marisma_root_candidates(), "MARISMa root")
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    outputs = [
        output_dir / "marisma_vectors_6000.npy",
        output_dir / "marisma_prediction_manifest.csv",
        output_dir / "marisma_unique_spectra_manifest.csv",
        output_dir / "marisma_preprocess_report.json",
    ]
    if args.skip_existing and all(path.exists() for path in outputs):
        log(f"[preprocess] outputs already exist in {output_dir}; skipping")
        return

    log(f"[preprocess] reading AMR labels: {amr_csv}")
    df = pd.read_csv(amr_csv, low_memory=False)
    targets = selected_targets(args.include_target)
    label_manifest = build_label_manifest(df, targets, args.min_s, args.min_r)
    if label_manifest.empty:
        raise ValueError("No eligible MARISMa label rows found")

    label_manifest["spot_dir"] = label_manifest["path"].map(lambda p: str(resolve_spot_path(marisma_root, p)))
    unique_spots = label_manifest[["spot_dir", "path", "isolate_id", "spot_id", "year", "organism", "species", "microorganism"]]
    unique_spots = unique_spots.drop_duplicates(subset=["spot_dir"]).sort_values(["year", "organism", "spot_dir"])
    if args.max_spectra is not None:
        unique_spots = unique_spots.head(args.max_spectra).copy()
        log(f"[preprocess] max_spectra cap active: {args.max_spectra}")

    vectors: list[np.ndarray] = []
    successes: list[dict] = []
    failures: list[dict] = []
    spot_to_index: dict[str, int] = {}

    log(f"[preprocess] vectorizing {len(unique_spots):,} unique spectra")
    for attempt_index, rec in enumerate(unique_spots.to_dict("records"), start=1):
        spot_dir = Path(rec["spot_dir"])
        try:
            mz, intensity = read_bruker_spot(spot_dir)
            vector = vectorize_spectrum(mz, intensity)
            vector_index = len(vectors)
            vectors.append(vector)
            spot_to_index[str(spot_dir)] = vector_index
            successes.append(
                {
                    **rec,
                    "vector_index": vector_index,
                    "mz_min": float(mz.min()),
                    "mz_max": float(mz.max()),
                    "raw_max": float(intensity.max()),
                    "vector_mean": float(vector.mean()),
                    "vector_std": float(vector.std()),
                }
            )
        except Exception as exc:  # noqa: BLE001 - preserve per-spectrum errors in CSV
            failures.append({**rec, "error": repr(exc)})
        if attempt_index == 1 or attempt_index % args.progress_every == 0:
            log(f"[preprocess] {attempt_index:,}/{len(unique_spots):,} attempted; ok={len(vectors):,} fail={len(failures):,}")

    if not vectors:
        raise ValueError("No MARISMa spectra were successfully vectorized")

    vector_array = np.stack(vectors).astype(np.float32)
    label_manifest["vector_index"] = label_manifest["spot_dir"].map(spot_to_index)
    label_manifest = label_manifest[label_manifest["vector_index"].notna()].copy()
    label_manifest["vector_index"] = label_manifest["vector_index"].astype(int)

    success_df = pd.DataFrame(successes).sort_values("vector_index")
    failure_df = pd.DataFrame(failures)

    np.save(output_dir / "marisma_vectors_6000.npy", vector_array)
    label_manifest.to_csv(output_dir / "marisma_prediction_manifest.csv", index=False)
    success_df.to_csv(output_dir / "marisma_unique_spectra_manifest.csv", index=False)
    failure_df.to_csv(output_dir / "marisma_preprocess_failures.csv", index=False)

    summary = {
        "amr_csv": str(amr_csv),
        "marisma_root": str(marisma_root),
        "n_amr_rows": int(len(df)),
        "n_label_rows_initial": int(len(build_label_manifest(df, targets, args.min_s, args.min_r))),
        "n_label_rows_with_vectors": int(len(label_manifest)),
        "n_unique_spectra_attempted": int(len(unique_spots)),
        "n_unique_spectra_vectorized": int(vector_array.shape[0]),
        "n_failures": int(len(failure_df)),
        "vector_shape": list(vector_array.shape),
        "target_mz_min": TARGET_MZ_MIN,
        "target_mz_max": TARGET_MZ_MAX,
        "n_bins": N_BINS,
        "targets": [target.__dict__ for target in targets],
    }
    (output_dir / "marisma_preprocess_report.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_markdown_report(output_dir, summary, label_manifest, failure_df)
    log(f"[preprocess] wrote outputs to {output_dir}")
    log(json.dumps(summary, indent=2))


def load_mega_module(path: Path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Mega_Model.py not found: {path}")
    spec = importlib.util.spec_from_file_location("Mega_Model_marisma_predict", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import Mega_Model.py from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["Mega_Model_marisma_predict"] = module
    spec.loader.exec_module(module)
    return module


def resolve_run_dir(run_dir: Path | None) -> Path:
    candidates = []
    if run_dir is not None:
        candidates.append(Path(run_dir))
    candidates.extend(
        [
            Path("/kaggle/input/newruns/runs/exp_ecoli_mechanism6_drugid_mae30"),
            Path("/kaggle/input/private-data-source/runs/exp_ecoli_mechanism6_drugid_mae30"),
            Path("/kaggle/input/datasets/bfdf121/newruns/runs/exp_ecoli_mechanism6_drugid_mae30"),
            Path("/kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30"),
        ]
    )
    checked: list[Path] = []
    for candidate in candidates:
        checked.append(candidate)
        if (candidate / "config.json").exists():
            return candidate
    checked_text = "\n  ".join(str(path) for path in checked)
    raise FileNotFoundError(f"Mega run config.json not found. Checked:\n  {checked_text}")


def resolve_prediction_inputs(args: argparse.Namespace) -> tuple[Path, Path]:
    vectors = args.vectors_npy or (args.output_dir / "marisma_vectors_6000.npy")
    manifest = args.manifest_csv or (args.output_dir / "marisma_prediction_manifest.csv")
    if not vectors.exists():
        raise FileNotFoundError(f"Missing MARISMa vector matrix: {vectors}")
    if not manifest.exists():
        raise FileNotFoundError(f"Missing MARISMa prediction manifest: {manifest}")
    return vectors, manifest


def resolve_checkpoints(config: dict, run_dir: Path) -> list[Path]:
    configured = config.get("ckpt_dir")
    checkpoint_dirs = []
    if configured:
        checkpoint_dirs.append(Path(configured))
        configured_path = Path(configured)
        if configured_path.is_absolute():
            checkpoint_dirs.append(run_dir / "models")
    else:
        checkpoint_dirs.append(run_dir / "models")
    checkpoint_dirs.append(run_dir / "models")

    selected = config.get("selected_seed_indices") or []
    for ckpt_dir in checkpoint_dirs:
        if not ckpt_dir.exists():
            continue
        if selected:
            checkpoints = [ckpt_dir / f"maldi_amr_seed{int(seed)}.pt" for seed in selected]
        else:
            checkpoints = sorted(ckpt_dir.glob("maldi_amr_seed*.pt"))
        checkpoints = [path for path in checkpoints if path.exists()]
        if checkpoints:
            return checkpoints

    checked = "\n  ".join(str(path) for path in checkpoint_dirs)
    raise FileNotFoundError(f"No Mega checkpoints found. Checked:\n  {checked}")


def load_mega_models(mega, config: dict, run_dir: Path):
    import torch

    drug_conditioning = config.get("drug_conditioning") or "task_id"
    checkpoints = resolve_checkpoints(config, run_dir)
    models = []
    for checkpoint in checkpoints:
        state = torch.load(checkpoint, map_location=mega.DEVICE)
        if any(str(key).startswith("module.") for key in state):
            state = {str(key).removeprefix("module."): value for key, value in state.items()}
        model = mega.create_maldi_model(
            n_sites=mega.N_SITES,
            n_organisms=mega.N_ORGANISMS,
            drug_conditioning=drug_conditioning,
        ).to(mega.DEVICE)
        model.load_state_dict(state)
        model.eval()
        models.append(model)
        log(f"[predict] loaded checkpoint: {checkpoint}")
    return models


def augment_for_tta(x: np.ndarray, mega) -> np.ndarray:
    """Apply Mega's test-time augmentation when requested."""
    try:
        return mega.augment(np.asarray(x, dtype=np.float32))
    except Exception:
        return np.asarray(x, dtype=np.float32)


def predict_vectors_for_rows(
    mega,
    models,
    vectors: np.ndarray,
    rows: pd.DataFrame,
    temperature: float,
    batch_size: int,
    tta_passes: int,
) -> np.ndarray:
    import torch

    probs_by_model = []
    vector_indices = rows["vector_index"].astype(int).to_numpy()
    org_ids = rows["org_id"].astype(int).to_numpy()

    for model_index, model in enumerate(models, start=1):
        pass_probs = []
        model.eval()
        for pass_index in range(max(1, tta_passes)):
            out_chunks = []
            for start in range(0, len(rows), batch_size):
                end = min(start + batch_size, len(rows))
                x_np = vectors[vector_indices[start:end]]
                if tta_passes > 1:
                    x_np = np.stack([augment_for_tta(row, mega) for row in x_np]).astype(np.float32)
                x = torch.from_numpy(np.asarray(x_np, dtype=np.float32)).unsqueeze(1).to(mega.DEVICE)
                org = torch.from_numpy(org_ids[start:end].astype(np.int64)).to(mega.DEVICE)
                with torch.no_grad():
                    logits = model(x, org) / temperature
                    out_chunks.append(torch.sigmoid(logits).detach().cpu().numpy())
            pass_probs.append(np.concatenate(out_chunks))
            log(
                f"[predict] model {model_index}/{len(models)} "
                f"pass {pass_index + 1}/{max(1, tta_passes)} complete"
            )
        probs_by_model.append(np.mean(pass_probs, axis=0))
    return np.mean(probs_by_model, axis=0)


def predict(args: argparse.Namespace) -> None:
    run_dir = resolve_run_dir(args.run_dir)
    vectors_path, manifest_path = resolve_prediction_inputs(args)
    prediction_csv = args.prediction_csv or (args.output_dir / "marisma_mega_predictions_long.csv")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log(f"[predict] using Mega run: {run_dir}")
    log(f"[predict] vectors: {vectors_path}")
    log(f"[predict] manifest: {manifest_path}")

    config = json.loads((run_dir / "config.json").read_text())
    mega = load_mega_module(args.mega_model_path)
    pair_profile = config.get("pair_profile", "ecoli_mechanism6")
    mega.init_config(pair_profile)

    active_pairs = [tuple(row) for row in config.get("active_pairs", [])]
    if not active_pairs:
        active_pairs = [(i, org, drug) for i, (org, drug) in enumerate(mega.ORGANISM_DRUG_PAIRS)]
    pair_to_org_id = {(organism, drug): int(org_id) for org_id, organism, drug in active_pairs}
    log(f"[predict] active pairs in checkpoint: {active_pairs}")

    manifest = pd.read_csv(manifest_path)
    keep_rows = []
    for (organism, drug), org_id in pair_to_org_id.items():
        pair_rows = manifest[manifest["organism"].eq(organism) & manifest["drug"].eq(drug)].copy()
        if pair_rows.empty:
            log(f"[predict] no MARISMa rows for checkpoint pair: {organism} / {drug}")
            continue
        pair_rows["org_id"] = org_id
        keep_rows.append(pair_rows)
        n_r = int(pair_rows["label"].sum())
        log(f"[predict] matched {organism} / {drug}: n={len(pair_rows):,} R={n_r:,}")

    if not keep_rows:
        raise ValueError("No MARISMa manifest rows overlap the Mega checkpoint active pairs")
    pred_rows = pd.concat(keep_rows, ignore_index=True)
    if args.max_prediction_rows is not None:
        pred_rows = pred_rows.head(args.max_prediction_rows).copy()
        log(f"[predict] max_prediction_rows cap active: {len(pred_rows):,}")

    vectors = np.load(vectors_path, mmap_mode="r")
    if vectors.ndim != 2 or vectors.shape[1] != N_BINS:
        raise ValueError(f"Expected vectors with shape [n, {N_BINS}], got {vectors.shape}")
    if int(pred_rows["vector_index"].max()) >= vectors.shape[0]:
        raise ValueError("Manifest vector_index exceeds vector matrix rows")

    models = load_mega_models(mega, config, run_dir)
    temperature = float(config.get("temperature", 1.0))
    tta_passes = int(args.tta_passes if args.tta_passes is not None else config.get("tta_passes", 1))
    log(
        f"[predict] scoring {len(pred_rows):,} rows using {len(models)} checkpoints, "
        f"temperature={temperature:.4f}, tta_passes={tta_passes}"
    )

    probs = predict_vectors_for_rows(
        mega=mega,
        models=models,
        vectors=vectors,
        rows=pred_rows,
        temperature=temperature,
        batch_size=args.batch_size,
        tta_passes=tta_passes,
    )

    out = pred_rows.copy()
    out.insert(0, "model_name", args.model_name)
    out["prob"] = probs.astype(float)
    out["score"] = out["prob"]
    output_columns = [
        "model_name",
        "site",
        "year",
        "isolate_id",
        "spot_id",
        "organism",
        "drug",
        "paper_drug",
        "marisma_drug",
        "drug_relationship",
        "ecology_block",
        "label",
        "prob",
        "score",
        "vector_index",
        "path",
    ]
    background_columns = [col for col in out.columns if col.startswith("background__")]
    output_columns.extend(background_columns)
    output_columns = [col for col in output_columns if col in out.columns]
    out[output_columns].to_csv(prediction_csv, index=False)

    report = {
        "run_dir": str(run_dir),
        "pair_profile": pair_profile,
        "n_prediction_rows": int(len(out)),
        "n_unique_vectors_scored": int(out["vector_index"].nunique()),
        "n_models": int(len(models)),
        "temperature": temperature,
        "tta_passes": tta_passes,
        "prediction_csv": str(prediction_csv),
        "pairs": (
            out.groupby(["organism", "drug", "marisma_drug"])
            .agg(n=("label", "size"), resistant=("label", "sum"))
            .reset_index()
            .to_dict("records")
        ),
    }
    (args.output_dir / "marisma_prediction_report.json").write_text(json.dumps(report, indent=2) + "\n")
    log(f"[predict] wrote {prediction_csv} ({len(out):,} rows)")
    log(json.dumps(report, indent=2))


def aggregate_predictions_to_isolate_drug(prediction_csv: Path, output_csv: Path, report_json: Path) -> pd.DataFrame:
    df = pd.read_csv(prediction_csv)
    required = {"site", "year", "isolate_id", "organism", "drug", "label", "prob"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{prediction_csv} is missing required columns: {', '.join(missing)}")

    key_cols = ["site", "year", "isolate_id", "organism", "drug"]
    exact_cols = [col for col in [*key_cols, "spot_id", "label", "prob"] if col in df.columns]
    duplicate_key_rows = int(df.duplicated(key_cols).sum())
    exact_duplicate_rows = int(df.duplicated(exact_cols).sum()) if exact_cols else 0

    rows: list[dict[str, object]] = []
    conflicts: list[dict[str, object]] = []
    for key, group in df.groupby(key_cols, dropna=False, sort=False):
        labels = sorted({int(float(value)) for value in group["label"].dropna().unique()})
        if len(labels) != 1:
            conflicts.append({col: value for col, value in zip(key_cols, key)} | {"labels": labels})
            continue
        probs = pd.to_numeric(group["prob"], errors="coerce").dropna()
        if probs.empty:
            continue
        row = {col: value for col, value in zip(key_cols, key)}
        row["label"] = labels[0]
        row["prob"] = float(probs.mean())
        row["n_prediction_rows"] = int(len(group))
        row["n_unique_spots"] = int(group["spot_id"].nunique()) if "spot_id" in group.columns else int(len(group))
        row["prob_sd"] = float(probs.std(ddof=0)) if len(probs) > 1 else 0.0
        if "model_name" in group.columns:
            row["model_name"] = str(group["model_name"].dropna().iloc[0]) if group["model_name"].notna().any() else ""
        rows.append(row)

    aggregated = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(output_csv, index=False)
    conflict_csv = output_csv.with_name("marisma_conflicting_isolate_drug_rows.csv")
    if conflicts:
        pd.DataFrame(conflicts).to_csv(conflict_csv, index=False)
    elif conflict_csv.exists():
        conflict_csv.unlink()
    report = {
        "input_prediction_csv": str(prediction_csv),
        "output_prediction_csv": str(output_csv),
        "conflict_csv": str(conflict_csv) if conflicts else "",
        "input_rows": int(len(df)),
        "isolate_drug_rows": int(len(aggregated)),
        "conflicting_isolate_drug_groups_excluded": int(len(conflicts)),
        "duplicate_site_year_isolate_drug_rows": duplicate_key_rows,
        "exact_duplicate_rows": exact_duplicate_rows,
        "max_prediction_rows_per_isolate_drug": int(aggregated["n_prediction_rows"].max()) if not aggregated.empty else 0,
        "max_unique_spots_per_isolate_drug": int(aggregated["n_unique_spots"].max()) if not aggregated.empty else 0,
    }
    report_json.write_text(json.dumps(report, indent=2) + "\n")
    return aggregated


def audit(args: argparse.Namespace) -> None:
    prediction_csv = args.prediction_csv or (args.output_dir / "marisma_mega_predictions_long.csv")
    if not prediction_csv.exists():
        raise FileNotFoundError(f"Missing MARISMa prediction CSV: {prediction_csv}")

    audit_output_dir = args.audit_output_dir or (args.output_dir / "marisma_isolate_background_audit")
    audit_output_dir.mkdir(parents=True, exist_ok=True)
    aggregate_csv = audit_output_dir / "marisma_isolate_level_predictions.csv"
    duplicate_report = audit_output_dir / "marisma_duplicate_handling_report.json"

    aggregated = aggregate_predictions_to_isolate_drug(prediction_csv, aggregate_csv, duplicate_report)
    log(
        "[audit] aggregated MARISMa predictions: "
        f"{prediction_csv} -> {aggregate_csv} ({len(aggregated):,} isolate/drug rows)"
    )

    audit_script = args.audit_script
    if not audit_script.is_absolute():
        audit_script = Path(__file__).resolve().parents[1] / audit_script
    cmd = [
        sys.executable,
        str(audit_script),
        "--predictions-csv",
        str(aggregate_csv),
        "--output-dir",
        str(audit_output_dir),
        "--id-col",
        "isolate_id",
        "--model-name",
        args.model_name,
        "--bootstrap-n",
        str(args.bootstrap_n),
        "--permutation-n",
        str(args.permutation_n),
    ]
    log("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    args = build_parser().parse_args()
    if args.stage == "preprocess":
        preprocess(args)
    elif args.stage == "all":
        preprocess(args)
        predict(args)
    elif args.stage == "predict":
        predict(args)
    elif args.stage == "audit":
        audit(args)


if __name__ == "__main__":
    main()
