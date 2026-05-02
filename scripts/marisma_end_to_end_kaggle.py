#!/usr/bin/env python3
"""MARISMa v2 preprocessing pipeline for background-matched MALDI-AMR audits.

This script is intentionally staged. The current reliable stage is
``preprocess``: it reads MARISMa Bruker spectra, converts them into the
6000-bin representation used by the DRIAMS/Mega experiments, and writes a
long manifest of organism-drug labels. The later ``predict`` and ``audit``
stages are placeholders until the locked Mega checkpoints/config are mounted.

Example Kaggle use:

    python /kaggle/working/marisma_end_to_end_kaggle.py \\
      --stage preprocess \\
      --amr-csv /kaggle/input/datasets/bfdf121/marisma/AMR.csv \\
      --marisma-root /kaggle/input/datasets/bfdf121/marisma/MARISMa \\
      --output-dir /kaggle/working/marisma_preprocessed
"""

from __future__ import annotations

import argparse
import json
import math
import re
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
    p = argparse.ArgumentParser(description="Preprocess MARISMa v2 spectra for MALDI-AMR background audits.")
    p.add_argument(
        "--stage",
        choices=["preprocess", "predict", "audit", "all"],
        default="preprocess",
        help="Currently only preprocess is implemented; predict/audit are guarded placeholders.",
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


def guarded_not_implemented(stage: str) -> None:
    raise NotImplementedError(
        f"Stage '{stage}' is intentionally not implemented yet. First generate "
        "marisma_vectors_6000.npy and marisma_prediction_manifest.csv, then mount "
        "the locked Mega checkpoints/config so prediction can write "
        "marisma_mega_predictions_long.csv."
    )


def main() -> None:
    args = build_parser().parse_args()
    if args.stage == "preprocess":
        preprocess(args)
    elif args.stage == "all":
        preprocess(args)
        guarded_not_implemented("predict")
    else:
        guarded_not_implemented(args.stage)


if __name__ == "__main__":
    main()
