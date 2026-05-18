#!/usr/bin/env python3
"""Compare a Weis-code rerun against upstream stored Weis result JSON files."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


METRICS = ["auroc", "auprc", "accuracy"]


def safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number):
        return None
    return number


def key_from_result(row: dict) -> tuple[str, str, str, str, str, str]:
    return (
        str(row.get("model", "")),
        str(row.get("train_site", "")),
        str(row.get("test_site", "")),
        str(row.get("species", "")),
        str(row.get("drug", row.get("antibiotic", ""))),
        str(row.get("seed", "")),
    )


def load_reference_results(root: Path) -> dict[tuple[str, str, str, str, str, str], dict]:
    index: dict[tuple[str, str, str, str, str, str], dict] = {}
    for path in sorted(root.rglob("*.json")):
        try:
            data = json.loads(path.read_text())
        except Exception:
            continue
        if not {"model", "train_site", "test_site", "species", "antibiotic", "seed"}.issubset(data):
            continue
        key = key_from_result(data)
        data["_reference_path"] = str(path)
        index[key] = data
    return index


def compare(raw_results: list[dict], reference_index: dict, tolerance: float) -> list[dict]:
    rows = []
    for raw in raw_results:
        key = key_from_result(raw)
        ref = reference_index.get(key)
        base = {
            "model": key[0],
            "train_site": key[1],
            "test_site": key[2],
            "species": key[3],
            "drug": key[4],
            "seed": key[5],
            "reference_found": bool(ref),
            "reference_path": "" if ref is None else ref.get("_reference_path", ""),
        }
        max_abs_diff = 0.0
        all_within = bool(ref)
        for metric in METRICS:
            raw_value = safe_float(raw.get(metric))
            ref_value = safe_float(ref.get(metric)) if ref is not None else None
            diff = None if raw_value is None or ref_value is None else raw_value - ref_value
            abs_diff = None if diff is None else abs(diff)
            base[f"raw_{metric}"] = raw_value
            base[f"reference_{metric}"] = ref_value
            base[f"diff_{metric}"] = diff
            if abs_diff is None:
                all_within = False
            else:
                max_abs_diff = max(max_abs_diff, abs_diff)
                if abs_diff > tolerance:
                    all_within = False
        base["max_abs_diff"] = max_abs_diff if ref is not None else None
        base["within_tolerance"] = all_within
        rows.append(base)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "model",
        "train_site",
        "test_site",
        "species",
        "drug",
        "seed",
        "reference_found",
        "reference_path",
        "raw_auroc",
        "reference_auroc",
        "diff_auroc",
        "raw_auprc",
        "reference_auprc",
        "diff_auprc",
        "raw_accuracy",
        "reference_accuracy",
        "diff_accuracy",
        "max_abs_diff",
        "within_tolerance",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict], tolerance: float) -> None:
    n = len(rows)
    found = sum(1 for row in rows if row["reference_found"])
    passing = sum(1 for row in rows if row["within_tolerance"])
    max_diff = max(
        (float(row["max_abs_diff"]) for row in rows if row["max_abs_diff"] is not None),
        default=None,
    )
    lines = [
        "# Weis Raw-Metric Parity Check",
        "",
        f"- Rows compared: {n}",
        f"- Reference rows found: {found}",
        f"- Rows within tolerance ({tolerance:g}): {passing}",
        f"- Maximum absolute metric difference: {'' if max_diff is None else f'{max_diff:.6g}'}",
        "",
    ]
    if passing == n and n > 0:
        lines.append("All rerun rows matched the upstream stored metrics within tolerance.")
    else:
        lines.append("Do not describe this run as an exact Weis et al. replication until all rows match.")
    path.write_text("\n".join(lines) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weis-raw-results", type=Path, required=True)
    parser.add_argument("--reference-results-root", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--summary-md", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    raw_results = json.loads(args.weis_raw_results.read_text())
    if not isinstance(raw_results, list):
        raise TypeError("--weis-raw-results must contain a JSON list")
    rows = compare(raw_results, load_reference_results(args.reference_results_root), args.tolerance)
    write_csv(args.output_csv, rows)
    if args.summary_md:
        write_summary(args.summary_md, rows, args.tolerance)


if __name__ == "__main__":
    main()
