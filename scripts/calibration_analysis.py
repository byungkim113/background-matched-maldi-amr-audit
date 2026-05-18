#!/usr/bin/env python3
"""Compute calibration and threshold metrics for long MALDI-AMR predictions."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_PREDICTIONS = str(
    Path(__file__).resolve().parents[1]
    / "outputs" / "analysis_outputs" / "background_matched_predictions.csv"
)
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/calibration_analysis"

FIELDS = [
    "site",
    "organism",
    "drug",
    "n",
    "n_r",
    "n_s",
    "prevalence",
    "mean_prob",
    "auc",
    "brier",
    "expected_calibration_error",
    "maximum_calibration_error",
    "threshold",
    "sensitivity",
    "specificity",
    "ppv",
    "npv",
    "accuracy",
    "balanced_accuracy",
    "calibration_label",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute calibration metrics for long prediction CSVs.")
    parser.add_argument("--predictions-csv", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-bins", type=int, default=10)
    parser.add_argument("--threshold", type=float, default=0.5)
    return parser.parse_args()


def normalize_label(value: object) -> int | None:
    text = str(value).strip().upper()
    if text in {"1", "R", "RESISTANT", "TRUE", "T"}:
        return 1
    if text in {"0", "S", "SUSCEPTIBLE", "FALSE", "F"}:
        return 0
    return None


def parse_float(value: object) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def safe_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    pairs = [(int(y), float(s)) for y, s in zip(labels, scores) if math.isfinite(float(s))]
    if not pairs:
        return math.nan
    n_pos = sum(y for y, _ in pairs)
    n_neg = len(pairs) - n_pos
    if n_pos == 0 or n_neg == 0:
        return math.nan

    order = sorted(range(len(pairs)), key=lambda idx: pairs[idx][1])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pairs[order[j + 1]][1] == pairs[order[i]][1]:
            j += 1
        rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = rank
        i = j + 1

    pos_rank_sum = sum(rank for rank, (label, _) in zip(ranks, pairs) if label == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def mean(values: Iterable[float]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else math.nan


def brier_score(labels: Sequence[int], scores: Sequence[float]) -> float:
    return mean((float(score) - int(label)) ** 2 for label, score in zip(labels, scores))


def calibration_errors(labels: Sequence[int], scores: Sequence[float], n_bins: int) -> tuple[float, float]:
    bins: list[list[tuple[int, float]]] = [[] for _ in range(max(1, n_bins))]
    for label, score in zip(labels, scores):
        idx = min(len(bins) - 1, max(0, int(float(score) * len(bins))))
        bins[idx].append((int(label), float(score)))

    total = len(labels)
    weighted_error = 0.0
    max_error = 0.0
    for bucket in bins:
        if not bucket:
            continue
        bucket_labels = [label for label, _ in bucket]
        bucket_scores = [score for _, score in bucket]
        error = abs(mean(bucket_scores) - mean(bucket_labels))
        weighted_error += (len(bucket) / total) * error
        max_error = max(max_error, error)
    return weighted_error, max_error


def safe_div(num: float, den: float) -> float:
    return num / den if den else math.nan


def threshold_metrics(labels: Sequence[int], scores: Sequence[float], threshold: float) -> dict:
    tp = fp = tn = fn = 0
    for label, score in zip(labels, scores):
        pred = 1 if float(score) >= threshold else 0
        if label == 1 and pred == 1:
            tp += 1
        elif label == 0 and pred == 1:
            fp += 1
        elif label == 0 and pred == 0:
            tn += 1
        elif label == 1 and pred == 0:
            fn += 1

    sensitivity = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    return {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "ppv": safe_div(tp, tp + fp),
        "npv": safe_div(tn, tn + fn),
        "accuracy": safe_div(tp + tn, tp + tn + fp + fn),
        "balanced_accuracy": mean(v for v in [sensitivity, specificity] if math.isfinite(v)),
    }


def calibration_label(ece: float, brier: float) -> str:
    if math.isfinite(ece) and math.isfinite(brier) and ece <= 0.15 and brier <= 0.20:
        return "well_calibrated"
    if math.isfinite(ece) and math.isfinite(brier) and ece <= 0.25 and brier <= 0.30:
        return "borderline_calibration"
    return "poorly_calibrated"


def normalized_rows(records: Iterable[dict]) -> list[dict]:
    rows = []
    for record in records:
        label = normalize_label(record.get("label"))
        prob = parse_float(record.get("prob"))
        if label is None or prob is None:
            continue
        row = dict(record)
        row["label"] = label
        row["prob"] = min(1.0, max(0.0, prob))
        rows.append(row)
    return rows


def build_calibration_rows(records: Iterable[dict], n_bins: int = 10, threshold: float = 0.5) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in normalized_rows(records):
        grouped[(row.get("site", ""), row.get("organism", ""), row.get("drug", ""))].append(row)

    output = []
    for (site, organism, drug), group in sorted(grouped.items()):
        labels = [int(row["label"]) for row in group]
        scores = [float(row["prob"]) for row in group]
        ece, mce = calibration_errors(labels, scores, n_bins)
        metrics = threshold_metrics(labels, scores, threshold)
        n_r = sum(labels)
        brier = brier_score(labels, scores)
        output.append(
            {
                "site": site,
                "organism": organism,
                "drug": drug,
                "n": len(group),
                "n_r": n_r,
                "n_s": len(group) - n_r,
                "prevalence": mean(labels),
                "mean_prob": mean(scores),
                "auc": safe_auc(labels, scores),
                "brier": brier,
                "expected_calibration_error": ece,
                "maximum_calibration_error": mce,
                "threshold": threshold,
                **metrics,
                "calibration_label": calibration_label(ece, brier),
            }
        )
    return output


def format_value(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.6g}"
    return value


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in FIELDS})


def markdown_table(rows: Sequence[dict], limit: int = 80) -> str:
    display = list(rows[:limit])
    if not display:
        return "_No rows._"
    headers = FIELDS
    rendered = [{field: str(format_value(row.get(field, ""))) for field in headers} for row in display]
    widths = {field: max(len(field), *(len(row[field]) for row in rendered)) for field in headers}
    lines = [
        "| " + " | ".join(field.ljust(widths[field]) for field in headers) + " |",
        "| " + " | ".join("-" * widths[field] for field in headers) + " |",
    ]
    for row in rendered:
        lines.append("| " + " | ".join(row[field].ljust(widths[field]) for field in headers) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    with Path(args.predictions_csv).open() as f:
        rows = build_calibration_rows(csv.DictReader(f), n_bins=args.n_bins, threshold=args.threshold)
    write_csv(output_dir / "calibration_summary.csv", rows)
    md = [
        "# Calibration Analysis",
        "",
        "Grouped Brier score, expected calibration error, and threshold metrics for long MALDI-AMR predictions.",
        "",
        f"- Prediction rows: `{args.predictions_csv}`",
        f"- Bins: {args.n_bins}",
        f"- Threshold: {args.threshold}",
        "",
        markdown_table(rows),
        "",
    ]
    (output_dir / "calibration_summary.md").write_text("\n".join(md))
    print(f"Wrote {len(rows)} calibration rows to {output_dir}")


if __name__ == "__main__":
    main()
