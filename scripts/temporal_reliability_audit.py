#!/usr/bin/env python3
"""Monitor MALDI-AMR reliability across chronological prediction periods."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.calibration_analysis import (
    brier_score,
    calibration_errors,
    mean,
    normalize_label,
    parse_float,
    safe_auc,
)


DEFAULT_PREDICTIONS = str(
    ROOT / "outputs" / "analysis_outputs" / "background_matched_predictions.csv"
)
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/temporal_reliability_audit"

FIELDS = [
    "site",
    "organism",
    "drug",
    "period",
    "n_periods_observed",
    "n",
    "n_r",
    "n_s",
    "prevalence",
    "mean_prob",
    "auc",
    "brier",
    "expected_calibration_error",
    "mean_background_resistant_count",
    "baseline_period",
    "auc_delta_from_baseline",
    "prevalence_delta_from_baseline",
    "ece_delta_from_baseline",
    "mean_background_resistant_delta",
    "support_label",
    "reliability_status",
    "recommended_action",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assess temporal reliability and drift in prediction rows.")
    parser.add_argument("--predictions-csv", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--min-pos", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=3)
    parser.add_argument("--auc-drop-alert", type=float, default=0.10)
    parser.add_argument("--ece-alert", type=float, default=0.10)
    parser.add_argument("--prevalence-shift-alert", type=float, default=0.15)
    parser.add_argument("--burden-shift-alert", type=float, default=0.75)
    return parser.parse_args()


def support_label(n: int, n_pos: int, n_neg: int, min_n: int, min_pos: int, min_neg: int) -> str:
    if n < min_n:
        return "underpowered_low_n"
    if n_pos < min_pos:
        return "underpowered_low_resistant"
    if n_neg < min_neg:
        return "underpowered_low_susceptible"
    return "adequate"


def normalized_rows(records: Iterable[dict]) -> list[dict]:
    rows = []
    for record in records:
        label = normalize_label(record.get("label"))
        prob = parse_float(record.get("prob"))
        if label is None or prob is None:
            continue
        burden = parse_float(record.get("background_resistant_count", ""))
        row = dict(record)
        row["label"] = label
        row["prob"] = min(1.0, max(0.0, prob))
        row["background_resistant_count"] = burden if burden is not None else math.nan
        row["year"] = str(record.get("year", "") or "UNKNOWN")
        rows.append(row)
    return rows


def period_sort_key(value: str) -> tuple[int, str]:
    text = str(value)
    digits = "".join(ch for ch in text if ch.isdigit())
    return (int(digits) if digits else 10**9, text)


def summarize_period(group: Sequence[dict], min_n: int, min_pos: int, min_neg: int) -> dict:
    labels = [int(row["label"]) for row in group]
    scores = [float(row["prob"]) for row in group]
    burdens = [float(row["background_resistant_count"]) for row in group if math.isfinite(float(row["background_resistant_count"]))]
    ece, _ = calibration_errors(labels, scores, n_bins=10)
    n_pos = sum(labels)
    return {
        "n": len(group),
        "n_r": n_pos,
        "n_s": len(group) - n_pos,
        "prevalence": mean(labels),
        "mean_prob": mean(scores),
        "auc": safe_auc(labels, scores),
        "brier": brier_score(labels, scores),
        "expected_calibration_error": ece,
        "mean_background_resistant_count": mean(burdens),
        "support_label": support_label(len(group), n_pos, len(group) - n_pos, min_n, min_pos, min_neg),
    }


def finite_delta(current: float, baseline: float) -> float:
    if math.isfinite(current) and math.isfinite(baseline):
        return current - baseline
    return math.nan


def classify_period(
    *,
    period_count: int,
    support: str,
    auc_delta: float,
    prevalence_delta: float,
    ece_delta: float,
    burden_delta: float,
    auc_drop_alert: float,
    prevalence_shift_alert: float,
    ece_alert: float,
    burden_shift_alert: float,
) -> tuple[str, str]:
    if support != "adequate":
        return "insufficient_support", "collect_more_labels_before_interpreting"
    if period_count < 2:
        return "insufficient_periods", "collect_future_periods"
    alerts = [
        math.isfinite(auc_delta) and auc_delta <= -abs(auc_drop_alert),
        math.isfinite(prevalence_delta) and abs(prevalence_delta) >= abs(prevalence_shift_alert),
        math.isfinite(ece_delta) and ece_delta >= abs(ece_alert),
        math.isfinite(burden_delta) and abs(burden_delta) >= abs(burden_shift_alert),
    ]
    if any(alerts):
        return "drift_alert", "recalibration_or_retraining_review"
    return "stable_within_monitored_periods", "continue_locked_monitoring"


def build_temporal_rows(
    records: Iterable[dict],
    *,
    min_n: int = 30,
    min_pos: int = 3,
    min_neg: int = 3,
    auc_drop_alert: float = 0.10,
    prevalence_shift_alert: float = 0.15,
    ece_alert: float = 0.10,
    burden_shift_alert: float = 0.75,
) -> list[dict]:
    period_groups: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in normalized_rows(records):
        period_groups[(row.get("site", ""), row.get("organism", ""), row.get("drug", ""), row["year"])].append(row)

    by_pair: dict[tuple[str, str, str], list[tuple[str, dict]]] = defaultdict(list)
    for (site, organism, drug, period), group in period_groups.items():
        by_pair[(site, organism, drug)].append((period, summarize_period(group, min_n, min_pos, min_neg)))

    output = []
    for (site, organism, drug), periods in sorted(by_pair.items()):
        periods = sorted(periods, key=lambda item: period_sort_key(item[0]))
        baseline_period, baseline = periods[0]
        for period, summary in periods:
            auc_delta = finite_delta(summary["auc"], baseline["auc"])
            prevalence_delta = finite_delta(summary["prevalence"], baseline["prevalence"])
            ece_delta = finite_delta(summary["expected_calibration_error"], baseline["expected_calibration_error"])
            burden_delta = finite_delta(summary["mean_background_resistant_count"], baseline["mean_background_resistant_count"])
            status, action = classify_period(
                period_count=len(periods),
                support=summary["support_label"],
                auc_delta=auc_delta,
                prevalence_delta=prevalence_delta,
                ece_delta=ece_delta,
                burden_delta=burden_delta,
                auc_drop_alert=auc_drop_alert,
                prevalence_shift_alert=prevalence_shift_alert,
                ece_alert=ece_alert,
                burden_shift_alert=burden_shift_alert,
            )
            output.append(
                {
                    "site": site,
                    "organism": organism,
                    "drug": drug,
                    "period": period,
                    "n_periods_observed": len(periods),
                    **summary,
                    "baseline_period": baseline_period,
                    "auc_delta_from_baseline": auc_delta,
                    "prevalence_delta_from_baseline": prevalence_delta,
                    "ece_delta_from_baseline": ece_delta,
                    "mean_background_resistant_delta": burden_delta,
                    "reliability_status": status,
                    "recommended_action": action,
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
    rendered = [{field: str(format_value(row.get(field, ""))) for field in FIELDS} for row in display]
    widths = {field: max(len(field), *(len(row[field]) for row in rendered)) for field in FIELDS}
    lines = [
        "| " + " | ".join(field.ljust(widths[field]) for field in FIELDS) + " |",
        "| " + " | ".join("-" * widths[field] for field in FIELDS) + " |",
    ]
    for row in rendered:
        lines.append("| " + " | ".join(row[field].ljust(widths[field]) for field in FIELDS) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    with Path(args.predictions_csv).open() as f:
        rows = build_temporal_rows(
            csv.DictReader(f),
            min_n=args.min_n,
            min_pos=args.min_pos,
            min_neg=args.min_neg,
            auc_drop_alert=args.auc_drop_alert,
            prevalence_shift_alert=args.prevalence_shift_alert,
            ece_alert=args.ece_alert,
            burden_shift_alert=args.burden_shift_alert,
        )
    write_csv(output_dir / "temporal_reliability.csv", rows)
    md = [
        "# Temporal Reliability Audit",
        "",
        "Period-wise model reliability and drift monitor. If only one period is present, duration of reliability cannot be estimated.",
        "",
        f"- Prediction rows: `{args.predictions_csv}`",
        f"- Minimum support: n>={args.min_n}, R>={args.min_pos}, S>={args.min_neg}",
        "",
        markdown_table(rows),
        "",
    ]
    (output_dir / "temporal_reliability.md").write_text("\n".join(md))
    print(f"Wrote {len(rows)} temporal reliability rows to {output_dir}")


if __name__ == "__main__":
    main()
