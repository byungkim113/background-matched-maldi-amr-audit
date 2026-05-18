#!/usr/bin/env python3
"""Stratify focal-drug performance by observable co-resistance background.

This is a companion to the background-matched audit. It does not claim causal
control; it describes where model signal is concentrated across the non-focal
AST ecology available in a prediction table.
"""

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
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/co_resistance_stratification"

FIELDS = [
    "site",
    "organism",
    "drug",
    "stratum_type",
    "stratum_value",
    "n",
    "n_r",
    "n_s",
    "auc",
    "mean_prob",
    "mean_prob_r",
    "mean_prob_s",
    "background_known_count",
    "background_resistant_count",
    "adequacy_label",
    "interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute focal-drug AUC inside co-resistance burden and exact-background strata."
    )
    parser.add_argument("--predictions-csv", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--min-pos", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=3)
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


def parse_background_signature(signature: str) -> tuple[int, int]:
    known = 0
    resistant = 0
    for part in str(signature or "").split("|"):
        if "=" not in part:
            continue
        value = part.rsplit("=", 1)[-1].strip().upper()
        if value in {"S", "R", "0", "1"}:
            known += 1
            resistant += int(value in {"R", "1"})
    return known, resistant


def get_background_counts(row: dict) -> tuple[int, int]:
    known = parse_float(row.get("background_known_count", ""))
    resistant = parse_float(row.get("background_resistant_count", ""))
    if known is not None and resistant is not None:
        return int(known), int(resistant)
    return parse_background_signature(str(row.get("background_signature", "")))


def burden_bin(known: int, resistant: int) -> str:
    if known <= 0:
        return "unknown"
    if resistant == 0:
        return "none"
    if resistant == 1:
        return "single"
    return "multi"


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


def adequacy_label(n: int, n_pos: int, n_neg: int, min_n: int, min_pos: int, min_neg: int) -> str:
    if n < min_n:
        return "underpowered_low_n"
    if n_pos < min_pos:
        return "underpowered_low_resistant"
    if n_neg < min_neg:
        return "underpowered_low_susceptible"
    return "interpretable"


def interpretation_for(auc: float, adequacy: str) -> str:
    if adequacy != "interpretable" or not math.isfinite(auc):
        return "do_not_overinterpret"
    if auc >= 0.60:
        return "signal_retained_in_stratum"
    if auc < 0.55:
        return "weak_or_reversed_in_stratum"
    return "borderline_in_stratum"


def normalized_rows(records: Iterable[dict]) -> list[dict]:
    rows = []
    for record in records:
        label = normalize_label(record.get("label"))
        prob = parse_float(record.get("prob"))
        if label is None or prob is None:
            continue
        known, resistant = get_background_counts(record)
        row = dict(record)
        row["label"] = label
        row["prob"] = prob
        row["background_known_count"] = known
        row["background_resistant_count"] = resistant
        row["background_signature"] = str(row.get("background_signature", "")).strip() or "NO_BACKGROUND_SIGNATURE"
        rows.append(row)
    return rows


def summarize_group(
    group: Sequence[dict],
    *,
    stratum_type: str,
    stratum_value: str,
    min_n: int,
    min_pos: int,
    min_neg: int,
) -> dict:
    labels = [int(row["label"]) for row in group]
    scores = [float(row["prob"]) for row in group]
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    auc = safe_auc(labels, scores)
    adequacy = adequacy_label(len(group), n_pos, n_neg, min_n, min_pos, min_neg)
    first = group[0]
    return {
        "site": first.get("site", ""),
        "organism": first.get("organism", ""),
        "drug": first.get("drug", ""),
        "stratum_type": stratum_type,
        "stratum_value": stratum_value,
        "n": len(group),
        "n_r": n_pos,
        "n_s": n_neg,
        "auc": auc,
        "mean_prob": mean(scores),
        "mean_prob_r": mean(score for label, score in zip(labels, scores) if label == 1),
        "mean_prob_s": mean(score for label, score in zip(labels, scores) if label == 0),
        "background_known_count": first.get("background_known_count", ""),
        "background_resistant_count": first.get("background_resistant_count", ""),
        "adequacy_label": adequacy,
        "interpretation": interpretation_for(auc, adequacy),
    }


def build_stratification_rows(records: Iterable[dict], min_n: int = 30, min_pos: int = 3, min_neg: int = 3) -> list[dict]:
    rows = normalized_rows(records)
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        base = (row.get("site", ""), row.get("organism", ""), row.get("drug", ""))
        grouped[base + ("background_resistance_burden", str(row["background_resistant_count"]))].append(row)
        grouped[base + ("background_resistance_burden_bin", burden_bin(row["background_known_count"], row["background_resistant_count"]))].append(row)
        grouped[base + ("exact_background_signature", row["background_signature"])].append(row)

    summaries = []
    for key, group in grouped.items():
        site, organism, drug, stratum_type, stratum_value = key
        summaries.append(
            summarize_group(
                group,
                stratum_type=stratum_type,
                stratum_value=str(stratum_value),
                min_n=min_n,
                min_pos=min_pos,
                min_neg=min_neg,
            )
        )
    return sorted(
        summaries,
        key=lambda row: (
            row["site"],
            row["organism"],
            row["drug"],
            row["stratum_type"],
            str(row["stratum_value"]),
        ),
    )


def format_value(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.6g}"
    return value


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in FIELDS})


def markdown_table(rows: Sequence[dict], limit: int = 120) -> str:
    display = list(rows[:limit])
    if not display:
        return "_No rows._"
    rendered = []
    for row in display:
        rendered.append({field: str(format_value(row.get(field, ""))) for field in FIELDS})
    widths = {field: max(len(field), *(len(row[field]) for row in rendered)) for field in FIELDS}
    header = "| " + " | ".join(field.ljust(widths[field]) for field in FIELDS) + " |"
    sep = "| " + " | ".join("-" * widths[field] for field in FIELDS) + " |"
    body = [
        "| " + " | ".join(row[field].ljust(widths[field]) for field in FIELDS) + " |"
        for row in rendered
    ]
    suffix = []
    if len(rows) > limit:
        suffix = ["", f"_Showing first {limit} of {len(rows)} rows. See CSV for complete output._"]
    return "\n".join([header, sep] + body + suffix)


def write_markdown(path: Path, rows: Sequence[dict], args: argparse.Namespace) -> None:
    interpretable = sum(1 for row in rows if row["adequacy_label"] == "interpretable")
    lines = [
        "# Co-Resistance Stratification",
        "",
        "This companion analysis reports focal-drug performance inside observable co-resistance strata.",
        "It is descriptive: co-resistance strata can mix lineage, plasmid background, selection ecology, and other biological context.",
        "",
        f"- Prediction rows: `{args.predictions_csv}`",
        f"- Minimum support: n>={args.min_n}, R>={args.min_pos}, S>={args.min_neg}",
        f"- Strata reported: {len(rows)} total; {interpretable} interpretable",
        "",
        markdown_table(rows),
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    predictions_csv = Path(args.predictions_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with predictions_csv.open(newline="") as f:
        rows = build_stratification_rows(
            csv.DictReader(f),
            min_n=args.min_n,
            min_pos=args.min_pos,
            min_neg=args.min_neg,
        )

    csv_path = output_dir / "co_resistance_stratification.csv"
    md_path = output_dir / "co_resistance_stratification.md"
    write_csv(csv_path, rows)
    write_markdown(md_path, rows, args)
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
