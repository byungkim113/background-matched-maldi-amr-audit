#!/usr/bin/env python3
"""Falsification controls for background-aware MALDI-AMR prediction tables."""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.calibration_analysis import mean, normalize_label, parse_float, safe_auc


DEFAULT_PREDICTIONS = str(
    ROOT / "outputs" / "analysis_outputs" / "background_matched_predictions.csv"
)
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/falsification_controls"

FIELDS = [
    "site",
    "organism",
    "drug",
    "n",
    "n_r",
    "n_s",
    "observed_auc",
    "background_burden_auc",
    "observed_minus_burden_auc",
    "score_background_burden_correlation",
    "shuffle_null_mean_auc",
    "shuffle_null_sd_auc",
    "shuffle_empirical_p_ge_observed",
    "observed_minus_shuffle_null_auc",
    "permutations",
    "adequacy_label",
    "control_interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run background-only and within-background shuffle controls.")
    parser.add_argument("--predictions-csv", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--min-pos", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=3)
    return parser.parse_args()


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    pairs = [(float(x), float(y)) for x, y in zip(xs, ys) if math.isfinite(float(x)) and math.isfinite(float(y))]
    if len(pairs) < 2:
        return math.nan
    x_vals = [x for x, _ in pairs]
    y_vals = [y for _, y in pairs]
    x_mean = mean(x_vals)
    y_mean = mean(y_vals)
    num = sum((x - x_mean) * (y - y_mean) for x, y in pairs)
    den_x = math.sqrt(sum((x - x_mean) ** 2 for x in x_vals))
    den_y = math.sqrt(sum((y - y_mean) ** 2 for y in y_vals))
    return num / (den_x * den_y) if den_x and den_y else math.nan


def adequacy_label(n: int, n_pos: int, n_neg: int, min_n: int, min_pos: int, min_neg: int) -> str:
    if n < min_n:
        return "underpowered_low_n"
    if n_pos < min_pos:
        return "underpowered_low_resistant"
    if n_neg < min_neg:
        return "underpowered_low_susceptible"
    return "interpretable"


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
        row["background_signature"] = str(record.get("background_signature", "") or "NO_BACKGROUND_SIGNATURE")
        rows.append(row)
    return rows


def shuffled_auc_within_background(group: Sequence[dict], rng: random.Random) -> float:
    by_background: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(group):
        by_background[row["background_signature"]].append(idx)

    labels = [int(row["label"]) for row in group]
    shuffled = list(labels)
    for indices in by_background.values():
        values = [shuffled[idx] for idx in indices]
        rng.shuffle(values)
        for idx, value in zip(indices, values):
            shuffled[idx] = value
    return safe_auc(shuffled, [float(row["prob"]) for row in group])


def shuffle_null(group: Sequence[dict], n_permutations: int, seed: int) -> tuple[float, float, float]:
    observed = safe_auc([int(row["label"]) for row in group], [float(row["prob"]) for row in group])
    values = []
    rng = random.Random(seed)
    for _ in range(max(0, n_permutations)):
        auc = shuffled_auc_within_background(group, rng)
        if math.isfinite(auc):
            values.append(auc)
    if not values:
        return math.nan, math.nan, math.nan
    null_mean = mean(values)
    null_sd = math.sqrt(mean((value - null_mean) ** 2 for value in values))
    p_ge = (sum(value >= observed - 1e-12 for value in values) + 1) / (len(values) + 1)
    return null_mean, null_sd, p_ge


def interpretation(observed_auc: float, burden_auc: float, null_mean: float, adequacy: str) -> str:
    if adequacy != "interpretable" or not math.isfinite(observed_auc):
        return "do_not_overinterpret"
    observed_minus_burden = observed_auc - burden_auc if math.isfinite(burden_auc) else math.nan
    observed_minus_null = observed_auc - null_mean if math.isfinite(null_mean) else math.nan
    if (
        math.isfinite(observed_minus_burden)
        and math.isfinite(observed_minus_null)
        and observed_minus_burden >= 0.05
        and observed_minus_null >= 0.05
    ):
        return "focal_score_exceeds_controls"
    if math.isfinite(observed_minus_null) and observed_minus_null >= 0.05:
        return "score_exceeds_shuffle_but_burden_competitive"
    return "background_control_competitive"


def build_falsification_rows(
    records: Iterable[dict],
    *,
    n_permutations: int = 500,
    min_n: int = 30,
    min_pos: int = 3,
    min_neg: int = 3,
    seed: int = 17,
) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in normalized_rows(records):
        grouped[(row.get("site", ""), row.get("organism", ""), row.get("drug", ""))].append(row)

    output = []
    for group_index, ((site, organism, drug), group) in enumerate(sorted(grouped.items())):
        labels = [int(row["label"]) for row in group]
        scores = [float(row["prob"]) for row in group]
        burdens = [float(row["background_resistant_count"]) for row in group]
        n_pos = sum(labels)
        n_neg = len(group) - n_pos
        observed = safe_auc(labels, scores)
        burden_auc = safe_auc(labels, burdens)
        null_mean, null_sd, p_ge = shuffle_null(group, n_permutations, seed + group_index)
        adequacy = adequacy_label(len(group), n_pos, n_neg, min_n, min_pos, min_neg)
        output.append(
            {
                "site": site,
                "organism": organism,
                "drug": drug,
                "n": len(group),
                "n_r": n_pos,
                "n_s": n_neg,
                "observed_auc": observed,
                "background_burden_auc": burden_auc,
                "observed_minus_burden_auc": observed - burden_auc if math.isfinite(observed) and math.isfinite(burden_auc) else math.nan,
                "score_background_burden_correlation": pearson(scores, burdens),
                "shuffle_null_mean_auc": null_mean,
                "shuffle_null_sd_auc": null_sd,
                "shuffle_empirical_p_ge_observed": p_ge,
                "observed_minus_shuffle_null_auc": observed - null_mean if math.isfinite(observed) and math.isfinite(null_mean) else math.nan,
                "permutations": n_permutations,
                "adequacy_label": adequacy,
                "control_interpretation": interpretation(observed, burden_auc, null_mean, adequacy),
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
        rows = build_falsification_rows(
            csv.DictReader(f),
            n_permutations=args.permutations,
            min_n=args.min_n,
            min_pos=args.min_pos,
            min_neg=args.min_neg,
            seed=args.seed,
        )
    write_csv(output_dir / "falsification_controls.csv", rows)
    md = [
        "# Falsification Controls",
        "",
        "Compares observed model AUC with background-burden-only AUC and a within-background label-shuffle null.",
        "",
        f"- Prediction rows: `{args.predictions_csv}`",
        f"- Permutations: {args.permutations}",
        "",
        markdown_table(rows),
        "",
    ]
    (output_dir / "falsification_controls.md").write_text("\n".join(md))
    print(f"Wrote {len(rows)} falsification rows to {output_dir}")


if __name__ == "__main__":
    main()
