#!/usr/bin/env python3
"""Exact co-resistance-only baselines for MALDI-AMR prediction tables."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.calibration_analysis import mean, normalize_label, parse_float, safe_auc
from scripts.falsification_controls import pearson


DEFAULT_PREDICTIONS = (
    "/Users/byungkim/Downloads/results-2/background_matched_contrastive/"
    "background_matched_predictions.csv"
)
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/co_resistance_only_baseline"

FIELDS = [
    "site",
    "organism",
    "drug",
    "n",
    "n_r",
    "n_s",
    "prevalence",
    "observed_auc",
    "background_burden_auc",
    "exact_background_auc",
    "observed_minus_burden_auc",
    "observed_minus_exact_background_auc",
    "exact_minus_burden_auc",
    "score_exact_background_correlation",
    "score_background_burden_correlation",
    "unique_backgrounds",
    "singleton_backgrounds",
    "smoothing_alpha",
    "adequacy_label",
    "baseline_interpretation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict each focal drug from non-focal AST background only, with no spectra."
    )
    parser.add_argument("--predictions-csv", default=DEFAULT_PREDICTIONS)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoothing-alpha", type=float, default=2.0)
    parser.add_argument("--min-n", type=int, default=30)
    parser.add_argument("--min-pos", type=int, default=3)
    parser.add_argument("--min-neg", type=int, default=3)
    return parser.parse_args()


def parse_background_resistant_count(signature: str) -> int:
    resistant = 0
    for part in str(signature or "").split("|"):
        if "=" not in part:
            continue
        value = part.rsplit("=", 1)[-1].strip().upper()
        resistant += int(value in {"R", "1", "RESISTANT"})
    return resistant


def normalized_rows(records: Iterable[dict]) -> list[dict]:
    rows = []
    for record in records:
        label = normalize_label(record.get("label"))
        prob = parse_float(record.get("prob"))
        if label is None or prob is None:
            continue
        signature = str(record.get("background_signature", "") or "NO_BACKGROUND_SIGNATURE").strip()
        burden = parse_float(record.get("background_resistant_count", ""))
        if burden is None:
            burden = float(parse_background_resistant_count(signature))
        row = dict(record)
        row["label"] = label
        row["prob"] = min(1.0, max(0.0, prob))
        row["background_signature"] = signature
        row["background_resistant_count"] = float(burden)
        rows.append(row)
    return rows


def adequacy_label(n: int, n_pos: int, n_neg: int, min_n: int, min_pos: int, min_neg: int) -> str:
    if n < min_n:
        return "underpowered_low_n"
    if n_pos < min_pos:
        return "underpowered_low_resistant"
    if n_neg < min_neg:
        return "underpowered_low_susceptible"
    return "interpretable"


def exact_background_leave_one_out_scores(group: Sequence[dict], smoothing_alpha: float) -> list[float]:
    labels = [int(row["label"]) for row in group]
    global_prev = mean(labels)
    by_signature: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(group):
        by_signature[str(row["background_signature"])].append(idx)

    scores = [global_prev] * len(group)
    alpha = max(0.0, float(smoothing_alpha))
    for indices in by_signature.values():
        n_bg = len(indices)
        pos_bg = sum(labels[idx] for idx in indices)
        for idx in indices:
            n_other = n_bg - 1
            pos_other = pos_bg - labels[idx]
            denominator = n_other + alpha
            if denominator > 0:
                scores[idx] = (pos_other + alpha * global_prev) / denominator
            else:
                scores[idx] = global_prev
    return scores


def interpretation(observed_auc: float, exact_auc: float, burden_auc: float, adequacy: str) -> str:
    if adequacy != "interpretable" or not math.isfinite(observed_auc):
        return "do_not_overinterpret"
    observed_minus_exact = observed_auc - exact_auc if math.isfinite(exact_auc) else math.nan
    observed_minus_burden = observed_auc - burden_auc if math.isfinite(burden_auc) else math.nan
    if math.isfinite(observed_minus_exact) and observed_minus_exact >= 0.05:
        return "maldi_exceeds_exact_background"
    if math.isfinite(observed_minus_burden) and observed_minus_burden >= 0.05:
        return "maldi_exceeds_burden_but_not_exact_background"
    if math.isfinite(exact_auc) and exact_auc >= observed_auc - 0.03:
        return "exact_background_competitive"
    return "background_only_weak_or_ambiguous"


def summarize_group(
    site: str,
    organism: str,
    drug: str,
    group: Sequence[dict],
    *,
    smoothing_alpha: float,
    min_n: int,
    min_pos: int,
    min_neg: int,
) -> dict:
    labels = [int(row["label"]) for row in group]
    scores = [float(row["prob"]) for row in group]
    burdens = [float(row["background_resistant_count"]) for row in group]
    exact_scores = exact_background_leave_one_out_scores(group, smoothing_alpha)
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    observed = safe_auc(labels, scores)
    burden_auc = safe_auc(labels, burdens)
    exact_auc = safe_auc(labels, exact_scores)
    adequacy = adequacy_label(len(group), n_pos, n_neg, min_n, min_pos, min_neg)
    counts = Counter(str(row["background_signature"]) for row in group)

    return {
        "site": site,
        "organism": organism,
        "drug": drug,
        "n": len(group),
        "n_r": n_pos,
        "n_s": n_neg,
        "prevalence": mean(labels),
        "observed_auc": observed,
        "background_burden_auc": burden_auc,
        "exact_background_auc": exact_auc,
        "observed_minus_burden_auc": observed - burden_auc if math.isfinite(observed) and math.isfinite(burden_auc) else math.nan,
        "observed_minus_exact_background_auc": observed - exact_auc if math.isfinite(observed) and math.isfinite(exact_auc) else math.nan,
        "exact_minus_burden_auc": exact_auc - burden_auc if math.isfinite(exact_auc) and math.isfinite(burden_auc) else math.nan,
        "score_exact_background_correlation": pearson(scores, exact_scores),
        "score_background_burden_correlation": pearson(scores, burdens),
        "unique_backgrounds": len(counts),
        "singleton_backgrounds": sum(1 for count in counts.values() if count == 1),
        "smoothing_alpha": smoothing_alpha,
        "adequacy_label": adequacy,
        "baseline_interpretation": interpretation(observed, exact_auc, burden_auc, adequacy),
    }


def build_baseline_rows(
    records: Iterable[dict],
    *,
    smoothing_alpha: float = 2.0,
    min_n: int = 30,
    min_pos: int = 3,
    min_neg: int = 3,
) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in normalized_rows(records):
        grouped[(row.get("site", ""), row.get("organism", ""), row.get("drug", ""))].append(row)

    output = []
    for (site, organism, drug), group in sorted(grouped.items()):
        output.append(
            summarize_group(
                site,
                organism,
                drug,
                group,
                smoothing_alpha=smoothing_alpha,
                min_n=min_n,
                min_pos=min_pos,
                min_neg=min_neg,
            )
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
    if len(rows) > limit:
        lines.extend(["", f"_Showing first {limit} of {len(rows)} rows. See CSV for complete output._"])
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    with Path(args.predictions_csv).open(newline="") as f:
        rows = build_baseline_rows(
            csv.DictReader(f),
            smoothing_alpha=args.smoothing_alpha,
            min_n=args.min_n,
            min_pos=args.min_pos,
            min_neg=args.min_neg,
        )
    write_csv(output_dir / "co_resistance_only_baseline.csv", rows)
    md = [
        "# Co-Resistance-Only Baseline",
        "",
        "Predicts each focal drug from the non-focal AST background signature only, without MALDI spectra.",
        "Exact-background scores use leave-one-out smoothed prevalence, so an isolate's own focal label is not used in its score.",
        "",
        f"- Prediction rows: `{args.predictions_csv}`",
        f"- Smoothing alpha: {args.smoothing_alpha}",
        f"- Minimum support: n>={args.min_n}, R>={args.min_pos}, S>={args.min_neg}",
        "",
        markdown_table(rows),
        "",
    ]
    (output_dir / "co_resistance_only_baseline.md").write_text("\n".join(md))
    print(f"Wrote {len(rows)} co-resistance-only baseline rows to {output_dir}")


if __name__ == "__main__":
    main()
