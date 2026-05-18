#!/usr/bin/env python3
"""Revision-specific audit checks for the manuscript critique.

This script answers two reviewer-facing questions that are not visible from the
main summary tables:

1. For ciprofloxacin, how many valid background strata contain norfloxacin-R,
   norfloxacin-S, or norfloxacin-unknown backgrounds?
2. Does the primary drug contrast persist when valid strata require at least
   5 or 10 resistant and susceptible isolates rather than 3?
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run_background_audit_framework as audit


DEFAULT_INPUT = Path("/Users/byungkim/Downloads/results-2/background_matched_contrastive/background_matched_predictions.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run manuscript revision audit checks.")
    parser.add_argument("--predictions-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/analysis_outputs/revision_checks"))
    parser.add_argument("--thresholds", default="3,5,10")
    return parser.parse_args()


def parse_thresholds(text: str) -> list[int]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    return sorted(set(values))


def signature_value(signature: str, drug: str) -> str:
    prefix = f"{drug}="
    for part in str(signature).split("|"):
        if part.startswith(prefix):
            value = part[len(prefix) :].strip()
            return value or "missing"
    return "absent"


def write_csv(path: Path, rows: Iterable[dict]) -> None:
    rows = list(rows)
    fields = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def cipro_norfloxacin_composition(rows: list[dict], thresholds: list[int]) -> tuple[list[dict], list[dict]]:
    cipro = [row for row in rows if row["organism"] == "Escherichia coli" and row["drug"] == "Ciprofloxacin"]
    by_site_signature = defaultdict(list)
    for row in cipro:
        by_site_signature[(row["site"], row["background_signature"])].append(row)

    detail_rows = []
    summary_rows = []
    for threshold in thresholds:
        summary = defaultdict(Counter)
        for (site, signature), group in sorted(by_site_signature.items()):
            n_r = sum(int(row["label"]) == 1 for row in group)
            n_s = sum(int(row["label"]) == 0 for row in group)
            if n_r < threshold or n_s < threshold:
                continue
            norfloxacin_state = signature_value(signature, "Norfloxacin")
            n_total = len(group)
            detail_rows.append(
                {
                    "threshold": threshold,
                    "site": site,
                    "background_signature": signature,
                    "n_total": n_total,
                    "n_cipro_r": n_r,
                    "n_cipro_s": n_s,
                    "norfloxacin_background": norfloxacin_state,
                }
            )
            summary[(threshold, site)]["valid_strata"] += 1
            summary[(threshold, site)]["n_matched"] += n_total
            summary[(threshold, site)][f"strata_norfloxacin_{norfloxacin_state}"] += 1
            summary[(threshold, site)][f"rows_norfloxacin_{norfloxacin_state}"] += n_total
            summary[(threshold, site)][f"cipro_r_norfloxacin_{norfloxacin_state}"] += n_r
            summary[(threshold, site)][f"cipro_s_norfloxacin_{norfloxacin_state}"] += n_s

        for (threshold_value, site), counts in sorted(summary.items()):
            valid = counts["valid_strata"]
            matched = counts["n_matched"]
            row = {
                "threshold": threshold_value,
                "site": site,
                "valid_strata": valid,
                "n_matched": matched,
                "strata_norfloxacin_S": counts["strata_norfloxacin_S"],
                "strata_norfloxacin_R": counts["strata_norfloxacin_R"],
                "strata_norfloxacin_U": counts["strata_norfloxacin_U"],
                "strata_norfloxacin_absent": counts["strata_norfloxacin_absent"],
                "pct_strata_norfloxacin_S": counts["strata_norfloxacin_S"] / valid if valid else math.nan,
                "pct_strata_norfloxacin_R": counts["strata_norfloxacin_R"] / valid if valid else math.nan,
                "rows_norfloxacin_S": counts["rows_norfloxacin_S"],
                "rows_norfloxacin_R": counts["rows_norfloxacin_R"],
                "rows_norfloxacin_U": counts["rows_norfloxacin_U"],
                "rows_norfloxacin_absent": counts["rows_norfloxacin_absent"],
                "pct_rows_norfloxacin_S": counts["rows_norfloxacin_S"] / matched if matched else math.nan,
                "pct_rows_norfloxacin_R": counts["rows_norfloxacin_R"] / matched if matched else math.nan,
            }
            summary_rows.append(row)
    return summary_rows, detail_rows


def threshold_sensitivity(rows: list[dict], thresholds: list[int]) -> list[dict]:
    output = []
    for threshold in thresholds:
        summary, _ = audit.compute_background_matched_summary(
            rows,
            min_pos_per_stratum=threshold,
            min_neg_per_stratum=threshold,
            bootstrap_n=0,
            permutation_n=0,
        )
        for row in summary:
            output.append({"threshold": threshold, **row})
    return output


def write_summary(path: Path, cipro_rows: list[dict], sensitivity_rows: list[dict]) -> None:
    focus_drugs = {"Ciprofloxacin", "Amoxicillin-Clavulanic acid"}
    lines = [
        "# Revision Audit Checks",
        "",
        "## Ciprofloxacin/Norfloxacin Background Composition",
        "",
        "| Threshold | Site | Valid strata | n matched | % strata Norflox-S | % rows Norflox-S | % strata Norflox-R | % rows Norflox-R |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in cipro_rows:
        lines.append(
            f"| {row['threshold']} | {row['site']} | {row['valid_strata']} | {row['n_matched']} | "
            f"{100 * row['pct_strata_norfloxacin_S']:.1f} | {100 * row['pct_rows_norfloxacin_S']:.1f} | "
            f"{100 * row['pct_strata_norfloxacin_R']:.1f} | {100 * row['pct_rows_norfloxacin_R']:.1f} |"
        )
    lines.extend(
        [
            "",
            "## Threshold Sensitivity For Primary Contrast",
            "",
            "| Threshold | Site | Drug | Raw AUC | Centered AUC | Retention | n matched | Valid strata | Adequacy |",
            "| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in sensitivity_rows:
        if row["organism"] != "Escherichia coli" or row["drug"] not in focus_drugs:
            continue
        centered = row["stratum_centered_auc"]
        centered_text = "" if math.isnan(centered) else f"{centered:.3f}"
        lines.append(
            f"| {row['threshold']} | {row['site']} | {row['drug']} | {row['raw_auc']:.3f} | "
            f"{centered_text} | {row['matched_retention']:.3f} | {row['n_matched']} | "
            f"{row['n_valid_strata']} | {row['adequacy_label']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = parse_thresholds(args.thresholds)
    rows = audit.read_long_predictions(
        args.predictions_csv,
        id_col="uid",
        site_col="site",
        year_col="year",
        organism_col="organism",
        drug_col="drug",
        label_col="label",
        prob_col="prob",
        background_signature_col="background_signature",
        model_name="Mega/CNN",
    )
    rows = audit.add_background_signatures(rows)

    cipro_summary, cipro_detail = cipro_norfloxacin_composition(rows, thresholds)
    sensitivity = threshold_sensitivity(rows, thresholds)

    write_csv(args.output_dir / "cipro_norfloxacin_stratum_composition.csv", cipro_summary)
    write_csv(args.output_dir / "cipro_norfloxacin_valid_strata_detail.csv", cipro_detail)
    write_csv(args.output_dir / "primary_threshold_sensitivity.csv", sensitivity)
    write_summary(args.output_dir / "revision_audit_checks.md", cipro_summary, sensitivity)
    print(f"Wrote revision checks to {args.output_dir}")


if __name__ == "__main__":
    main()
