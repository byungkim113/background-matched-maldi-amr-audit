#!/usr/bin/env python3
"""Sensitivity analysis: vary minimum-stratum size and report audit stability.

Re-runs the background-matched audit at a range of minimum-stratum thresholds to
show that conclusions are robust to this hyperparameter choice.  The primary result
of interest is that macro_centered_auc and mean_delta are stable across thresholds;
if they change dramatically, the matching constraint may be too strict for the
available data.

Usage:
    python scripts/sensitivity_sweep.py \\
        --predictions-csv example_predictions.csv \\
        --output-dir outputs/sensitivity_sweep

    python scripts/sensitivity_sweep.py \\
        --predictions-csv example_predictions.csv \\
        --output-dir outputs/sensitivity_sweep \\
        --thresholds 2,3,5,10 \\
        --site-robustness   # also report per-site stability

Outputs:
    sensitivity_detail.csv   — one row per (threshold × site × organism × drug)
    sensitivity_summary.csv  — macro statistics per threshold
    sensitivity_summary.md   — markdown table for the paper
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run_background_audit_framework as fw


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sensitivity sweep over minimum-stratum thresholds."
    )
    p.add_argument("--predictions-csv", required=True,
                   help="Long-format predictions CSV (see SCHEMA.md).")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--thresholds", default="2,3,5,10",
                   help="Comma-separated min pos/neg per stratum values (default: 2,3,5,10).")
    p.add_argument("--id-col",       default="isolate_id")
    p.add_argument("--site-col",     default="site")
    p.add_argument("--year-col",     default="year")
    p.add_argument("--organism-col", default="organism")
    p.add_argument("--drug-col",     default="drug")
    p.add_argument("--label-col",    default="label")
    p.add_argument("--prob-col",     default="prob")
    p.add_argument("--background-drugs", default="",
                   help="Optional comma-separated drug list to use as background.")
    p.add_argument("--match-year", action="store_true",
                   help="Also match on collection year within each stratum.")
    p.add_argument("--site-robustness", action="store_true",
                   help="Add a per-site stability table to the markdown output.")
    return p.parse_args()


def _stratum_centered_auc(pair_rows: list[dict], min_t: int,
                           match_year: bool = False) -> dict:
    """Compute matched and stratum-centered AUC for one (site, organism, drug) slice."""
    labels_all = [r["label"] for r in pair_rows]
    probs_all  = [r["prob"]  for r in pair_rows]
    raw_auc    = fw.safe_auc(labels_all, probs_all)
    raw_aupr   = fw.safe_aupr(labels_all, probs_all)

    # Build strata
    def sig_key(r: dict) -> tuple:
        base = (r["background_signature"],)
        return base + ((r["year"],) if match_year else ())

    strata: dict[tuple, list[dict]] = defaultdict(list)
    for r in pair_rows:
        strata[sig_key(r)].append(r)

    valid_strata: dict[tuple, list[dict]] = {}
    for key, group in strata.items():
        n_pos = sum(1 for r in group if r["label"] == 1)
        n_neg = sum(1 for r in group if r["label"] == 0)
        if n_pos >= min_t and n_neg >= min_t:
            valid_strata[key] = group

    valid_rows = [r for rows in valid_strata.values() for r in rows]
    n_total  = len(pair_rows)
    n_r      = sum(1 for r in pair_rows if r["label"] == 1)
    n_matched   = len(valid_rows)
    n_matched_r = sum(1 for r in valid_rows if r["label"] == 1)
    matched_retention = n_matched / n_total if n_total else math.nan

    if valid_rows:
        strata_means = {
            key: sum(r["prob"] for r in rows) / len(rows)
            for key, rows in valid_strata.items()
        }
        matched_labels  = [r["label"] for r in valid_rows]
        matched_probs   = [r["prob"]  for r in valid_rows]
        centered_probs  = [r["prob"] - strata_means[sig_key(r)] for r in valid_rows]
        matched_auc     = fw.safe_auc(matched_labels, matched_probs)
        centered_auc    = fw.safe_auc(matched_labels, centered_probs)

        # Pairwise accuracy within strata
        wins = 0.0
        total_pairs = 0
        for rows in valid_strata.values():
            pos_p = [r["prob"] for r in rows if r["label"] == 1]
            neg_p = [r["prob"] for r in rows if r["label"] == 0]
            if not pos_p or not neg_p:
                continue
            for pp in pos_p:
                for np_ in neg_p:
                    total_pairs += 1
                    wins += 1.0 if pp > np_ else (0.5 if pp == np_ else 0.0)
        pairwise_acc = wins / total_pairs if total_pairs else math.nan
    else:
        matched_auc   = math.nan
        centered_auc  = math.nan
        pairwise_acc  = math.nan
        total_pairs   = 0

    return dict(
        raw_auc=raw_auc,
        raw_aupr=raw_aupr,
        matched_auc=matched_auc,
        centered_auc=centered_auc,
        pairwise_accuracy=pairwise_acc,
        pairwise_comparisons=total_pairs,
        matched_retention=matched_retention,
        n_total=n_total,
        n_r=n_r,
        n_matched=n_matched,
        n_matched_r=n_matched_r,
        n_valid_strata=len(valid_strata),
    )


def run_sweep(rows_with_sigs: list[dict], thresholds: list[int],
              match_year: bool = False) -> tuple[list[dict], list[dict]]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows_with_sigs:
        grouped[(r["site"], r["organism"], r["drug"])].append(r)

    detail_rows: list[dict] = []
    for t in thresholds:
        for (site, organism, drug), pair_rows in sorted(grouped.items()):
            result = _stratum_centered_auc(pair_rows, t, match_year=match_year)
            detail_rows.append(dict(
                min_stratum=t,
                site=site,
                organism=organism,
                drug=drug,
                raw_auc=_fmt(result["raw_auc"]),
                matched_auc=_fmt(result["matched_auc"]),
                centered_auc=_fmt(result["centered_auc"]),
                pairwise_accuracy=_fmt(result["pairwise_accuracy"]),
                pairwise_comparisons=result["pairwise_comparisons"],
                matched_retention=_fmt(result["matched_retention"]),
                n_total=result["n_total"],
                n_r=result["n_r"],
                n_matched=result["n_matched"],
                n_matched_r=result["n_matched_r"],
                n_valid_strata=result["n_valid_strata"],
            ))

    summary_rows = [_summarize(t, [r for r in detail_rows if r["min_stratum"] == t])
                    for t in thresholds]
    return detail_rows, summary_rows


def _fmt(v: float, decimals: int = 4) -> str:
    return "" if not math.isfinite(v) else f"{v:.{decimals}f}"


def _mean_str(vals: list[float]) -> str:
    finite = [v for v in vals if math.isfinite(v)]
    return _fmt(sum(finite) / len(finite)) if finite else ""


def _summarize(t: int, detail: list[dict]) -> dict:
    adequate = [r for r in detail if r["centered_auc"] != ""]
    raw_aucs  = [float(r["raw_auc"])  for r in detail    if r["raw_auc"]  != ""]
    cen_aucs  = [float(r["centered_auc"]) for r in adequate]
    retentions = [float(r["matched_retention"]) for r in detail if r["matched_retention"] != ""]
    deltas = [float(r["raw_auc"]) - float(r["centered_auc"])
              for r in adequate if r["raw_auc"] != ""]
    return dict(
        min_stratum=t,
        n_drug_site_pairs=len(detail),
        n_adequate=len(adequate),
        mean_matched_retention=_mean_str(retentions),
        macro_raw_auc=_mean_str(raw_aucs),
        macro_centered_auc=_mean_str(cen_aucs),
        mean_delta=_mean_str(deltas),
    )


def _site_robustness(detail_rows: list[dict], thresholds: list[int]) -> list[dict]:
    """Per-site summary across thresholds — checks if any site is anomalous."""
    rows = []
    sites = sorted({r["site"] for r in detail_rows})
    for site in sites:
        for t in thresholds:
            subset = [r for r in detail_rows if r["site"] == site and r["min_stratum"] == t]
            adequate = [r for r in subset if r["centered_auc"] != ""]
            raw_aucs = [float(r["raw_auc"]) for r in subset if r["raw_auc"] != ""]
            cen_aucs = [float(r["centered_auc"]) for r in adequate]
            rows.append(dict(
                site=site,
                min_stratum=t,
                n_pairs=len(subset),
                n_adequate=len(adequate),
                macro_raw_auc=_mean_str(raw_aucs),
                macro_centered_auc=_mean_str(cen_aucs),
            ))
    return rows


def _markdown_table(rows: list[dict], cols: list[str], title: str = "") -> str:
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join("---" for _ in cols) + " |"
    lines  = ([f"## {title}\n"] if title else []) + [header, sep]
    for r in rows:
        lines.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(lines)


DETAIL_COLS = [
    "min_stratum", "site", "organism", "drug",
    "raw_auc", "matched_auc", "centered_auc", "pairwise_accuracy", "pairwise_comparisons",
    "matched_retention", "n_total", "n_r", "n_matched", "n_matched_r", "n_valid_strata",
]
SUMMARY_COLS = [
    "min_stratum", "n_drug_site_pairs", "n_adequate",
    "mean_matched_retention", "macro_raw_auc", "macro_centered_auc", "mean_delta",
]
SITE_COLS = ["site", "min_stratum", "n_pairs", "n_adequate",
             "macro_raw_auc", "macro_centered_auc"]


def main() -> None:
    args = parse_args()
    thresholds = [int(t.strip()) for t in args.thresholds.split(",") if t.strip()]
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = fw.read_long_predictions(
        args.predictions_csv,
        id_col=args.id_col,
        site_col=args.site_col,
        year_col=args.year_col,
        organism_col=args.organism_col,
        drug_col=args.drug_col,
        label_col=args.label_col,
        prob_col=args.prob_col,
    )
    bg_drugs = [d.strip() for d in args.background_drugs.split(",") if d.strip()] or None
    rows_with_sigs = fw.add_background_signatures(rows, background_drugs=bg_drugs)

    detail_rows, summary_rows = run_sweep(rows_with_sigs, thresholds,
                                          match_year=args.match_year)

    with (out_dir / "sensitivity_detail.csv").open("w", newline="") as f:
        csv.DictWriter(f, fieldnames=DETAIL_COLS).writeheader()
        csv.DictWriter(f, fieldnames=DETAIL_COLS).writerows(detail_rows)

    with (out_dir / "sensitivity_summary.csv").open("w", newline="") as f:
        csv.DictWriter(f, fieldnames=SUMMARY_COLS).writeheader()
        csv.DictWriter(f, fieldnames=SUMMARY_COLS).writerows(summary_rows)

    md_parts = [
        "# Sensitivity Analysis: Minimum Stratum Size\n",
        "The audit is re-run at each minimum-stratum threshold. "
        "Stable `macro_centered_auc` and `mean_delta` across thresholds shows "
        "conclusions are not sensitive to this hyperparameter.\n",
        "\n",
        _markdown_table(summary_rows, SUMMARY_COLS, "Macro summary across thresholds"),
    ]

    if args.site_robustness:
        site_rows = _site_robustness(detail_rows, thresholds)
        with (out_dir / "sensitivity_site_robustness.csv").open("w", newline="") as f:
            csv.DictWriter(f, fieldnames=SITE_COLS).writeheader()
            csv.DictWriter(f, fieldnames=SITE_COLS).writerows(site_rows)
        md_parts += ["\n\n", _markdown_table(site_rows, SITE_COLS, "Per-site robustness")]
        print(f"Wrote {out_dir / 'sensitivity_site_robustness.csv'}")

    (out_dir / "sensitivity_summary.md").write_text("\n".join(md_parts) + "\n")

    print(f"Wrote {out_dir / 'sensitivity_detail.csv'}")
    print(f"Wrote {out_dir / 'sensitivity_summary.csv'}")
    print(f"Wrote {out_dir / 'sensitivity_summary.md'}")
    print()
    print(f"{'min':>5}  {'adequate':>10}  {'retention':>10}  {'raw_auc':>9}  {'cen_auc':>9}  {'delta':>7}")
    for r in summary_rows:
        print(f"  {r['min_stratum']:>3}  {r['n_adequate']:>3}/{r['n_drug_site_pairs']:<6}  "
              f"{r['mean_matched_retention']:>10}  "
              f"{r['macro_raw_auc']:>9}  {r['macro_centered_auc']:>9}  {r['mean_delta']:>7}")


if __name__ == "__main__":
    main()
