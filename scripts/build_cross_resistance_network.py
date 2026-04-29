#!/usr/bin/env python3
"""Build co-resistance network artifacts from background-audit predictions.

The background-matched audit prediction CSV has one row per isolate, site, and
focal drug. The focal row contains the S/R label for that drug, so pivoting those
rows reconstructs a compact AST matrix for the evaluated panel. This script turns
that matrix into drug prevalences, pairwise co-resistance edges, and simple SVG
figures that can be used as paper/table inputs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_INPUT = (
    "/Users/byungkim/Downloads/results-2/background_matched_contrastive/"
    "background_matched_predictions.csv"
)
DEFAULT_OUTPUT_DIR = "analysis_outputs/cross_resistance_network"

DRUG_ORDER = [
    "Ciprofloxacin",
    "Norfloxacin",
    "Amoxicillin-Clavulanic acid",
    "Ceftriaxone",
    "Ceftazidime",
    "Cefepime",
]

SHORT = {
    "Ciprofloxacin": "Cipro",
    "Norfloxacin": "Norflox",
    "Amoxicillin-Clavulanic acid": "Amox-Clav",
    "Ceftriaxone": "CRO",
    "Ceftazidime": "CAZ",
    "Cefepime": "FEP",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions-csv", default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--min-edge-n",
        type=int,
        default=30,
        help="Minimum isolates with both drug labels before reporting an edge.",
    )
    parser.add_argument(
        "--network-threshold",
        type=float,
        default=0.15,
        help="Absolute phi threshold for drawing an edge in the SVG network.",
    )
    return parser.parse_args()


def read_prediction_rows(path: Path) -> Tuple[List[dict], List[str]]:
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    drugs = [d for d in DRUG_ORDER if any(r.get("drug") == d for r in rows)]
    for r in rows:
        if r.get("drug") not in SHORT:
            SHORT.setdefault(r.get("drug", ""), r.get("drug", ""))
    return rows, drugs


def isolate_key(row: dict) -> Tuple[str, str, str, str]:
    return (
        row.get("site", ""),
        row.get("raw_site", ""),
        row.get("year", ""),
        row.get("uid", ""),
    )


def build_label_matrix(rows: Sequence[dict]) -> List[dict]:
    by_iso: Dict[Tuple[str, str, str, str], dict] = {}
    for r in rows:
        drug = r.get("drug", "")
        label = r.get("label", "")
        if label not in {"0", "1", 0, 1}:
            continue
        key = isolate_key(r)
        rec = by_iso.setdefault(
            key,
            {
                "site": key[0],
                "raw_site": key[1],
                "year": key[2],
                "uid": key[3],
            },
        )
        rec[drug] = int(label)
    return list(by_iso.values())


def safe_float(x: float | None) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return f"{x:.6g}"


def phi_from_counts(n00: int, n01: int, n10: int, n11: int) -> float:
    denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    if denom == 0:
        return math.nan
    return (n11 * n00 - n10 * n01) / denom


def pair_metrics(records: Sequence[dict], drug_a: str, drug_b: str) -> dict:
    n00 = n01 = n10 = n11 = 0
    for rec in records:
        if drug_a not in rec or drug_b not in rec:
            continue
        a = rec[drug_a]
        b = rec[drug_b]
        if a == 0 and b == 0:
            n00 += 1
        elif a == 0 and b == 1:
            n01 += 1
        elif a == 1 and b == 0:
            n10 += 1
        elif a == 1 and b == 1:
            n11 += 1
    n = n00 + n01 + n10 + n11
    if n == 0:
        return {}
    prev_a = (n10 + n11) / n
    prev_b = (n01 + n11) / n
    expected_rr = prev_a * prev_b
    observed_rr = n11 / n
    lift = observed_rr / expected_rr if expected_rr > 0 else math.nan
    log2_lift = math.log(lift, 2) if lift > 0 else math.nan
    odds_ratio = ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5))
    return {
        "drug_a": drug_a,
        "drug_b": drug_b,
        "n_both_known": n,
        "n_ss": n00,
        "n_sr": n01,
        "n_rs": n10,
        "n_rr": n11,
        "prevalence_a": prev_a,
        "prevalence_b": prev_b,
        "p_rr_observed": observed_rr,
        "p_rr_expected_independent": expected_rr,
        "rr_lift": lift,
        "log2_rr_lift": log2_lift,
        "odds_ratio_haldane": odds_ratio,
        "log2_odds_ratio": math.log(odds_ratio, 2),
        "phi": phi_from_counts(n00, n01, n10, n11),
        "resistant_jaccard": n11 / (n11 + n10 + n01) if (n11 + n10 + n01) else math.nan,
    }


def build_edges(records: Sequence[dict], drugs: Sequence[str], min_n: int) -> List[dict]:
    edges: List[dict] = []
    by_site: Dict[str, List[dict]] = defaultdict(list)
    by_site["ALL"] = list(records)
    for rec in records:
        by_site[rec["site"]].append(rec)
    for site, site_records in sorted(by_site.items()):
        for i, drug_a in enumerate(drugs):
            for drug_b in drugs[i + 1 :]:
                m = pair_metrics(site_records, drug_a, drug_b)
                if not m or m["n_both_known"] < min_n:
                    continue
                m = {"site": site, **m}
                edges.append(m)
    return edges


def build_prevalence(records: Sequence[dict], drugs: Sequence[str]) -> List[dict]:
    rows = []
    by_site: Dict[str, List[dict]] = defaultdict(list)
    by_site["ALL"] = list(records)
    for rec in records:
        by_site[rec["site"]].append(rec)
    for site, site_records in sorted(by_site.items()):
        for drug in drugs:
            vals = [rec[drug] for rec in site_records if drug in rec]
            if not vals:
                continue
            rows.append(
                {
                    "site": site,
                    "drug": drug,
                    "n_known": len(vals),
                    "n_resistant": sum(vals),
                    "prevalence": sum(vals) / len(vals),
                }
            )
    return rows


def write_csv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: safe_float(r.get(k)) if isinstance(r.get(k), float) else r.get(k, "") for k in fields})


def color_for_phi(phi: float) -> str:
    if math.isnan(phi):
        return "#f2f2f2"
    phi = max(-1.0, min(1.0, phi))
    if phi >= 0:
        intensity = int(255 - 130 * phi)
        return f"rgb({intensity},{intensity},255)"
    intensity = int(255 - 130 * abs(phi))
    return f"rgb(255,{intensity},{intensity})"


def write_heatmap_svg(path: Path, edges: Sequence[dict], drugs: Sequence[str], site: str) -> None:
    site_edges = {(e["drug_a"], e["drug_b"]): e for e in edges if e["site"] == site}
    cell = 64
    label_w = 132
    top = 104
    width = label_w + cell * len(drugs) + 20
    height = top + cell * len(drugs) + 36
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{label_w}" y="32" font-family="Arial" font-size="18" font-weight="700">Co-resistance phi matrix: {site}</text>',
        f'<text x="{label_w}" y="56" font-family="Arial" font-size="12" fill="#555">blue = positive co-resistance, red = inverse association</text>',
    ]
    for i, drug in enumerate(drugs):
        x = label_w + i * cell + cell / 2
        parts.append(
            f'<text x="{x}" y="92" font-family="Arial" font-size="11" text-anchor="end" transform="rotate(-45 {x} 92)">{SHORT[drug]}</text>'
        )
        y = top + i * cell + cell / 2 + 4
        parts.append(f'<text x="{label_w - 8}" y="{y}" font-family="Arial" font-size="12" text-anchor="end">{SHORT[drug]}</text>')
    for i, drug_a in enumerate(drugs):
        for j, drug_b in enumerate(drugs):
            x = label_w + j * cell
            y = top + i * cell
            if i == j:
                fill = "#e8e8e8"
                text = "1.00"
            else:
                key = (drug_a, drug_b) if drug_a < drug_b else (drug_b, drug_a)
                # Preserve the original drug ordering rather than lexical ordering.
                e = None
                for cand in ((drug_a, drug_b), (drug_b, drug_a)):
                    if cand in site_edges:
                        e = site_edges[cand]
                        break
                phi = float(e["phi"]) if e else math.nan
                fill = color_for_phi(phi)
                text = "" if math.isnan(phi) else f"{phi:.2f}"
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="#fff"/>')
            parts.append(
                f'<text x="{x + cell / 2}" y="{y + cell / 2 + 4}" font-family="Arial" font-size="13" text-anchor="middle">{text}</text>'
            )
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def write_network_svg(path: Path, edges: Sequence[dict], prevalence: Sequence[dict], drugs: Sequence[str], site: str, threshold: float) -> None:
    site_edges = [e for e in edges if e["site"] == site and abs(float(e["phi"])) >= threshold]
    prev = {(r["site"], r["drug"]): float(r["prevalence"]) for r in prevalence}
    width = 820
    height = 560
    cx = width / 2
    cy = height / 2 + 20
    radius = 190
    positions = {}
    for i, drug in enumerate(drugs):
        angle = -math.pi / 2 + 2 * math.pi * i / len(drugs)
        positions[drug] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="32" y="36" font-family="Arial" font-size="20" font-weight="700">Co-resistance network: {site}</text>',
        f'<text x="32" y="60" font-family="Arial" font-size="12" fill="#555">Edges shown when |phi| >= {threshold:.2f}; blue positive, red inverse. Node size = resistance prevalence.</text>',
    ]
    for e in site_edges:
        a, b = e["drug_a"], e["drug_b"]
        x1, y1 = positions[a]
        x2, y2 = positions[b]
        phi = float(e["phi"])
        color = "#3b5bdb" if phi >= 0 else "#d9480f"
        sw = 1.5 + 7 * min(abs(phi), 1)
        parts.append(
            f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{sw:.2f}" stroke-opacity="0.62"/>'
        )
        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        parts.append(f'<text x="{mx}" y="{my}" font-family="Arial" font-size="11" text-anchor="middle" fill="#333">{phi:.2f}</text>')
    for drug in drugs:
        x, y = positions[drug]
        p = prev.get((site, drug), math.nan)
        r = 19 if math.isnan(p) else 16 + 38 * p
        parts.append(f'<circle cx="{x}" cy="{y}" r="{r:.1f}" fill="#f8f9fa" stroke="#222" stroke-width="1.5"/>')
        parts.append(f'<text x="{x}" y="{y - 3}" font-family="Arial" font-size="13" font-weight="700" text-anchor="middle">{SHORT[drug]}</text>')
        label = "NA" if math.isnan(p) else f"R {p:.0%}"
        parts.append(f'<text x="{x}" y="{y + 13}" font-family="Arial" font-size="11" text-anchor="middle" fill="#555">{label}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def write_markdown_summary(path: Path, edges: Sequence[dict], prevalence: Sequence[dict], drugs: Sequence[str]) -> None:
    all_edges = [e for e in edges if e["site"] == "ALL"]
    top_pos = sorted(all_edges, key=lambda e: float(e["phi"]), reverse=True)[:8]
    top_lift = sorted(all_edges, key=lambda e: float(e["log2_rr_lift"]), reverse=True)[:8]
    prev_all = [r for r in prevalence if r["site"] == "ALL"]
    lines = [
        "# Cross-Resistance Network Summary",
        "",
        "Built from the same isolate/drug label matrix used by the background-matched audit.",
        "",
        "## Overall Resistance Prevalence",
        "",
        "| Drug | Known n | Resistant n | Prevalence |",
        "| --- | ---: | ---: | ---: |",
    ]
    for r in prev_all:
        lines.append(f"| {r['drug']} | {r['n_known']} | {r['n_resistant']} | {float(r['prevalence']):.3f} |")
    lines += [
        "",
        "## Strongest Positive Phi Edges",
        "",
        "| Drug A | Drug B | n | RR observed | RR expected | Lift | Phi | Resistant Jaccard |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for e in top_pos:
        lines.append(
            f"| {e['drug_a']} | {e['drug_b']} | {e['n_both_known']} | "
            f"{float(e['p_rr_observed']):.3f} | {float(e['p_rr_expected_independent']):.3f} | "
            f"{float(e['rr_lift']):.2f} | {float(e['phi']):.3f} | {float(e['resistant_jaccard']):.3f} |"
        )
    lines += [
        "",
        "## Strongest Co-Resistance Lift Edges",
        "",
        "| Drug A | Drug B | n | Lift | Phi | n RR |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for e in top_lift:
        lines.append(
            f"| {e['drug_a']} | {e['drug_b']} | {e['n_both_known']} | "
            f"{float(e['rr_lift']):.2f} | {float(e['phi']):.3f} | {e['n_rr']} |"
        )
    lines += [
        "",
        "## Paper Use",
        "",
        "Use this as the ecological/background layer underneath the audit: drugs connected by strong co-resistance edges share resistant subpopulation structure, so high raw AUC may reflect that shared background rather than focal-drug biology.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    pred_path = Path(args.predictions_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows, drugs = read_prediction_rows(pred_path)
    records = build_label_matrix(rows)
    if not records:
        raise ValueError("No isolate label records could be reconstructed.")

    prevalence = build_prevalence(records, drugs)
    edges = build_edges(records, drugs, args.min_edge_n)

    prevalence_fields = ["site", "drug", "n_known", "n_resistant", "prevalence"]
    edge_fields = [
        "site",
        "drug_a",
        "drug_b",
        "n_both_known",
        "n_ss",
        "n_sr",
        "n_rs",
        "n_rr",
        "prevalence_a",
        "prevalence_b",
        "p_rr_observed",
        "p_rr_expected_independent",
        "rr_lift",
        "log2_rr_lift",
        "odds_ratio_haldane",
        "log2_odds_ratio",
        "phi",
        "resistant_jaccard",
    ]
    write_csv(output_dir / "cross_resistance_prevalence.csv", prevalence, prevalence_fields)
    write_csv(output_dir / "cross_resistance_edges.csv", edges, edge_fields)
    write_markdown_summary(output_dir / "cross_resistance_network_summary.md", edges, prevalence, drugs)
    for site in ["ALL", "A-2018", "DRIAMS-B", "DRIAMS-C", "DRIAMS-D"]:
        if any(e["site"] == site for e in edges):
            write_heatmap_svg(output_dir / f"cross_resistance_phi_heatmap_{site}.svg", edges, drugs, site)
            write_network_svg(
                output_dir / f"cross_resistance_network_{site}.svg",
                edges,
                prevalence,
                drugs,
                site,
                args.network_threshold,
            )
    print(f"Reconstructed {len(records)} isolate label records from {pred_path}")
    print(f"Wrote {len(edges)} pairwise edges and {len(prevalence)} prevalence rows to {output_dir}")


if __name__ == "__main__":
    main()
