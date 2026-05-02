#!/usr/bin/env python3
"""Feasibility audit for MARISMa v2 AMR labels.

This metadata-first pass answers three questions before downloading or
preprocessing the full MARISMa spectrum archive:

1. Which organism-drug pairs have enough S/R labels for external validation?
2. Which DRIAMS paper pairs need MARISMa-specific drug aliases?
3. Which spectrum paths should be extracted if we proceed to model deployment?
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_OUTPUT_DIR = Path("outputs/analysis_outputs/marisma_v2_metadata_audit")


@dataclass(frozen=True)
class TargetPair:
    organism: str
    paper_drug: str
    marisma_drug: str
    relationship: str
    ecology_block: str
    prediction: str
    rationale: str


TARGET_PAIRS = [
    TargetPair(
        "Escherichia coli",
        "Ciprofloxacin",
        "Ciprofloxacin",
        "exact",
        "fluoroquinolone",
        "transfer_retained_or_partial",
        "Fluoroquinolone resistance should retain signal if conserved resistant lineages are shared.",
    ),
    TargetPair(
        "Escherichia coli",
        "Norfloxacin",
        "Norfloxacin",
        "exact",
        "fluoroquinolone",
        "transfer_retained_or_partial",
        "Norfloxacin is expected to track the same fluoroquinolone block as ciprofloxacin.",
    ),
    TargetPair(
        "Escherichia coli",
        "Levofloxacin",
        "Levofloxacin",
        "external_extension",
        "fluoroquinolone",
        "transfer_retained_or_partial",
        "Levofloxacin provides an additional MARISMa fluoroquinolone validation target.",
    ),
    TargetPair(
        "Escherichia coli",
        "Amoxicillin-Clavulanic acid",
        "Amoxicillin/Clavulanic acid",
        "spelling_alias",
        "beta-lactam/inhibitor",
        "background_sensitive_or_weak",
        "Prior DRIAMS audits suggest Amox-Clav is more ecology/background-driven.",
    ),
    TargetPair(
        "Escherichia coli",
        "Ceftriaxone",
        "Cefotaxime",
        "third_generation_cephalosporin_analog",
        "third-generation cephalosporin",
        "background_sensitive_or_mixed",
        "MARISMa has E. coli cefotaxime rather than ceftriaxone labels; use as an analogous ESBL/AmpC phenotype.",
    ),
    TargetPair(
        "Escherichia coli",
        "Ceftazidime",
        "Ceftazidime",
        "exact",
        "third-generation cephalosporin",
        "background_sensitive_or_mixed",
        "Expected to sit in the same cephalosporin co-resistance block as cefotaxime/cefepime.",
    ),
    TargetPair(
        "Escherichia coli",
        "Cefepime",
        "Cefepime",
        "exact",
        "fourth-generation cephalosporin",
        "background_sensitive_or_mixed",
        "Expected to show raw signal but potentially lower background-centered signal.",
    ),
    TargetPair(
        "Escherichia coli",
        "Cotrimoxazole",
        "Trimethoprim/Sulfamethoxazole",
        "name_alias",
        "folate-pathway/MDR block",
        "exploratory",
        "Useful exploratory E. coli MDR/background target with strong label coverage.",
    ),
    TargetPair(
        "Staphylococcus aureus",
        "Oxacillin",
        "Oxacillin",
        "exact",
        "MRSA/beta-lactam",
        "transfer_retained_or_partial",
        "MRSA-related spectral structure should be externally testable.",
    ),
    TargetPair(
        "Staphylococcus epidermidis",
        "Erythromycin",
        "Erythromycin",
        "exact",
        "macrolide/opportunistic background",
        "background_sensitive_or_weak",
        "Prior DRIAMS clinical4 results suggest weak generalization for this pair.",
    ),
]


PANEL_BY_ORGANISM = {
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
    "Klebsiella pneumoniae": [
        "Ciprofloxacin",
        "Norfloxacin",
        "Levofloxacin",
        "Amoxicillin/Clavulanic acid",
        "Cefotaxime",
        "Ceftazidime",
        "Cefepime",
        "Trimethoprim/Sulfamethoxazole",
    ],
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Audit MARISMa v2 AMR metadata before spectrum processing.")
    p.add_argument("--amr-csv", type=Path, required=True, help="Path to MARISMa v2 AMR.csv.")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--min-s", type=int, default=50)
    p.add_argument("--min-r", type=int, default=50)
    p.add_argument(
        "--min-background-labels",
        type=int,
        default=3,
        help="Minimum number of other S/R labels needed to consider a row useful for background matching.",
    )
    return p


def clean_label(series: pd.Series) -> pd.Series:
    labels = series.astype("string").str.strip().str.upper()
    labels = labels.replace({"": pd.NA, "NAN": pd.NA, "NA": pd.NA})
    return labels


def sr_counts(labels: pd.Series) -> dict[str, int]:
    vc = clean_label(labels).value_counts(dropna=True)
    return {
        "n_s": int(vc.get("S", 0)),
        "n_i": int(vc.get("I", 0)),
        "n_r": int(vc.get("R", 0) + vc.get("R*", 0) + vc.get("ESBL", 0)),
        "n_r_plain": int(vc.get("R", 0)),
        "n_r_star": int(vc.get("R*", 0)),
        "n_esbl": int(vc.get("ESBL", 0)),
        "n_other_nonmissing": int(vc.sum() - vc.get("S", 0) - vc.get("I", 0) - vc.get("R", 0) - vc.get("R*", 0) - vc.get("ESBL", 0)),
    }


def binary_label(series: pd.Series) -> pd.Series:
    labels = clean_label(series)
    return labels.map({"S": 0, "R": 1, "R*": 1, "ESBL": 1})


def phi_from_counts(n00: int, n01: int, n10: int, n11: int) -> float:
    denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    if denom == 0:
        return math.nan
    return (n11 * n00 - n10 * n01) / denom


def co_resistance_edges(df: pd.DataFrame, organism: str, drugs: list[str], min_n: int = 30) -> pd.DataFrame:
    sub = df[df["Species"].eq(organism)].copy()
    available = [d for d in drugs if d in sub.columns]
    rows: list[dict] = []
    for drug_a, drug_b in combinations(available, 2):
        a = binary_label(sub[drug_a])
        b = binary_label(sub[drug_b])
        valid = a.notna() & b.notna()
        if int(valid.sum()) < min_n:
            continue
        av = a[valid].astype(int)
        bv = b[valid].astype(int)
        n00 = int(((av == 0) & (bv == 0)).sum())
        n01 = int(((av == 0) & (bv == 1)).sum())
        n10 = int(((av == 1) & (bv == 0)).sum())
        n11 = int(((av == 1) & (bv == 1)).sum())
        n = n00 + n01 + n10 + n11
        prev_a = (n10 + n11) / n
        prev_b = (n01 + n11) / n
        expected_rr = prev_a * prev_b
        observed_rr = n11 / n
        rows.append(
            {
                "organism": organism,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "n_both_known": n,
                "n_rr": n11,
                "prevalence_a": prev_a,
                "prevalence_b": prev_b,
                "rr_observed": observed_rr,
                "rr_expected_independent": expected_rr,
                "rr_lift": observed_rr / expected_rr if expected_rr else math.nan,
                "phi": phi_from_counts(n00, n01, n10, n11),
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["organism", "phi"], ascending=[True, False])


def count_background_labels(sub: pd.DataFrame, focal: str, panel: list[str]) -> pd.Series:
    others = [d for d in panel if d != focal and d in sub.columns]
    if not others:
        return pd.Series(0, index=sub.index)
    sr = pd.DataFrame({d: binary_label(sub[d]) for d in others})
    return sr.notna().sum(axis=1)


def feasibility_table(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    rows: list[dict] = []
    for target in TARGET_PAIRS:
        if target.marisma_drug not in df.columns:
            rows.append(
                {
                    "organism": target.organism,
                    "paper_drug": target.paper_drug,
                    "marisma_drug": target.marisma_drug,
                    "relationship": target.relationship,
                    "ecology_block": target.ecology_block,
                    "prediction": target.prediction,
                    "n_rows": 0,
                    "n_s": 0,
                    "n_i": 0,
                    "n_r": 0,
                    "n_sr": 0,
                    "r_prevalence": math.nan,
                    "n_with_min_background_labels": 0,
                    "eligible_metadata_only": False,
                    "rationale": target.rationale,
                }
            )
            continue
        sub = df[df["Species"].eq(target.organism)].copy()
        counts = sr_counts(sub[target.marisma_drug])
        n_sr = counts["n_s"] + counts["n_r"]
        panel = PANEL_BY_ORGANISM.get(target.organism, [])
        focal_binary = binary_label(sub[target.marisma_drug]).notna()
        bg_counts = count_background_labels(sub, target.marisma_drug, panel)
        n_bg = int((focal_binary & (bg_counts >= args.min_background_labels)).sum())
        rows.append(
            {
                "organism": target.organism,
                "paper_drug": target.paper_drug,
                "marisma_drug": target.marisma_drug,
                "relationship": target.relationship,
                "ecology_block": target.ecology_block,
                "prediction": target.prediction,
                "n_rows": int(len(sub)),
                **counts,
                "n_sr": n_sr,
                "r_prevalence": counts["n_r"] / n_sr if n_sr else math.nan,
                "n_with_min_background_labels": n_bg,
                "eligible_metadata_only": bool(
                    counts["n_s"] >= args.min_s and counts["n_r"] >= args.min_r and n_bg >= args.min_s + args.min_r
                ),
                "rationale": target.rationale,
            }
        )
    return pd.DataFrame(rows)


def spectrum_manifest(df: pd.DataFrame, feasible: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for rec in feasible[feasible["eligible_metadata_only"]].to_dict("records"):
        organism = rec["organism"]
        drug = rec["marisma_drug"]
        sub = df[df["Species"].eq(organism)].copy()
        labels = binary_label(sub[drug])
        keep = labels.notna()
        if not keep.any():
            continue
        out = sub.loc[keep, ["Identifier", "target_position", "Year", "Path", "Species", "Microorganism", "Sample"]].copy()
        out["paper_drug"] = rec["paper_drug"]
        out["marisma_drug"] = drug
        out["label"] = labels.loc[keep].astype(int).values
        rows.append(out)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).drop_duplicates()


def markdown_table(rows: Iterable[dict] | pd.DataFrame, columns: list[str]) -> str:
    """Render a small markdown table without pandas' optional tabulate dependency."""
    if isinstance(rows, pd.DataFrame):
        records = rows[columns].to_dict("records") if not rows.empty else []
    else:
        records = list(rows)
    header = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = []
    for rec in records:
        values = []
        for col in columns:
            val = rec.get(col, "")
            if isinstance(val, float):
                val = "" if math.isnan(val) else f"{val:.3f}"
            values.append(str(val))
        body.append("| " + " | ".join(values) + " |")
    return "\n".join([header, divider, *body])


def write_markdown(
    output_dir: Path,
    feasible: pd.DataFrame,
    edges: pd.DataFrame,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    eligible = feasible[feasible["eligible_metadata_only"]].copy()
    top_edges = edges.head(20).copy() if not edges.empty else edges
    unique_paths = int(manifest["Path"].nunique()) if not manifest.empty and "Path" in manifest.columns else 0
    unique_identifiers = int(manifest["Identifier"].nunique()) if not manifest.empty and "Identifier" in manifest.columns else 0
    lines = [
        "# MARISMa v2 Metadata Feasibility Audit",
        "",
        f"Input: `{args.amr_csv}`",
        "",
        "## Main Takeaways",
        "",
        f"- Eligible organism-drug targets by metadata: **{len(eligible)} / {len(feasible)}**.",
        f"- Target-level manifest rows: **{len(manifest):,}**.",
        f"- Unique spectrum paths to extract for all eligible targets: **{unique_paths:,}**.",
        f"- Unique isolate identifiers represented in selected targets: **{unique_identifiers:,}**.",
        "- `E. coli / Ceftriaxone` is mapped to MARISMa `E. coli / Cefotaxime` as an analogous third-generation cephalosporin phenotype.",
        "- This is a metadata-only audit; model AUC requires extracting/preprocessing spectra for the manifest rows.",
        "",
        "## Feasibility Table",
        "",
        markdown_table(
            feasible,
            [
                "organism",
                "paper_drug",
                "marisma_drug",
                "relationship",
                "n_s",
                "n_i",
                "n_r",
                "n_sr",
                "r_prevalence",
                "n_with_min_background_labels",
                "eligible_metadata_only",
                "prediction",
            ],
        ),
        "",
        "## Top Co-Resistance Edges",
        "",
        markdown_table(
            top_edges,
            [
                "organism",
                "drug_a",
                "drug_b",
                "n_both_known",
                "n_rr",
                "rr_lift",
                "phi",
            ],
        )
        if not top_edges.empty
        else "_No reportable edges._",
        "",
        "## Next Step",
        "",
        "Extract MARISMa spectra listed in `marisma_spectrum_manifest_for_selected_targets.csv`, preprocess them to the DRIAMS 6000-bin representation, export model predictions, and run the background-matched audit.",
    ]
    (output_dir / "marisma_v2_metadata_feasibility_audit.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.amr_csv, low_memory=False)
    feasible = feasibility_table(df, args)

    edge_frames: list[pd.DataFrame] = []
    for organism, panel in PANEL_BY_ORGANISM.items():
        if organism in set(df["Species"].dropna()):
            edge_frames.append(co_resistance_edges(df, organism, panel))
    edges = pd.concat([e for e in edge_frames if not e.empty], ignore_index=True) if edge_frames else pd.DataFrame()
    if not edges.empty:
        edges = edges.sort_values(["phi", "rr_lift"], ascending=[False, False])

    manifest = spectrum_manifest(df, feasible)
    unique_paths = (
        manifest[["Path", "Identifier", "Year", "Species", "Microorganism"]]
        .drop_duplicates(subset=["Path"])
        .sort_values(["Species", "Year", "Path"])
        if not manifest.empty
        else pd.DataFrame(columns=["Path", "Identifier", "Year", "Species", "Microorganism"])
    )

    feasible.to_csv(args.output_dir / "marisma_pair_feasibility.csv", index=False)
    edges.to_csv(args.output_dir / "marisma_cross_resistance_edges.csv", index=False)
    manifest.to_csv(args.output_dir / "marisma_spectrum_manifest_for_selected_targets.csv", index=False)
    unique_paths.to_csv(args.output_dir / "marisma_unique_spectrum_paths_for_selected_targets.csv", index=False)

    summary = {
        "input_csv": str(args.amr_csv),
        "n_rows": int(len(df)),
        "n_unique_identifiers": int(df["Identifier"].nunique()) if "Identifier" in df.columns else None,
        "n_unique_paths": int(df["Path"].nunique()) if "Path" in df.columns else None,
        "n_eligible_targets": int(feasible["eligible_metadata_only"].sum()),
        "n_manifest_rows": int(len(manifest)),
        "n_manifest_unique_paths": int(manifest["Path"].nunique()) if not manifest.empty and "Path" in manifest.columns else 0,
        "n_manifest_unique_identifiers": int(manifest["Identifier"].nunique()) if not manifest.empty and "Identifier" in manifest.columns else 0,
        "min_s": args.min_s,
        "min_r": args.min_r,
        "min_background_labels": args.min_background_labels,
    }
    (args.output_dir / "marisma_metadata_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_markdown(args.output_dir, feasible, edges, manifest, args)

    print(f"Wrote MARISMa metadata audit to {args.output_dir}")
    print(feasible[["organism", "paper_drug", "marisma_drug", "n_s", "n_i", "n_r", "n_with_min_background_labels", "eligible_metadata_only"]].to_string(index=False))


if __name__ == "__main__":
    main()
