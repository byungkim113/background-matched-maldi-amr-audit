#!/usr/bin/env python3
"""Proteomic cross-reference for WGS-linked Basel UPEC MALDI peaks.

This script is intentionally separate from the main validation script. It uses
the public WGS-linked UPEC metadata plus Bruker median-peak table to identify
discriminative MALDI peak bins for lineage/resistance targets and checks whether
those bins overlap published LC-MS/MS-annotated ST131 MALDI biomarkers.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter
from pathlib import Path
from statistics import mean, pstdev


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "data_manifests" / "upec_master_metadata.tsv"
DEFAULT_MEDIAN_PEAKS = ROOT / "data_manifests" / "Bruker_csv_medianpeaks_df.csv"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "analysis_outputs" / "updated_proteomic_overlap_outputs"

PUBLISHED_ST131_PEAKS = [
    {
        "published_mz": 3236,
        "protein": "HdeA multivalent ion",
        "annotation": "multivalent ion of HdeA m/z 9710",
    },
    {
        "published_mz": 4176,
        "protein": "YjbJ multivalent ion",
        "annotation": "multivalent ion of YjbJ m/z 8351",
    },
    {
        "published_mz": 4857,
        "protein": "HdeA multivalent ion",
        "annotation": "multivalent ion of HdeA m/z 9710",
    },
    {
        "published_mz": 5381,
        "protein": "unidentified",
        "annotation": "ST131-specific peak not identified by LC-MS/MS",
    },
    {
        "published_mz": 6827,
        "protein": "unidentified",
        "annotation": "ST131-specific peak not identified by LC-MS/MS",
    },
    {
        "published_mz": 7655,
        "protein": "YahO",
        "annotation": "uncharacterized protein YahO",
    },
    {
        "published_mz": 8351,
        "protein": "YjbJ",
        "annotation": "UPF0337 protein YjbJ",
    },
    {
        "published_mz": 8448,
        "protein": "YnfD",
        "annotation": "uncharacterized protein YnfD",
    },
    {
        "published_mz": 9710,
        "protein": "HdeA",
        "annotation": "acid stress chaperone HdeA",
    },
    {
        "published_mz": 11783,
        "protein": "cytochrome b562",
        "annotation": "soluble cytochrome b562",
    },
]
LITERATURE_SOURCE = "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z"

# Confidence tiers for mass-matching based on MALDI-TOF linear-mode accuracy.
# At 4–12 kDa, linear mode is typically ±0.1–0.3% (~5–25 Da worst case).
# ≤10 Da: within tight accuracy window — confirmed match.
# 10–20 Da: borderline, within worst-case linear accuracy — putative.
# >20 Da: exceeds typical accuracy; enrichment test remains valid but individual
#          protein identities should not be claimed.
CONFIDENCE_THRESHOLDS = [
    (10.0, "high_confidence"),
    (20.0, "putative"),
    (float("inf"), "loose"),
]


def assign_confidence_tier(delta_da: float) -> str:
    for threshold, label in CONFIDENCE_THRESHOLDS:
        if delta_da <= threshold:
            return label
    return "loose"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: str) -> int | None:
    text = str(value).strip().lower()
    if text == "true":
        return 1
    if text == "false":
        return 0
    return None


def normalize_ast_call(value: str) -> int | None:
    text = str(value).strip().lower()
    if text in {"sensibel", "susceptible", "s"}:
        return 0
    if text in {"resistent", "resistant", "r"}:
        return 1
    return None


def get_binary_labels(
    metadata_rows: list[dict[str, str]],
    column: str,
    parser,
) -> dict[str, int]:
    labels: dict[str, int] = {}
    for row in metadata_rows:
        label = parser(row.get(column, ""))
        if label is not None:
            labels[row["TGNR"]] = label
    return labels


def build_peak_features(
    median_peaks: Path,
    allowed_tgnrs: set[str],
    mz_min: float,
    mz_max: float,
    bin_width: float,
) -> tuple[list[str], list[list[float]], list[str]]:
    n_bins = int(math.ceil((mz_max - mz_min) / bin_width))
    features = {tgnr: [0.0] * n_bins for tgnr in allowed_tgnrs}
    with median_peaks.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            tgnr = row.get("TGNR", "").strip().replace("_", "-")
            if tgnr not in features:
                continue
            try:
                mass = float(row["mass"])
                intensity = float(row["intensity"])
            except (KeyError, ValueError):
                continue
            if mass < mz_min or mass >= mz_max or intensity <= 0:
                continue
            idx = int((mass - mz_min) // bin_width)
            features[tgnr][idx] += intensity

    tgnrs: list[str] = []
    matrix: list[list[float]] = []
    for tgnr in sorted(features):
        row = features[tgnr]
        total = sum(row)
        if total <= 0:
            continue
        tgnrs.append(tgnr)
        matrix.append([math.log1p(value / total * 1_000_000.0) for value in row])

    names = [f"mz_{mz_min + i * bin_width:.0f}_{mz_min + (i + 1) * bin_width:.0f}" for i in range(n_bins)]
    return tgnrs, matrix, names


def feature_center(feature_name: str) -> float:
    _, low, high = feature_name.split("_")
    return (float(low) + float(high)) / 2.0


def match_published_peak(mz: float, tolerance_da: float) -> dict[str, object]:
    best = min(PUBLISHED_ST131_PEAKS, key=lambda row: abs(float(row["published_mz"]) - mz))
    delta = abs(float(best["published_mz"]) - mz)
    if delta > tolerance_da:
        return {"published_mz": "", "protein": "", "annotation": "", "delta_da": "", "source": ""}
    return {
        "published_mz": best["published_mz"],
        "protein": best["protein"],
        "annotation": best["annotation"],
        "delta_da": round(delta, 4),
        "confidence_tier": assign_confidence_tier(delta),
        "source": LITERATURE_SOURCE,
    }


def has_published_peak_match(mz: float, tolerance_da: float) -> bool:
    return match_published_peak(mz, tolerance_da)["published_mz"] != ""


def mz_stratum(mz: float, mz_min: float, mz_max: float, n_strata: int) -> int:
    if n_strata <= 1:
        return 0
    width = (mz_max - mz_min) / n_strata
    if width <= 0:
        return 0
    return max(0, min(n_strata - 1, int((mz - mz_min) // width)))


def sample_mass_matched_centers(
    rng: random.Random,
    observed_centers: list[float],
    feature_centers: list[float],
    mz_min: float,
    mz_max: float,
    n_strata: int,
) -> list[float]:
    centers_by_stratum: dict[int, list[float]] = {}
    for center in feature_centers:
        centers_by_stratum.setdefault(mz_stratum(center, mz_min, mz_max, n_strata), []).append(center)

    observed_by_stratum = Counter(mz_stratum(center, mz_min, mz_max, n_strata) for center in observed_centers)
    sampled: list[float] = []
    for stratum, count in observed_by_stratum.items():
        candidates = centers_by_stratum.get(stratum, feature_centers)
        if len(candidates) >= count:
            sampled.extend(rng.sample(candidates, count))
        else:
            sampled.extend(rng.choice(candidates) for _ in range(count))
    return sampled


def permutation_enrichment(
    target: str,
    observed_peak_rows: list[dict[str, object]],
    feature_names: list[str],
    tolerance_da: float,
    permutations: int,
    seed: int,
    mz_min: float,
    mz_max: float,
    n_strata: int,
) -> dict[str, object]:
    observed_centers = [float(row["mz_center"]) for row in observed_peak_rows]
    feature_centers = [feature_center(name) for name in feature_names]
    observed_count = sum(has_published_peak_match(center, tolerance_da) for center in observed_centers)

    rng = random.Random(seed)
    null_counts: list[int] = []
    for _ in range(permutations):
        sampled_centers = sample_mass_matched_centers(
            rng,
            observed_centers,
            feature_centers,
            mz_min,
            mz_max,
            n_strata,
        )
        null_counts.append(sum(has_published_peak_match(center, tolerance_da) for center in sampled_centers))

    expected = mean(null_counts) if null_counts else 0.0
    sd = pstdev(null_counts) if len(null_counts) > 1 else 0.0
    p_ge = (1 + sum(count >= observed_count for count in null_counts)) / (len(null_counts) + 1)
    fold = (observed_count / expected) if expected > 0 else ""
    z_score = ((observed_count - expected) / sd) if sd > 0 else ""
    return {
        "target": target,
        "observed_overlap_count": observed_count,
        "top_peak_count": len(observed_centers),
        "null_mean_overlap_count": round(expected, 6),
        "null_sd_overlap_count": round(sd, 6),
        "fold_enrichment": round(fold, 6) if isinstance(fold, float) else fold,
        "z_score": round(z_score, 6) if isinstance(z_score, float) else z_score,
        "empirical_p_ge_observed": round(p_ge, 6),
        "permutations": permutations,
        "mz_strata": n_strata,
        "tolerance_da": tolerance_da,
    }


def peak_effect_table(
    tgnrs: list[str],
    matrix: list[list[float]],
    feature_names: list[str],
    labels: dict[str, int],
    target: str,
    top_n: int,
) -> list[dict[str, object]]:
    usable = [(row, labels[tgnr]) for tgnr, row in zip(tgnrs, matrix) if tgnr in labels]
    counts = Counter(label for _, label in usable)
    if len(counts) != 2:
        return []
    pos = [row for row, label in usable if label == 1]
    neg = [row for row, label in usable if label == 0]
    rows: list[dict[str, object]] = []
    for idx, feature in enumerate(feature_names):
        pos_mean = mean(row[idx] for row in pos)
        neg_mean = mean(row[idx] for row in neg)
        effect = pos_mean - neg_mean
        rows.append(
            {
                "target": target,
                "feature": feature,
                "mz_center": round(feature_center(feature), 4),
                "n_pos": len(pos),
                "n_neg": len(neg),
                "mean_pos": round(pos_mean, 6),
                "mean_neg": round(neg_mean, 6),
                "effect_pos_minus_neg": round(effect, 6),
                "abs_effect": round(abs(effect), 6),
            }
        )
    rows.sort(key=lambda row: float(row["abs_effect"]), reverse=True)
    return rows[:top_n]


def run(args: argparse.Namespace) -> None:
    metadata = read_tsv(args.metadata)
    bruker_rows = [row for row in metadata if row.get("has_bruker_maldi") == "True"]
    allowed_tgnrs = {row["TGNR"] for row in bruker_rows}
    tgnrs, matrix, feature_names = build_peak_features(
        args.median_peaks,
        allowed_tgnrs,
        args.mz_min,
        args.mz_max,
        args.bin_width,
    )

    targets = {
        "ST131": get_binary_labels(metadata, "is_ST131", parse_bool),
        "Ciprofloxacin_R": get_binary_labels(metadata, "vitek_CIP_result", normalize_ast_call),
        "Ceftriaxone_R": get_binary_labels(metadata, "vitek_CTR_result", normalize_ast_call),
    }

    peak_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    rows_by_target: dict[str, list[dict[str, object]]] = {}
    for target, labels in targets.items():
        rows = peak_effect_table(tgnrs, matrix, feature_names, labels, target, args.top_peaks)
        rows_by_target[target] = rows
        peak_rows.extend(rows)
        for row in rows:
            match = match_published_peak(float(row["mz_center"]), args.tolerance_da)
            if match["published_mz"] != "":
                overlap_rows.append(
                    {
                        "target": target,
                        "feature": row["feature"],
                        "mz_center": row["mz_center"],
                        "effect_pos_minus_neg": row["effect_pos_minus_neg"],
                        "published_mz": match["published_mz"],
                        "delta_da": match["delta_da"],
                        "confidence_tier": match["confidence_tier"],
                        "protein": match["protein"],
                        "annotation": match["annotation"],
                        "source": match["source"],
                    }
                )

    enrichment_rows = [
        permutation_enrichment(
            target,
            rows,
            feature_names,
            args.tolerance_da,
            args.permutations,
            args.permutation_seed + idx,
            args.mz_min,
            args.mz_max,
            args.mz_strata,
        )
        for idx, (target, rows) in enumerate(rows_by_target.items())
    ]
    all_rows = [row for rows in rows_by_target.values() for row in rows]
    enrichment_rows.append(
        permutation_enrichment(
            "ALL_TARGETS",
            all_rows,
            feature_names,
            args.tolerance_da,
            args.permutations,
            args.permutation_seed + 10_000,
            args.mz_min,
            args.mz_max,
            args.mz_strata,
        )
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.output_dir / "updated_top_discriminative_peak_bins.csv",
        peak_rows,
        [
            "target",
            "feature",
            "mz_center",
            "n_pos",
            "n_neg",
            "mean_pos",
            "mean_neg",
            "effect_pos_minus_neg",
            "abs_effect",
        ],
    )
    write_csv(
        args.output_dir / "updated_published_st131_proteomic_overlap.csv",
        overlap_rows,
        [
            "target",
            "feature",
            "mz_center",
            "effect_pos_minus_neg",
            "published_mz",
            "delta_da",
            "confidence_tier",
            "protein",
            "annotation",
            "source",
        ],
    )
    write_csv(
        args.output_dir / "updated_proteomic_overlap_permutation_enrichment.csv",
        enrichment_rows,
        [
            "target",
            "observed_overlap_count",
            "top_peak_count",
            "null_mean_overlap_count",
            "null_sd_overlap_count",
            "fold_enrichment",
            "z_score",
            "empirical_p_ge_observed",
            "permutations",
            "mz_strata",
            "tolerance_da",
        ],
    )
    summary = {
        "metadata_rows": len(metadata),
        "bruker_maldi_rows": len(bruker_rows),
        "tgnrs_with_peak_features": len(tgnrs),
        "top_peaks_per_target": args.top_peaks,
        "bin_width": args.bin_width,
        "tolerance_da": args.tolerance_da,
        "permutations": args.permutations,
        "mz_strata": args.mz_strata,
        "overlap_rows": len(overlap_rows),
        "overlap_by_target": dict(Counter(row["target"] for row in overlap_rows)),
        "permutation_enrichment": enrichment_rows,
    }
    (args.output_dir / "updated_proteomic_overlap_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    tier_order = ["high_confidence", "putative", "loose"]
    tier_labels = {
        "high_confidence": "High confidence (≤10 Da): within typical MALDI-TOF linear-mode accuracy",
        "putative": "Putative (10–20 Da): borderline, individual protein identity uncertain",
        "loose": "Loose (>20 Da): exceeds typical accuracy — enrichment result valid but protein identity not claimed",
    }
    rows_by_tier: dict[str, list] = {t: [] for t in tier_order}
    for row in overlap_rows:
        rows_by_tier[row["confidence_tier"]].append(row)

    n_high = len(rows_by_tier["high_confidence"])
    n_putative = len(rows_by_tier["putative"])
    n_loose = len(rows_by_tier["loose"])

    lines = [
        "# Updated Proteomic Overlap Analysis",
        "",
        f"Bruker MALDI rows: {len(bruker_rows)}",
        f"TGNRs with peak features: {len(tgnrs)}",
        f"Published ST131 peak overlaps: {len(overlap_rows)} "
        f"({n_high} high-confidence, {n_putative} putative, {n_loose} loose)",
        "",
        "## Matches by Confidence Tier",
        "",
        "Mass-matching tolerance is 40 Da (=bin width). Confidence tiers reflect MALDI-TOF",
        "linear-mode accuracy at 4–12 kDa. Individual protein identities should only be",
        "claimed for high-confidence matches; the permutation enrichment result holds for all tiers.",
        "",
    ]
    for tier in tier_order:
        tier_rows = rows_by_tier[tier]
        if not tier_rows:
            continue
        lines.append(f"### {tier_labels[tier]}")
        lines.append("")
        for row in tier_rows:
            lines.append(
                f"- {row['target']} {row['feature']} center={row['mz_center']} matched "
                f"published m/z {row['published_mz']} ({row['protein']}, delta={row['delta_da']} Da)"
            )
        lines.append("")
    lines += [
        "",
        "## Mass-Matched Permutation Enrichment",
        "",
        (
            "The null model samples the same number of peak bins with the same coarse m/z-stratum "
            "distribution, then counts overlaps with the published ST131 MALDI biomarker list."
        ),
        "",
    ]
    for row in enrichment_rows:
        lines.append(
            f"- {row['target']}: observed={row['observed_overlap_count']}/"
            f"{row['top_peak_count']}, null_mean={row['null_mean_overlap_count']}, "
            f"fold={row['fold_enrichment']}, p_ge={row['empirical_p_ge_observed']}"
        )
    (args.output_dir / "updated_proteomic_overlap_summary.md").write_text("\n".join(lines) + "\n")

    print(f"Wrote outputs to {args.output_dir}")
    print(f"Published ST131 peak overlaps: {len(overlap_rows)}")
    print(f"Overlap by target: {dict(Counter(row['target'] for row in overlap_rows))}")
    print("Permutation enrichment:")
    for row in enrichment_rows:
        print(
            f"  {row['target']}: observed={row['observed_overlap_count']} "
            f"null_mean={row['null_mean_overlap_count']} "
            f"fold={row['fold_enrichment']} p={row['empirical_p_ge_observed']}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--median-peaks", type=Path, default=DEFAULT_MEDIAN_PEAKS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--mz-min", type=float, default=2000.0)
    parser.add_argument("--mz-max", type=float, default=20000.0)
    parser.add_argument("--bin-width", type=float, default=25.0)
    parser.add_argument("--top-peaks", type=int, default=75)
    parser.add_argument("--tolerance-da", type=float, default=40.0)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--permutation-seed", type=int, default=13)
    parser.add_argument("--mz-strata", type=int, default=10)
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
