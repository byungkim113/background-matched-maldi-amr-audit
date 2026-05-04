#!/usr/bin/env python3
"""Clone-controlled UPEC MALDI/WGS analysis.

This script asks whether MALDI peak features predict ciprofloxacin or
ceftriaxone resistance after controlling for WGS-derived lineage labels in the
public Basel/Cuenod UPEC cohort. It complements the DRIAMS background-matched
audit by replacing AST-derived co-resistance background with observed WGS
background: ST131 status, exact sequence type, and phylogroup.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "data_manifests" / "upec_master_metadata.tsv"
DEFAULT_MEDIAN_PEAKS = ROOT / "data_manifests" / "Bruker_csv_medianpeaks_df.csv"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "analysis_outputs" / "upec_clone_control_outputs"


def ast_to_binary(value: object) -> float:
    text = str(value).strip().lower()
    if text in {"sensibel", "susceptible", "s"}:
        return 0.0
    if text in {"resistent", "resistant", "r"}:
        return 1.0
    return np.nan


def bool_to_binary(value: object) -> float:
    text = str(value).strip().lower()
    if text == "true":
        return 1.0
    if text == "false":
        return 0.0
    return np.nan


def auc_ci(y: np.ndarray, scores: np.ndarray, n_boot: int, seed: int) -> tuple[float, float, float]:
    auc = float(binary_auc(y, scores))
    if n_boot <= 0:
        return auc, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = []
    indices = np.arange(len(y))
    for _ in range(n_boot):
        sample = rng.choice(indices, size=len(indices), replace=True)
        if len(np.unique(y[sample])) < 2:
            continue
        boot.append(float(binary_auc(y[sample], scores[sample])))
    if not boot:
        return auc, np.nan, np.nan
    return auc, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def binary_auc(y: np.ndarray, scores: np.ndarray) -> float:
    """Rank-based ROC AUC with average ranks for ties."""
    y = y.astype(int)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        ranks[order[i:j]] = avg_rank
        i = j
    pos_rank_sum = ranks[y == 1].sum()
    return float((pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def build_peak_features(
    metadata: pd.DataFrame,
    median_peaks_path: Path,
    mz_min: float,
    mz_max: float,
    bin_width: float,
) -> tuple[pd.DataFrame, list[str]]:
    allowed = set(metadata["TGNR"].astype(str))
    peaks = pd.read_csv(median_peaks_path, usecols=["TGNR", "mass", "intensity"])
    peaks["TGNR"] = peaks["TGNR"].astype(str).str.replace("_", "-", regex=False)
    peaks = peaks[
        peaks["TGNR"].isin(allowed)
        & peaks["mass"].ge(mz_min)
        & peaks["mass"].lt(mz_max)
        & peaks["intensity"].gt(0)
    ].copy()
    peaks["bin"] = np.floor((peaks["mass"] - mz_min) / bin_width).astype(int)
    n_bins = int(math.ceil((mz_max - mz_min) / bin_width))

    wide = (
        peaks.groupby(["TGNR", "bin"], observed=True)["intensity"]
        .sum()
        .unstack(fill_value=0.0)
        .reindex(columns=range(n_bins), fill_value=0.0)
    )
    totals = wide.sum(axis=1)
    wide = wide.loc[totals > 0].copy()
    totals = totals.loc[wide.index]
    wide = np.log1p(wide.div(totals, axis=0) * 1_000_000.0)
    feature_names = [f"mz_{mz_min + i * bin_width:.0f}_{mz_min + (i + 1) * bin_width:.0f}" for i in range(n_bins)]
    wide.columns = feature_names
    return wide, feature_names


def stratified_folds(y: np.ndarray, folds: int, seed: int) -> list[np.ndarray]:
    counts = np.bincount(y.astype(int), minlength=2)
    n_splits = int(max(2, min(folds, counts.min())))
    rng = np.random.default_rng(seed)
    fold_indices: list[list[int]] = [[] for _ in range(n_splits)]
    for label in [0, 1]:
        indices = np.where(y == label)[0]
        rng.shuffle(indices)
        for idx, sample_idx in enumerate(indices):
            fold_indices[idx % n_splits].append(int(sample_idx))
    return [np.asarray(sorted(fold), dtype=int) for fold in fold_indices]


def centroid_scores(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> np.ndarray:
    pos_mean = x_train[y_train == 1].mean(axis=0)
    neg_mean = x_train[y_train == 0].mean(axis=0)
    direction = pos_mean - neg_mean
    # Centering by the training grand mean keeps the dot-product score stable
    # without adding a dependency-heavy classifier.
    center = x_train.mean(axis=0)
    return (x_test - center) @ direction


def oof_centroid_scores(
    x: np.ndarray,
    y: np.ndarray,
    folds: int,
    seed: int,
) -> tuple[np.ndarray, int]:
    counts = np.bincount(y.astype(int), minlength=2)
    n_splits = int(max(2, min(folds, counts.min())))
    if counts.min() < 2:
        raise ValueError("Need at least two samples per class for cross-validated AUC")
    folds_idx = stratified_folds(y, folds, seed)
    scores = np.full(len(y), np.nan, dtype=float)
    all_indices = np.arange(len(y))
    for test_idx in folds_idx:
        train_idx = np.setdiff1d(all_indices, test_idx, assume_unique=True)
        scores[test_idx] = centroid_scores(x[train_idx], y[train_idx], x[test_idx])
    return scores, n_splits


def group_prevalence_scores_oof(groups: np.ndarray, y: np.ndarray, folds: int, seed: int) -> tuple[np.ndarray, int]:
    counts = np.bincount(y.astype(int), minlength=2)
    n_splits = int(max(2, min(folds, counts.min())))
    folds_idx = stratified_folds(y, folds, seed)
    scores = np.full(len(y), np.nan, dtype=float)
    all_indices = np.arange(len(y))
    for test_idx in folds_idx:
        train_idx = np.setdiff1d(all_indices, test_idx, assume_unique=True)
        train = pd.DataFrame({"group": groups[train_idx], "y": y[train_idx]})
        global_rate = float(train["y"].mean())
        group_rate = train.groupby("group", observed=True)["y"].mean().to_dict()
        scores[test_idx] = [group_rate.get(group, global_rate) for group in groups[test_idx]]
    return scores, n_splits


def centered_auc_by_group(
    y: np.ndarray,
    scores: np.ndarray,
    groups: np.ndarray,
    min_pos: int,
    min_neg: int,
) -> dict[str, object]:
    frame = pd.DataFrame({"y": y.astype(int), "score": scores, "group": groups.astype(str)})
    support = frame.groupby("group", observed=True)["y"].agg(["count", "sum"])
    support["neg"] = support["count"] - support["sum"]
    valid_groups = support[(support["sum"] >= min_pos) & (support["neg"] >= min_neg)].index
    retained = frame[frame["group"].isin(valid_groups)].copy()
    if retained.empty or retained["y"].nunique() < 2:
        return {
            "status": "not_interpretable_no_valid_groups",
            "retained_n": 0,
            "retention": 0.0,
            "valid_groups": 0,
            "centered_auc": np.nan,
            "pairwise_within_group_accuracy": np.nan,
            "pairwise_comparisons": 0,
        }
    retained["centered_score"] = retained["score"] - retained.groupby("group", observed=True)["score"].transform("mean")
    centered_auc = float(binary_auc(retained["y"].to_numpy(), retained["centered_score"].to_numpy()))

    correct = 0
    total = 0
    for _, group_df in retained.groupby("group", observed=True):
        pos = group_df.loc[group_df["y"] == 1, "score"].to_numpy()
        neg = group_df.loc[group_df["y"] == 0, "score"].to_numpy()
        if len(pos) == 0 or len(neg) == 0:
            continue
        diffs = pos[:, None] - neg[None, :]
        correct += float((diffs > 0).sum() + 0.5 * (diffs == 0).sum())
        total += int(diffs.size)
    pair_acc = correct / total if total else np.nan
    return {
        "status": "ok",
        "retained_n": int(len(retained)),
        "retention": float(len(retained) / len(frame)),
        "valid_groups": int(len(valid_groups)),
        "centered_auc": centered_auc,
        "pairwise_within_group_accuracy": pair_acc,
        "pairwise_comparisons": int(total),
    }


def prepare_analysis_frame(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    metadata = pd.read_csv(args.metadata, sep="\t")
    metadata = metadata[metadata["has_bruker_maldi"].astype(str).eq("True")].copy()
    metadata["ST"] = metadata["ST"].astype(str)
    metadata["is_ST131_binary"] = metadata["is_ST131"].map(bool_to_binary)
    metadata["ciprofloxacin_R"] = metadata["vitek_CIP_result"].map(ast_to_binary)
    metadata["ceftriaxone_R"] = metadata["vitek_CTR_result"].map(ast_to_binary)
    metadata["phylogroup_rough"] = metadata["phylogroup_rough"].fillna("Unknown").astype(str)

    features, feature_names = build_peak_features(metadata, args.median_peaks, args.mz_min, args.mz_max, args.bin_width)
    metadata = metadata.set_index("TGNR").loc[features.index].reset_index()
    return metadata, features, feature_names


def subset_auc_rows(
    metadata: pd.DataFrame,
    features: pd.DataFrame,
    target: str,
    subset_col: str | None,
    subset_value: object | None,
    subset_name: str,
    args: argparse.Namespace,
) -> dict[str, object]:
    mask = metadata[target].notna()
    if subset_col is not None:
        mask &= metadata[subset_col].eq(subset_value)
    sub = metadata.loc[mask].copy()
    if sub.empty:
        return {"target": target, "subset": subset_name, "status": "skipped_no_rows", "n": 0}
    y = sub[target].astype(int).to_numpy()
    counts = np.bincount(y, minlength=2)
    row = {
        "target": target,
        "subset": subset_name,
        "status": "ok",
        "n": int(len(sub)),
        "class_0": int(counts[0]),
        "class_1": int(counts[1]),
    }
    if counts.min() < 2:
        row.update({"status": "skipped_insufficient_classes", "auc": np.nan, "auc_ci_low": np.nan, "auc_ci_high": np.nan})
        return row
    x = features.loc[sub["TGNR"]].to_numpy(dtype=float)
    scores, n_folds = oof_centroid_scores(x, y, args.folds, args.seed)
    auc, low, high = auc_ci(y, scores, args.bootstrap_n, args.seed)
    row.update(
        {
            "auc": auc,
            "auc_ci_low": low,
            "auc_ci_high": high,
            "folds": n_folds,
            "model": "centroid_direction",
        }
    )
    return row


def lineage_control_rows(
    metadata: pd.DataFrame,
    features: pd.DataFrame,
    target: str,
    background: str,
    args: argparse.Namespace,
) -> dict[str, object]:
    sub = metadata[metadata[target].notna() & metadata[background].notna()].copy()
    y = sub[target].astype(int).to_numpy()
    groups = sub[background].astype(str).to_numpy()
    counts = np.bincount(y, minlength=2)
    row = {
        "target": target,
        "background_control": background,
        "n": int(len(sub)),
        "class_0": int(counts[0]),
        "class_1": int(counts[1]),
    }
    if counts.min() < 2:
        row.update({"status": "skipped_insufficient_classes"})
        return row
    x = features.loc[sub["TGNR"]].to_numpy(dtype=float)
    scores, n_folds = oof_centroid_scores(x, y, args.folds, args.seed)
    raw_auc, raw_low, raw_high = auc_ci(y, scores, args.bootstrap_n, args.seed)
    bg_scores, _ = group_prevalence_scores_oof(groups, y, args.folds, args.seed)
    bg_auc, bg_low, bg_high = auc_ci(y, bg_scores, args.bootstrap_n, args.seed)
    centered = centered_auc_by_group(y, scores, groups, args.min_pos_per_group, args.min_neg_per_group)
    row.update(
        {
            "status": centered["status"],
            "raw_auc": raw_auc,
            "raw_auc_ci_low": raw_low,
            "raw_auc_ci_high": raw_high,
            "background_only_auc": bg_auc,
            "background_only_auc_ci_low": bg_low,
            "background_only_auc_ci_high": bg_high,
            "centered_auc": centered["centered_auc"],
            "centered_auc_drop": raw_auc - centered["centered_auc"] if not np.isnan(centered["centered_auc"]) else np.nan,
            "retained_n": centered["retained_n"],
            "retention": centered["retention"],
            "valid_groups": centered["valid_groups"],
            "pairwise_within_group_accuracy": centered["pairwise_within_group_accuracy"],
            "pairwise_comparisons": centered["pairwise_comparisons"],
            "folds": n_folds,
        }
    )
    return row


def st_support_table(metadata: pd.DataFrame, target: str, min_total: int) -> pd.DataFrame:
    frame = metadata[metadata[target].notna()].copy()
    frame[target] = frame[target].astype(int)
    rows = []
    for st, group in frame.groupby("ST", observed=True):
        counts = group[target].value_counts().to_dict()
        n0 = int(counts.get(0, 0))
        n1 = int(counts.get(1, 0))
        if n0 + n1 < min_total:
            continue
        rows.append(
            {
                "ST": st,
                "n": n0 + n1,
                "susceptible": n0,
                "resistant": n1,
                "has_both_classes": bool(n0 > 0 and n1 > 0),
                "resistant_fraction": n1 / (n0 + n1),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "ST"], ascending=[False, True])


def write_markdown_report(
    output_dir: Path,
    metadata: pd.DataFrame,
    subset_df: pd.DataFrame,
    control_df: pd.DataFrame,
) -> None:
    lines = [
        "# UPEC Clone-Controlled MALDI/WGS Analysis",
        "",
        f"Bruker MALDI-linked isolates with nonzero peak features: {len(metadata)}",
        f"ST131 isolates: {int(metadata['is_ST131_binary'].sum())}",
        "",
        "## Subset AUCs",
        "",
    ]
    for _, row in subset_df.iterrows():
        auc = row.get("auc")
        auc_text = "" if pd.isna(auc) else f"{auc:.3f}"
        lines.append(
            f"- {row['target']} / {row['subset']}: n={int(row['n'])}, "
            f"class_1={int(row.get('class_1', 0))}, AUC={auc_text}, status={row['status']}"
        )
    lines += ["", "## Lineage-Controlled AUCs", ""]
    for _, row in control_df.iterrows():
        centered = row.get("centered_auc")
        centered_text = "" if pd.isna(centered) else f"{centered:.3f}"
        lines.append(
            f"- {row['target']} controlled by {row['background_control']}: "
            f"raw={row.get('raw_auc', np.nan):.3f}, background_only={row.get('background_only_auc', np.nan):.3f}, "
            f"centered={centered_text}, retention={row.get('retention', np.nan):.2f}, "
            f"valid_groups={row.get('valid_groups', '')}"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "- A high ST131 AUC means MALDI peaks encode lineage strongly.",
        "- If ciprofloxacin AUC drops after ST/ST131/phylogroup centering, part of the resistance signal is lineage/background-associated.",
        "- If ciprofloxacin AUC remains above chance within ST131 or outside ST131, that is evidence for residual focal resistance-associated signal beyond the dominant clone.",
        "- Exact-ST centering is the strictest control and may retain fewer isolates because many STs contain only susceptible or only resistant isolates.",
        "",
    ]
    (output_dir / "clone_control_summary.md").write_text("\n".join(lines) + "\n")


def make_figure(output_dir: Path, subset_df: pd.DataFrame, control_df: pd.DataFrame) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    plot_subset = subset_df[subset_df["status"].eq("ok")].copy()
    labels = [f"{row.target}\n{row.subset}" for row in plot_subset.itertuples()]
    axes[0].barh(range(len(plot_subset)), plot_subset["auc"], color="#4C78A8")
    axes[0].axvline(0.5, color="0.55", linestyle="--", linewidth=1)
    axes[0].set_yticks(range(len(plot_subset)))
    axes[0].set_yticklabels(labels, fontsize=8)
    axes[0].set_xlim(0.3, 1.0)
    axes[0].set_xlabel("Cross-validated AUC")
    axes[0].set_title("Subset prediction")

    cipro = control_df[control_df["target"].eq("ciprofloxacin_R")].copy()
    x = np.arange(len(cipro))
    width = 0.25
    axes[1].bar(x - width, cipro["raw_auc"], width, label="MALDI raw", color="#4C78A8")
    axes[1].bar(x, cipro["background_only_auc"], width, label="background only", color="#F58518")
    axes[1].bar(x + width, cipro["centered_auc"], width, label="background-centered", color="#54A24B")
    axes[1].axhline(0.5, color="0.55", linestyle="--", linewidth=1)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cipro["background_control"], rotation=25, ha="right")
    axes[1].set_ylim(0.3, 1.0)
    axes[1].set_ylabel("AUC")
    axes[1].set_title("Ciprofloxacin after lineage controls")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(output_dir / "clone_control_auc_summary.png", dpi=300)
    plt.close(fig)


def make_svg_figure(output_dir: Path, subset_df: pd.DataFrame, control_df: pd.DataFrame) -> None:
    """Write a dependency-free SVG summary figure."""
    subset = subset_df[subset_df["status"].eq("ok")].copy()
    cipro = control_df[control_df["target"].eq("ciprofloxacin_R")].copy()
    width, height = 1100, 520
    left_x, right_x = 80, 630
    top = 70
    bar_h = 24
    gap = 14
    axis_w = 380
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#222} .small{font-size:12px}.label{font-size:11px}.title{font-size:17px;font-weight:700}.axis{stroke:#777;stroke-width:1}.tick{stroke:#ddd;stroke-width:1}</style>',
        f'<text class="title" x="{left_x}" y="35">A. MALDI peak prediction within clone strata</text>',
        f'<text class="title" x="{right_x}" y="35">B. Ciprofloxacin after WGS background controls</text>',
    ]
    def xscale(val: float, origin: int) -> float:
        return origin + (val - 0.3) / 0.7 * axis_w
    for origin in [left_x, right_x]:
        for tick in [0.3, 0.5, 0.7, 0.9, 1.0]:
            x = xscale(tick, origin)
            svg.append(f'<line class="tick" x1="{x:.1f}" y1="55" x2="{x:.1f}" y2="{height-55}"/>')
            svg.append(f'<text class="small" x="{x-10:.1f}" y="{height-25}">{tick:.1f}</text>')
        svg.append(f'<line class="axis" x1="{xscale(0.5, origin):.1f}" y1="55" x2="{xscale(0.5, origin):.1f}" y2="{height-55}" stroke-dasharray="4 4"/>')

    for i, row in enumerate(subset.itertuples(index=False)):
        y = top + i * (bar_h + gap)
        auc = float(row.auc)
        x0 = xscale(0.3, left_x)
        x1 = xscale(auc, left_x)
        label = f"{str(row.target).replace('_', ' ')} / {row.subset}"
        svg.append(f'<text class="label" x="{left_x}" y="{y-4}">{label}</text>')
        svg.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(1, x1-x0):.1f}" height="{bar_h}" fill="#4C78A8"/>')
        svg.append(f'<text class="small" x="{x1+5:.1f}" y="{y+17:.1f}">{auc:.3f}</text>')

    colors = [("#4C78A8", "raw"), ("#F58518", "lineage only"), ("#54A24B", "lineage-centered")]
    grouped = [
        ("raw_auc", -bar_h - 2, colors[0]),
        ("background_only_auc", 0, colors[1]),
        ("centered_auc", bar_h + 2, colors[2]),
    ]
    group_gap = 88
    for i, row in enumerate(cipro.itertuples(index=False)):
        base_y = top + 20 + i * group_gap
        svg.append(f'<text class="label" x="{right_x}" y="{base_y-30}">{row.background_control}</text>')
        for col, offset, (color, _) in grouped:
            val = float(getattr(row, col))
            x0 = xscale(0.3, right_x)
            x1 = xscale(val, right_x)
            y = base_y + offset
            svg.append(f'<rect x="{x0:.1f}" y="{y:.1f}" width="{max(1, x1-x0):.1f}" height="{bar_h}" fill="{color}"/>')
            svg.append(f'<text class="small" x="{x1+5:.1f}" y="{y+17:.1f}">{val:.3f}</text>')
    leg_x, leg_y = right_x, height - 75
    for idx, (color, text) in enumerate(colors):
        y = leg_y + idx * 18
        svg.append(f'<rect x="{leg_x}" y="{y-10}" width="12" height="12" fill="{color}"/>')
        svg.append(f'<text class="small" x="{leg_x+18}" y="{y}">{text}</text>')
    svg.append(f'<text class="small" x="{left_x + axis_w/2 - 20}" y="{height-5}">AUC</text>')
    svg.append(f'<text class="small" x="{right_x + axis_w/2 - 20}" y="{height-5}">AUC</text>')
    svg.append("</svg>")
    (output_dir / "clone_control_auc_summary.svg").write_text("\n".join(svg) + "\n")


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata, features, feature_names = prepare_analysis_frame(args)

    subset_specs = [
        ("is_ST131_binary", None, None, "overall"),
        ("ciprofloxacin_R", None, None, "overall"),
        ("ciprofloxacin_R", "is_ST131_binary", 1.0, "ST131_only"),
        ("ciprofloxacin_R", "is_ST131_binary", 0.0, "non_ST131_only"),
        ("ceftriaxone_R", None, None, "overall"),
        ("ceftriaxone_R", "is_ST131_binary", 1.0, "ST131_only"),
        ("ceftriaxone_R", "is_ST131_binary", 0.0, "non_ST131_only"),
    ]
    subset_rows = [
        subset_auc_rows(metadata, features, target, subset_col, subset_value, name, args)
        for target, subset_col, subset_value, name in subset_specs
    ]
    subset_df = pd.DataFrame(subset_rows)
    subset_df.to_csv(args.output_dir / "clone_control_subset_auc.csv", index=False)

    control_rows = []
    for target in ["ciprofloxacin_R", "ceftriaxone_R"]:
        for background in ["is_ST131_binary", "ST", "phylogroup_rough"]:
            control_rows.append(lineage_control_rows(metadata, features, target, background, args))
    control_df = pd.DataFrame(control_rows)
    control_df.to_csv(args.output_dir / "lineage_controlled_auc.csv", index=False)

    st_support = st_support_table(metadata, "ciprofloxacin_R", args.min_st_table_n)
    st_support.to_csv(args.output_dir / "ciprofloxacin_st_support.csv", index=False)

    manifest = {
        "metadata": str(args.metadata),
        "median_peaks": str(args.median_peaks),
        "n_isolates_with_features": int(len(metadata)),
        "n_features": int(len(feature_names)),
        "mz_min": args.mz_min,
        "mz_max": args.mz_max,
        "bin_width": args.bin_width,
        "folds": args.folds,
        "bootstrap_n": args.bootstrap_n,
        "seed": args.seed,
        "min_pos_per_group": args.min_pos_per_group,
        "min_neg_per_group": args.min_neg_per_group,
    }
    (args.output_dir / "clone_control_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_markdown_report(args.output_dir, metadata, subset_df, control_df)
    make_svg_figure(args.output_dir, subset_df, control_df)
    make_figure(args.output_dir, subset_df, control_df)

    print(f"Wrote clone-control outputs to {args.output_dir}")
    print(subset_df.to_string(index=False))
    print()
    print(control_df.to_string(index=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--median-peaks", type=Path, default=DEFAULT_MEDIAN_PEAKS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--mz-min", type=float, default=2000.0)
    parser.add_argument("--mz-max", type=float, default=20000.0)
    parser.add_argument("--bin-width", type=float, default=50.0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--bootstrap-n", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--min-pos-per-group", type=int, default=1)
    parser.add_argument("--min-neg-per-group", type=int, default=1)
    parser.add_argument("--min-st-table-n", type=int, default=4)
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
