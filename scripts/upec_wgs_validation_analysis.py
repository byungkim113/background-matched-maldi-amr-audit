#!/usr/bin/env python3
"""Run first-pass WGS-linked MALDI validation analyses for the Basel UPEC cohort.

The script joins the public UPEC master metadata table to Bruker median-peak
MALDI data, builds TGNR-level binned spectral features, and quantifies whether
MALDI peak patterns can recover WGS lineage/background labels such as
phylogroup and ST131.

It writes small CSV/JSON/Markdown outputs and does not write raw spectra.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "data_manifests" / "upec_master_metadata.tsv"
DEFAULT_MEDIAN_PEAKS = ROOT / "data_manifests" / "Bruker_csv_medianpeaks_df.csv"
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "analysis_outputs" / "upec_wgs_validation_outputs"
PUBLISHED_ST131_PEAKS = [
    {
        "published_mz": 3236,
        "protein": "HdeA multivalent ion",
        "annotation": "multivalent ion of HdeA m/z 9710",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 4176,
        "protein": "YjbJ multivalent ion",
        "annotation": "multivalent ion of YjbJ m/z 8351",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 4857,
        "protein": "HdeA multivalent ion",
        "annotation": "multivalent ion of HdeA m/z 9710",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 5381,
        "protein": "unidentified",
        "annotation": "ST131-specific peak not identified by LC-MS/MS",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 6827,
        "protein": "unidentified",
        "annotation": "ST131-specific peak not identified by LC-MS/MS",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 7655,
        "protein": "YahO",
        "annotation": "uncharacterized protein YahO",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 8351,
        "protein": "YjbJ",
        "annotation": "UPF0337 protein YjbJ",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 8448,
        "protein": "YnfD",
        "annotation": "uncharacterized protein YnfD",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 9710,
        "protein": "HdeA",
        "annotation": "acid stress chaperone HdeA",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
    {
        "published_mz": 11783,
        "protein": "cytochrome b562",
        "annotation": "soluble cytochrome b562",
        "source": "Scientific Reports 2019; DOI 10.1038/s41598-019-45051-z",
    },
]


def parse_bool(value: str) -> int | None:
    text = str(value).strip().lower()
    if text == "true":
        return 1
    if text == "false":
        return 0
    return None


def normalize_ast_call(value: str) -> int | None:
    """Map German AST calls to binary labels: susceptible=0, resistant=1."""
    text = str(value).strip().lower()
    if text in {"sensibel", "susceptible", "s"}:
        return 0
    if text in {"resistent", "resistant", "r"}:
        return 1
    return None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def target_rows(
    metadata_rows: list[dict[str, str]],
    column: str,
    parser: Callable[[str], int | None],
) -> dict[str, int]:
    labels: dict[str, int] = {}
    for row in metadata_rows:
        label = parser(row.get(column, ""))
        if label is not None:
            labels[row["TGNR"]] = label
    return labels


def categorical_target_rows(
    metadata_rows: list[dict[str, str]],
    column: str,
    min_count: int,
) -> dict[str, str]:
    values = [row.get(column, "") for row in metadata_rows]
    counts = Counter(value for value in values if value and value not in {"NA", "nan"})
    keep = {value for value, count in counts.items() if count >= min_count}
    return {
        row["TGNR"]: row[column]
        for row in metadata_rows
        if row.get(column, "") in keep
    }


def binary_auc(y_true: list[int], scores: list[float]) -> float:
    """Compute rank-based ROC AUC without external dependencies."""
    if len(y_true) != len(scores):
        raise ValueError("y_true and scores must have the same length")
    n_pos = sum(y_true)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")

    order = sorted(range(len(scores)), key=lambda idx: scores[idx])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j

    pos_rank_sum = sum(rank for rank, label in zip(ranks, y_true) if label == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def build_peak_feature_matrix(
    median_peaks: Path,
    allowed_tgnrs: set[str],
    mz_min: float,
    mz_max: float,
    bin_width: float,
) -> tuple[list[str], list[list[float]], list[str]]:
    n_bins = int(math.ceil((mz_max - mz_min) / bin_width))
    features: dict[str, list[float]] = {tgnr: [0.0] * n_bins for tgnr in allowed_tgnrs}

    with median_peaks.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"mass", "intensity", "TGNR"}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Median peak file is missing required columns: {missing}")
        for row in reader:
            tgnr = row["TGNR"].strip().replace("_", "-")
            if tgnr not in features:
                continue
            try:
                mass = float(row["mass"])
                intensity = float(row["intensity"])
            except ValueError:
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

    feature_names = [f"mz_{mz_min + i * bin_width:.0f}_{mz_min + (i + 1) * bin_width:.0f}" for i in range(n_bins)]
    return tgnrs, matrix, feature_names


def feature_center(feature_name: str) -> float:
    _, low, high = feature_name.split("_")
    return (float(low) + float(high)) / 2.0


def match_literature_st131_peak(mz: float, tolerance_da: float) -> dict[str, object]:
    best = min(PUBLISHED_ST131_PEAKS, key=lambda row: abs(float(row["published_mz"]) - mz))
    delta = abs(float(best["published_mz"]) - mz)
    if delta <= tolerance_da:
        return {
            "published_mz": best["published_mz"],
            "protein": best["protein"],
            "annotation": best["annotation"],
            "source": best["source"],
            "delta_da": round(delta, 4),
        }
    return {
        "published_mz": "",
        "protein": "",
        "annotation": "",
        "source": "",
        "delta_da": "",
    }


def peak_effect_rows(
    tgnrs: list[str],
    x_matrix: list[list[float]],
    feature_names: list[str],
    labels_by_tgnr: dict[str, int],
    top_n: int,
) -> list[dict[str, object]]:
    usable = [(row, labels_by_tgnr[tgnr]) for tgnr, row in zip(tgnrs, x_matrix) if tgnr in labels_by_tgnr]
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


def score_simple_centroid_binary(
    x_train: list[list[float]], y_train: list[int], x_test: list[list[float]]
) -> list[float]:
    n_features = len(x_train[0])
    pos = [row for row, label in zip(x_train, y_train) if label == 1]
    neg = [row for row, label in zip(x_train, y_train) if label == 0]
    pos_mean = [mean(row[i] for row in pos) for i in range(n_features)]
    neg_mean = [mean(row[i] for row in neg) for i in range(n_features)]
    direction = [pos_mean[i] - neg_mean[i] for i in range(n_features)]
    return [sum(row[i] * direction[i] for i in range(n_features)) for row in x_test]


def cross_validated_binary_auc(
    tgnrs: list[str],
    x_matrix: list[list[float]],
    labels_by_tgnr: dict[str, int],
    folds: int,
) -> dict[str, object]:
    usable = [(tgnr, row, labels_by_tgnr[tgnr]) for tgnr, row in zip(tgnrs, x_matrix) if tgnr in labels_by_tgnr]
    counts = Counter(label for _, _, label in usable)
    if len(counts) != 2 or min(counts.values()) < 2:
        return {"status": "skipped_insufficient_classes", "n": len(usable), "auc": ""}

    # Deterministic stratified round-robin folds.
    by_class: dict[int, list[tuple[str, list[float], int]]] = defaultdict(list)
    for item in usable:
        by_class[item[2]].append(item)
    n_folds = max(2, min(folds, min(len(items) for items in by_class.values())))
    fold_items: list[list[tuple[str, list[float], int]]] = [[] for _ in range(n_folds)]
    for label in sorted(by_class):
        for idx, item in enumerate(by_class[label]):
            fold_items[idx % n_folds].append(item)

    y_all: list[int] = []
    score_all: list[float] = []
    for fold_idx in range(n_folds):
        test = fold_items[fold_idx]
        train = [item for idx, fold in enumerate(fold_items) if idx != fold_idx for item in fold]
        x_train = [item[1] for item in train]
        y_train = [item[2] for item in train]
        x_test = [item[1] for item in test]
        y_test = [item[2] for item in test]
        scores = score_simple_centroid_binary(x_train, y_train, x_test)
        y_all.extend(y_test)
        score_all.extend(scores)

    return {
        "status": "ok",
        "n": len(usable),
        "class_0": counts[0],
        "class_1": counts[1],
        "auc": round(binary_auc(y_all, score_all), 4),
        "folds": n_folds,
        "model": "centroid_direction",
    }


def run_sklearn_models(
    tgnrs: list[str],
    x_matrix: list[list[float]],
    targets: dict[str, dict[str, int | str]],
    folds: int,
) -> list[dict[str, object]]:
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import accuracy_score, roc_auc_score
        from sklearn.model_selection import StratifiedKFold
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
    except Exception as exc:  # pragma: no cover - depends on local environment
        return [{"target": "sklearn_models", "status": f"skipped_missing_sklearn:{exc}"}]

    x_all = np.asarray(x_matrix, dtype=float)
    tgnr_to_index = {tgnr: idx for idx, tgnr in enumerate(tgnrs)}
    rows: list[dict[str, object]] = []
    for target_name, labels in targets.items():
        usable = [(tgnr_to_index[tgnr], label) for tgnr, label in labels.items() if tgnr in tgnr_to_index]
        y = np.asarray([label for _, label in usable])
        if len(set(y.tolist())) < 2:
            rows.append({"target": target_name, "status": "skipped_one_class", "n": len(y)})
            continue
        counts = Counter(y.tolist())
        n_splits = max(2, min(folds, min(counts.values())))
        if n_splits < 2:
            rows.append({"target": target_name, "status": "skipped_insufficient_classes", "n": len(y)})
            continue
        x = x_all[[idx for idx, _ in usable]]
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=17)
        aucs: list[float] = []
        accs: list[float] = []
        for train_idx, test_idx in cv.split(x, y):
            clf = make_pipeline(
                StandardScaler(),
                LogisticRegression(max_iter=2000, class_weight="balanced", solver="liblinear"),
            )
            clf.fit(x[train_idx], y[train_idx])
            pred = clf.predict(x[test_idx])
            accs.append(float(accuracy_score(y[test_idx], pred)))
            if len(counts) == 2:
                prob = clf.predict_proba(x[test_idx])[:, 1]
                aucs.append(float(roc_auc_score(y[test_idx], prob)))
        rows.append(
            {
                "target": target_name,
                "status": "ok",
                "n": len(y),
                "classes": json.dumps(dict(counts), sort_keys=True),
                "mean_auc": round(mean(aucs), 4) if aucs else "",
                "mean_accuracy": round(mean(accs), 4),
                "folds": n_splits,
                "model": "logistic_regression",
            }
        )
    return rows


def association_table(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for phenotype_col, label in [("vitek_CIP_result", "Ciprofloxacin"), ("vitek_CTR_result", "Ceftriaxone")]:
        counts = Counter()
        for row in rows:
            if row.get("has_bruker_maldi") != "True":
                continue
            ast = normalize_ast_call(row.get(phenotype_col, ""))
            st131 = parse_bool(row.get("is_ST131", ""))
            if ast is None or st131 is None:
                continue
            counts[(st131, ast)] += 1
        a = counts[(1, 1)]
        b = counts[(1, 0)]
        c = counts[(0, 1)]
        d = counts[(0, 0)]
        odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        output.append(
            {
                "phenotype": label,
                "st131_resistant": a,
                "st131_susceptible": b,
                "nonst131_resistant": c,
                "nonst131_susceptible": d,
                "haldane_odds_ratio": round(odds_ratio, 4),
            }
        )
    return output


def run_analysis(args: argparse.Namespace) -> None:
    rows = read_tsv(args.metadata)
    bruker_rows = [row for row in rows if row.get("has_bruker_maldi") == "True"]
    allowed_tgnrs = {row["TGNR"] for row in bruker_rows}
    tgnrs, x_matrix, feature_names = build_peak_feature_matrix(
        args.median_peaks,
        allowed_tgnrs,
        args.mz_min,
        args.mz_max,
        args.bin_width,
    )

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    targets_binary = {
        "ST131": target_rows(rows, "is_ST131", parse_bool),
        "Ciprofloxacin_R": target_rows(rows, "vitek_CIP_result", normalize_ast_call),
        "Ceftriaxone_R": target_rows(rows, "vitek_CTR_result", normalize_ast_call),
    }
    phylogroup = categorical_target_rows(bruker_rows, "phylogroup_rough", args.min_class_count)
    phylogroup_encoded = {tgnr: label for tgnr, label in phylogroup.items()}

    centroid_rows = []
    for target_name, labels in targets_binary.items():
        result = cross_validated_binary_auc(tgnrs, x_matrix, labels, args.folds)
        result["target"] = target_name
        centroid_rows.append(result)
    write_csv(
        output_dir / "centroid_binary_cv_results.csv",
        centroid_rows,
        ["target", "status", "n", "class_0", "class_1", "auc", "folds", "model"],
    )

    sklearn_targets: dict[str, dict[str, int | str]] = dict(targets_binary)
    sklearn_targets["phylogroup_rough"] = phylogroup_encoded
    sklearn_rows = run_sklearn_models(tgnrs, x_matrix, sklearn_targets, args.folds)
    write_csv(
        output_dir / "sklearn_cv_results.csv",
        sklearn_rows,
        ["target", "status", "n", "classes", "mean_auc", "mean_accuracy", "folds", "model"],
    )

    associations = association_table(rows)
    write_csv(
        output_dir / "st131_resistance_associations.csv",
        associations,
        [
            "phenotype",
            "st131_resistant",
            "st131_susceptible",
            "nonst131_resistant",
            "nonst131_susceptible",
            "haldane_odds_ratio",
        ],
    )

    peak_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    for target_name, labels in targets_binary.items():
        target_peaks = peak_effect_rows(tgnrs, x_matrix, feature_names, labels, args.top_peaks)
        for row in target_peaks:
            row = dict(row)
            row["target"] = target_name
            peak_rows.append(row)
            match = match_literature_st131_peak(float(row["mz_center"]), args.literature_tolerance_da)
            if match["published_mz"] != "":
                overlap_rows.append(
                    {
                        "target": target_name,
                        "feature": row["feature"],
                        "mz_center": row["mz_center"],
                        "effect_pos_minus_neg": row["effect_pos_minus_neg"],
                        "published_mz": match["published_mz"],
                        "delta_da": match["delta_da"],
                        "protein": match["protein"],
                        "annotation": match["annotation"],
                        "source": match["source"],
                    }
                )
    write_csv(
        output_dir / "top_discriminative_peak_bins.csv",
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
        output_dir / "published_st131_peak_overlap.csv",
        overlap_rows,
        [
            "target",
            "feature",
            "mz_center",
            "effect_pos_minus_neg",
            "published_mz",
            "delta_da",
            "protein",
            "annotation",
            "source",
        ],
    )

    summary = {
        "metadata_rows": len(rows),
        "bruker_maldi_metadata_rows": len(bruker_rows),
        "feature_tgnrs_with_nonzero_peaks": len(tgnrs),
        "feature_count": len(feature_names),
        "mz_min": args.mz_min,
        "mz_max": args.mz_max,
        "bin_width": args.bin_width,
        "top_phylogroup_rough_counts": Counter(row["phylogroup_rough"] for row in bruker_rows).most_common(20),
        "top_st_counts": Counter(row["ST"] for row in bruker_rows).most_common(20),
        "ciprofloxacin_result_counts": Counter(row["vitek_CIP_result"] for row in bruker_rows).most_common(),
        "ceftriaxone_result_counts": Counter(row["vitek_CTR_result"] for row in bruker_rows).most_common(),
        "published_st131_peak_overlap_rows": len(overlap_rows),
    }
    (output_dir / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    report_lines = [
        "# UPEC WGS-Linked MALDI Validation Summary",
        "",
        f"Metadata rows: {len(rows)}",
        f"Bruker MALDI-linked metadata rows: {len(bruker_rows)}",
        f"TGNRs with nonzero binned peaks: {len(tgnrs)}",
        f"Feature bins: {len(feature_names)}",
        "",
        "## Centroid Binary CV",
        "",
    ]
    for row in centroid_rows:
        report_lines.append(
            f"- {row['target']}: status={row['status']}, n={row.get('n', '')}, AUC={row.get('auc', '')}"
        )
    report_lines += ["", "## ST131 Resistance Association", ""]
    for row in associations:
        report_lines.append(
            f"- {row['phenotype']}: OR={row['haldane_odds_ratio']}, "
            f"ST131 R/S={row['st131_resistant']}/{row['st131_susceptible']}, "
            f"non-ST131 R/S={row['nonst131_resistant']}/{row['nonst131_susceptible']}"
        )
    report_lines += ["", "## Published ST131 Proteomic Peak Overlap", ""]
    if overlap_rows:
        for row in overlap_rows[:20]:
            report_lines.append(
                f"- {row['target']} peak {row['mz_center']} matched published m/z "
                f"{row['published_mz']} ({row['protein']}, delta={row['delta_da']} Da)"
            )
    else:
        report_lines.append("- No top peak bins matched published ST131 peaks within the selected tolerance.")
    (output_dir / "validation_summary.md").write_text("\n".join(report_lines) + "\n")

    print(f"Wrote outputs to {output_dir}")
    print(f"TGNRs with nonzero peaks: {len(tgnrs)}")
    for row in centroid_rows:
        print(f"{row['target']}: {row['status']} n={row.get('n', '')} auc={row.get('auc', '')}")
    print(f"Published ST131 peak overlaps: {len(overlap_rows)}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--median-peaks", type=Path, default=DEFAULT_MEDIAN_PEAKS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--mz-min", type=float, default=2000.0)
    parser.add_argument("--mz-max", type=float, default=20000.0)
    parser.add_argument("--bin-width", type=float, default=50.0)
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--min-class-count", type=int, default=10)
    parser.add_argument("--top-peaks", type=int, default=30)
    parser.add_argument("--literature-tolerance-da", type=float, default=75.0)
    return parser


def main() -> None:
    run_analysis(build_parser().parse_args())


if __name__ == "__main__":
    main()
