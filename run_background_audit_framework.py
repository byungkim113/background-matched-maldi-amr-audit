#!/usr/bin/env python3
"""Model-agnostic Background-Matched Transfer Audit.

This script takes prediction rows from any MALDI-AMR model and asks whether the
prediction for a focal antibiotic survives after matching isolates on their
co-resistance background. It is intentionally independent of Mega_Model, Torch,
LightGBM, pandas, and DRIAMS-specific paths.

Minimum input is a long CSV with one row per isolate/drug prediction:

    isolate_id,site,year,organism,drug,label,prob

Column names are configurable so outputs from Weis-style scripts or other
models can be adapted without changing code.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


DEFAULT_BOOTSTRAP_N = 500
DEFAULT_PERMUTATION_N = 500
DEFAULT_RANDOM_SEED = 20260427

LABEL_CHAR = {0: "S", 1: "R"}
UNKNOWN_CHAR = "U"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a model-agnostic background-matched MALDI-AMR audit."
    )
    parser.add_argument("--predictions-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--id-col", default="isolate_id")
    parser.add_argument("--site-col", default="site")
    parser.add_argument("--year-col", default="year")
    parser.add_argument("--organism-col", default="organism")
    parser.add_argument("--drug-col", default="drug")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--prob-col", default="prob")
    parser.add_argument(
        "--background-signature-col",
        default="",
        help=(
            "Optional input column containing a precomputed co-resistance "
            "background signature. If supplied, the audit uses this column "
            "instead of deriving signatures from other prediction rows."
        ),
    )
    parser.add_argument("--model-name", default="")
    parser.add_argument(
        "--background-drugs",
        default="",
        help="Optional comma-separated drug list to use in background signatures.",
    )
    parser.add_argument("--min-pos-per-stratum", type=int, default=3)
    parser.add_argument("--min-neg-per-stratum", type=int, default=3)
    parser.add_argument("--match-year", action="store_true")
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP_N)
    parser.add_argument("--permutation-n", type=int, default=DEFAULT_PERMUTATION_N)
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED)
    parser.add_argument("--adequacy-min-n-matched", type=int, default=100)
    parser.add_argument("--adequacy-min-retention", type=float, default=0.10)
    parser.add_argument("--adequacy-min-pairwise", type=int, default=100)
    parser.add_argument(
        "--sensitivity-thresholds",
        default="2,3,5",
        help="Comma-separated min pos/neg stratum thresholds for sensitivity.",
    )
    parser.add_argument("--min-edge-n", type=int, default=30)
    return parser.parse_args()


def normalize_label(value) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    upper = text.upper()
    if upper in {"R", "RESISTANT", "TRUE", "T", "YES", "Y"}:
        return 1
    if upper in {"S", "SUSCEPTIBLE", "FALSE", "F", "NO", "N"}:
        return 0
    if upper in {"I", "INTERMEDIATE", "U", "UNKNOWN", "NA", "NAN", "NONE", "-", "."}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if math.isnan(number):
        return None
    if number == 1:
        return 1
    if number == 0:
        return 0
    return None


def parse_probability(value) -> float | None:
    try:
        prob = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if not math.isfinite(prob):
        return None
    return prob


def require_columns(fieldnames: Sequence[str], required: Sequence[str], path: Path | None = None) -> None:
    missing = [name for name in required if name not in fieldnames]
    if missing:
        prefix = f"{path}: " if path else ""
        raise ValueError(f"{prefix}missing required columns: {', '.join(missing)}")


def read_rows_from_records(
    records: Iterable[dict],
    *,
    id_col: str,
    site_col: str,
    year_col: str,
    organism_col: str,
    drug_col: str,
    label_col: str,
    prob_col: str,
    background_signature_col: str = "",
    model_name: str = "",
) -> List[dict]:
    rows: List[dict] = []
    for record in records:
        label = normalize_label(record.get(label_col))
        prob = parse_probability(record.get(prob_col))
        if label is None or prob is None:
            continue
        row = {
            "uid": str(record.get(id_col, "")).strip(),
            "site": str(record.get(site_col, "")).strip(),
            "year": str(record.get(year_col, "")).strip(),
            "organism": str(record.get(organism_col, "")).strip(),
            "drug": str(record.get(drug_col, "")).strip(),
            "label": label,
            "prob": prob,
            "model_name": str(model_name or record.get("model_name", "")).strip(),
        }
        if background_signature_col:
            row["background_signature"] = str(record.get(background_signature_col, "")).strip()
        rows.append(row)
    if not rows:
        raise ValueError("No usable prediction rows after label/probability normalization.")
    return rows


def read_long_predictions(
    path: Path | str,
    *,
    id_col: str,
    site_col: str,
    year_col: str,
    organism_col: str,
    drug_col: str,
    label_col: str,
    prob_col: str,
    background_signature_col: str = "",
    model_name: str = "",
) -> List[dict]:
    path = Path(path)
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        require_columns(
            reader.fieldnames or [],
            [
                id_col,
                site_col,
                year_col,
                organism_col,
                drug_col,
                label_col,
                prob_col,
                *([background_signature_col] if background_signature_col else []),
            ],
            path,
        )
        return read_rows_from_records(
            reader,
            id_col=id_col,
            site_col=site_col,
            year_col=year_col,
            organism_col=organism_col,
            drug_col=drug_col,
            label_col=label_col,
            prob_col=prob_col,
            background_signature_col=background_signature_col,
            model_name=model_name,
        )


def first_seen_drugs(rows: Sequence[dict]) -> List[str]:
    seen = set()
    drugs = []
    for row in rows:
        drug = row["drug"]
        if drug not in seen:
            seen.add(drug)
            drugs.append(drug)
    return drugs


def isolate_key(row: dict) -> Tuple[str, str, str, str]:
    return (row["site"], row["year"], row["organism"], row["uid"])


def label_matrix(rows: Sequence[dict]) -> Dict[Tuple[str, str, str, str], Dict[str, int]]:
    matrix: Dict[Tuple[str, str, str, str], Dict[str, int]] = defaultdict(dict)
    for row in rows:
        matrix[isolate_key(row)][row["drug"]] = int(row["label"])
    return matrix


def add_background_signatures(rows: Sequence[dict], background_drugs: Sequence[str] | None = None) -> List[dict]:
    drugs = list(background_drugs) if background_drugs else first_seen_drugs(rows)
    matrix = label_matrix(rows)
    enriched = []
    for row in rows:
        if row.get("background_signature"):
            new_row = dict(row)
            parts = [part for part in str(new_row["background_signature"]).split("|") if part]
            known = 0
            resistant = 0
            for part in parts:
                if "=" not in part:
                    continue
                _, char = part.rsplit("=", 1)
                char = char.strip().upper()
                if char in {"S", "R"}:
                    known += 1
                    resistant += int(char == "R")
            new_row["background_known_count"] = known
            new_row["background_resistant_count"] = resistant
            enriched.append(new_row)
            continue
        labels = matrix[isolate_key(row)]
        parts = []
        known = 0
        resistant = 0
        for drug in drugs:
            if drug == row["drug"]:
                continue
            label = labels.get(drug)
            char = LABEL_CHAR.get(label, UNKNOWN_CHAR)
            parts.append(f"{drug}={char}")
            if label in (0, 1):
                known += 1
                resistant += int(label == 1)
        new_row = dict(row)
        new_row["background_signature"] = "|".join(parts) if parts else "NO_BACKGROUND_DRUGS"
        new_row["background_known_count"] = known
        new_row["background_resistant_count"] = resistant
        enriched.append(new_row)
    return enriched


def safe_auc(labels: Sequence[int], scores: Sequence[float]) -> float:
    pairs = [(int(y), float(s)) for y, s in zip(labels, scores) if math.isfinite(float(s))]
    if not pairs:
        return math.nan
    labels = [p[0] for p in pairs]
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return math.nan

    order = sorted(range(len(pairs)), key=lambda i: pairs[i][1])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pairs[order[j + 1]][1] == pairs[order[i]][1]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        i = j + 1
    pos_rank_sum = sum(rank for rank, label in zip(ranks, labels) if label == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def safe_aupr(labels: Sequence[int], scores: Sequence[float]) -> float:
    pairs = [(int(y), float(s)) for y, s in zip(labels, scores) if math.isfinite(float(s))]
    n_pos = sum(y for y, _ in pairs)
    if n_pos == 0:
        return math.nan
    pairs.sort(key=lambda p: p[1], reverse=True)
    tp = 0
    precisions = []
    for idx, (label, _) in enumerate(pairs, start=1):
        if label == 1:
            tp += 1
            precisions.append(tp / idx)
    return sum(precisions) / n_pos


def bootstrap_metric_ci(
    labels: Sequence[int],
    scores: Sequence[float],
    metric_fn,
    *,
    n_boot: int,
    seed: int,
    ci: float = 0.95,
) -> Tuple[float, float]:
    if n_boot <= 0:
        return math.nan, math.nan
    pairs = [(int(y), float(s)) for y, s in zip(labels, scores) if math.isfinite(float(s))]
    if not pairs or sum(y for y, _ in pairs) == 0:
        return math.nan, math.nan
    rng = random.Random(seed)
    values = []
    for _ in range(n_boot):
        sample = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        value = metric_fn([y for y, _ in sample], [s for _, s in sample])
        if not math.isnan(value):
            values.append(value)
    if not values:
        return math.nan, math.nan
    values.sort()
    alpha = (1.0 - ci) / 2.0
    low_idx = min(len(values) - 1, max(0, int(math.floor(alpha * len(values)))))
    high_idx = min(len(values) - 1, max(0, int(math.ceil((1.0 - alpha) * len(values))) - 1))
    return values[low_idx], values[high_idx]


def group_by(rows: Iterable[dict], keys: Sequence[str]) -> Dict[Tuple, List[dict]]:
    grouped: Dict[Tuple, List[dict]] = defaultdict(list)
    for row in rows:
        grouped[tuple(row.get(key, "") for key in keys)].append(row)
    return grouped


def pairwise_accuracy_for_strata(rows: Sequence[dict], stratum_cols: Sequence[str], score_col: str = "prob") -> Tuple[float, int]:
    wins = 0.0
    total = 0
    for group in group_by(rows, stratum_cols).values():
        pos = [float(r[score_col]) for r in group if int(r["label"]) == 1]
        neg = [float(r[score_col]) for r in group if int(r["label"]) == 0]
        for p in pos:
            for n in neg:
                if p > n:
                    wins += 1.0
                elif p == n:
                    wins += 0.5
                total += 1
    return (wins / total if total else math.nan), total


def permutation_null_within_strata(
    rows: Sequence[dict],
    stratum_cols: Sequence[str],
    *,
    score_col: str,
    observed: float,
    n_perm: int,
    seed: int,
) -> Tuple[float, float, float]:
    if n_perm <= 0 or not rows or math.isnan(observed):
        return math.nan, math.nan, math.nan
    rng = random.Random(seed)
    groups = list(group_by(rows, stratum_cols).values())
    null_values = []
    for _ in range(n_perm):
        shuffled_rows = []
        for group in groups:
            labels = [int(r["label"]) for r in group]
            rng.shuffle(labels)
            for row, label in zip(group, labels):
                new_row = dict(row)
                new_row["label"] = label
                shuffled_rows.append(new_row)
        value = safe_auc(
            [int(r["label"]) for r in shuffled_rows],
            [float(r[score_col]) for r in shuffled_rows],
        )
        if not math.isnan(value):
            null_values.append(value)
    if not null_values:
        return math.nan, math.nan, math.nan
    p_value = (1.0 + sum(value >= observed for value in null_values)) / (1.0 + len(null_values))
    mean = sum(null_values) / len(null_values)
    if len(null_values) > 1:
        var = sum((value - mean) ** 2 for value in null_values) / (len(null_values) - 1)
        std = math.sqrt(var)
    else:
        std = 0.0
    return p_value, mean, std


def assign_adequacy_label(
    *,
    n_matched: int,
    matched_retention: float,
    n_valid_strata: int,
    pairwise_comparisons: int,
    min_n_matched: int,
    min_retention: float,
    min_pairwise: int,
) -> str:
    if n_matched == 0 or n_valid_strata == 0:
        return "not_interpretable_no_valid_strata"
    issues = []
    if n_matched < min_n_matched:
        issues.append("low_n_matched")
    if matched_retention < min_retention:
        issues.append("low_retention")
    if pairwise_comparisons < min_pairwise:
        issues.append("low_pairwise")
    if issues:
        return "caution_" + "_and_".join(issues)
    return "interpretable"


def interpretation_category(row: dict) -> str:
    adequacy = str(row.get("adequacy_label", ""))
    raw_auc = float(row.get("raw_auc", math.nan))
    centered = float(row.get("stratum_centered_auc", math.nan))
    if adequacy.startswith("not_interpretable"):
        return "insufficient_matched_overlap"
    if adequacy.startswith("caution"):
        return "caution_low_matched_support"
    if math.isnan(raw_auc) or raw_auc < 0.60:
        return "weak_raw_signal"
    if math.isnan(centered):
        return "insufficient_matched_overlap"
    if centered >= 0.60:
        return "focal_signal_retained"
    if centered >= 0.55:
        return "partially_retained_or_uncertain"
    return "background_driven_collapse"


def compute_background_matched_summary(
    rows: Sequence[dict],
    *,
    min_pos_per_stratum: int = 3,
    min_neg_per_stratum: int = 3,
    match_year: bool = False,
    bootstrap_n: int = DEFAULT_BOOTSTRAP_N,
    permutation_n: int = DEFAULT_PERMUTATION_N,
    random_seed: int = DEFAULT_RANDOM_SEED,
    adequacy_min_n_matched: int = 100,
    adequacy_min_retention: float = 0.10,
    adequacy_min_pairwise: int = 100,
) -> Tuple[List[dict], List[dict]]:
    summary = []
    retained_rows = []
    stratum_cols = ["background_signature"] + (["year"] if match_year else [])
    pair_groups = group_by(rows, ["site", "organism", "drug"])

    for group_idx, ((site, organism, drug), pair_rows) in enumerate(sorted(pair_groups.items())):
        labels = [int(r["label"]) for r in pair_rows]
        probs = [float(r["prob"]) for r in pair_rows]
        raw_auc = safe_auc(labels, probs)
        raw_aupr = safe_aupr(labels, probs)
        seed_base = int(random_seed) + group_idx * 100003
        raw_low, raw_high = bootstrap_metric_ci(labels, probs, safe_auc, n_boot=bootstrap_n, seed=seed_base + 1)

        matched = []
        for stratum_rows in group_by(pair_rows, stratum_cols).values():
            n_pos = sum(int(r["label"]) == 1 for r in stratum_rows)
            n_neg = sum(int(r["label"]) == 0 for r in stratum_rows)
            if n_pos >= min_pos_per_stratum and n_neg >= min_neg_per_stratum:
                matched.extend(stratum_rows)

        if matched:
            centered = []
            for stratum_rows in group_by(matched, stratum_cols).values():
                mean_prob = sum(float(r["prob"]) for r in stratum_rows) / len(stratum_rows)
                for row in stratum_rows:
                    new_row = dict(row)
                    new_row["centered_prob"] = float(row["prob"]) - mean_prob
                    new_row["matched_valid_stratum"] = True
                    centered.append(new_row)
            matched = centered
            matched_labels = [int(r["label"]) for r in matched]
            matched_probs = [float(r["prob"]) for r in matched]
            centered_probs = [float(r["centered_prob"]) for r in matched]
            matched_auc = safe_auc(matched_labels, matched_probs)
            matched_aupr = safe_aupr(matched_labels, matched_probs)
            centered_auc = safe_auc(matched_labels, centered_probs)
            matched_low, matched_high = bootstrap_metric_ci(
                matched_labels, matched_probs, safe_auc, n_boot=bootstrap_n, seed=seed_base + 2
            )
            centered_low, centered_high = bootstrap_metric_ci(
                matched_labels, centered_probs, safe_auc, n_boot=bootstrap_n, seed=seed_base + 3
            )
            permutation_p, null_mean, null_std = permutation_null_within_strata(
                matched,
                stratum_cols,
                score_col="centered_prob",
                observed=centered_auc,
                n_perm=permutation_n,
                seed=seed_base + 4,
            )
            pair_acc, pair_n = pairwise_accuracy_for_strata(matched, stratum_cols)
            matched_retention = len(matched) / len(pair_rows)
            n_valid_strata = len(group_by(matched, stratum_cols))
            retained_rows.extend(matched)
        else:
            matched_auc = matched_aupr = centered_auc = math.nan
            matched_low = matched_high = centered_low = centered_high = math.nan
            permutation_p = null_mean = null_std = math.nan
            pair_acc = math.nan
            pair_n = 0
            matched_retention = 0.0
            n_valid_strata = 0

        adequacy = assign_adequacy_label(
            n_matched=len(matched),
            matched_retention=matched_retention,
            n_valid_strata=n_valid_strata,
            pairwise_comparisons=pair_n,
            min_n_matched=adequacy_min_n_matched,
            min_retention=adequacy_min_retention,
            min_pairwise=adequacy_min_pairwise,
        )
        row = {
            "site": site,
            "organism": organism,
            "drug": drug,
            "raw_auc": raw_auc,
            "raw_auc_ci_low": raw_low,
            "raw_auc_ci_high": raw_high,
            "raw_aupr": raw_aupr,
            "matched_auc": matched_auc,
            "matched_auc_ci_low": matched_low,
            "matched_auc_ci_high": matched_high,
            "matched_aupr": matched_aupr,
            "stratum_centered_auc": centered_auc,
            "stratum_centered_auc_ci_low": centered_low,
            "stratum_centered_auc_ci_high": centered_high,
            "pairwise_accuracy": pair_acc,
            "pairwise_comparisons": pair_n,
            "permutation_p": permutation_p,
            "permutation_null_mean": null_mean,
            "permutation_null_std": null_std,
            "matched_retention": matched_retention,
            "adequacy_label": adequacy,
            "interpretation_category": "",
            "n_total": len(pair_rows),
            "n_r": sum(labels),
            "n_matched": len(matched),
            "n_matched_r": sum(int(r["label"]) for r in matched) if matched else 0,
            "n_valid_strata": n_valid_strata,
            "min_pos_per_stratum": min_pos_per_stratum,
            "min_neg_per_stratum": min_neg_per_stratum,
            "match_year": bool(match_year),
        }
        row["interpretation_category"] = interpretation_category(row)
        summary.append(row)
    return summary, retained_rows


def parse_thresholds(text: str) -> List[int]:
    values = []
    for part in str(text or "").split(","):
        part = part.strip()
        if not part:
            continue
        value = int(part)
        if value <= 0:
            raise ValueError("Sensitivity thresholds must be positive.")
        values.append(value)
    return sorted(set(values))


def phi_from_counts(n00: int, n01: int, n10: int, n11: int) -> float:
    denom = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    if denom == 0:
        return math.nan
    return (n11 * n00 - n10 * n01) / denom


def build_isolate_label_records(rows: Sequence[dict]) -> List[dict]:
    by_iso: Dict[Tuple[str, str, str, str], dict] = {}
    for row in rows:
        key = isolate_key(row)
        record = by_iso.setdefault(
            key,
            {"site": row["site"], "year": row["year"], "organism": row["organism"], "uid": row["uid"]},
        )
        record[row["drug"]] = int(row["label"])
    return list(by_iso.values())


def pair_metrics(records: Sequence[dict], drug_a: str, drug_b: str) -> dict | None:
    n00 = n01 = n10 = n11 = 0
    for record in records:
        if drug_a not in record or drug_b not in record:
            continue
        a = int(record[drug_a])
        b = int(record[drug_b])
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
        return None
    prev_a = (n10 + n11) / n
    prev_b = (n01 + n11) / n
    expected_rr = prev_a * prev_b
    observed_rr = n11 / n
    lift = observed_rr / expected_rr if expected_rr > 0 else math.nan
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
        "phi": phi_from_counts(n00, n01, n10, n11),
        "resistant_jaccard": n11 / (n11 + n10 + n01) if (n11 + n10 + n01) else math.nan,
    }


def build_cross_resistance_edges(rows: Sequence[dict], min_edge_n: int) -> Tuple[List[dict], List[dict]]:
    records = build_isolate_label_records(rows)
    drugs = first_seen_drugs(rows)
    by_site: Dict[str, List[dict]] = defaultdict(list)
    by_site["ALL"] = list(records)
    for record in records:
        by_site[record["site"]].append(record)

    edges = []
    prevalence = []
    for site, site_records in sorted(by_site.items()):
        for drug in drugs:
            vals = [record[drug] for record in site_records if drug in record]
            if vals:
                prevalence.append(
                    {
                        "site": site,
                        "drug": drug,
                        "n_known": len(vals),
                        "n_resistant": sum(vals),
                        "prevalence": sum(vals) / len(vals),
                    }
                )
        for i, drug_a in enumerate(drugs):
            for drug_b in drugs[i + 1 :]:
                metrics = pair_metrics(site_records, drug_a, drug_b)
                if metrics and metrics["n_both_known"] >= min_edge_n:
                    edges.append({"site": site, **metrics})
    return edges, prevalence


def strongest_partner_by_drug(edges: Sequence[dict]) -> Dict[str, dict]:
    partners = {}
    for edge in edges:
        if edge.get("site") != "ALL":
            continue
        for drug, other in [(edge["drug_a"], edge["drug_b"]), (edge["drug_b"], edge["drug_a"])]:
            current = partners.get(drug)
            if current is None or float(edge["phi"]) > float(current["phi"]):
                partners[drug] = {"partner": other, "phi": edge["phi"], "rr_lift": edge["rr_lift"]}
    return partners


def add_ecology_annotations(summary: Sequence[dict], edges: Sequence[dict]) -> List[dict]:
    partners = strongest_partner_by_drug(edges)
    annotated = []
    for row in summary:
        partner = partners.get(row["drug"], {})
        new_row = dict(row)
        new_row["strongest_network_partner"] = partner.get("partner", "")
        new_row["partner_phi"] = partner.get("phi", math.nan)
        new_row["partner_rr_lift"] = partner.get("rr_lift", math.nan)
        annotated.append(new_row)
    return annotated


def safe_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.6g}"
    return str(value)


def write_csv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str] | None = None) -> None:
    if fieldnames is None:
        fields = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    fields.append(key)
    else:
        fields = list(fieldnames)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: safe_text(row.get(field)) for field in fields})


def write_markdown_table(path: Path, title: str, rows: Sequence[dict], fields: Sequence[str]) -> None:
    lines = [f"# {title}", ""]
    if not rows:
        lines.append("No rows.")
    else:
        lines.append("| " + " | ".join(fields) + " |")
        lines.append("| " + " | ".join(["---"] * len(fields)) + " |")
        for row in rows:
            lines.append("| " + " | ".join(safe_text(row.get(field)) for field in fields) + " |")
    path.write_text("\n".join(lines) + "\n")


def write_report(path: Path, summary: Sequence[dict], edges: Sequence[dict]) -> None:
    top_edges = sorted(
        [edge for edge in edges if edge["site"] == "ALL"],
        key=lambda edge: float(edge["phi"]),
        reverse=True,
    )[:10]
    lines = [
        "# Background-Matched Transfer Audit Report",
        "",
        "## What This Tests",
        "",
        "The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.",
        "Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.",
        "",
        "## Strongest Co-Resistance Edges",
        "",
        "| Drug A | Drug B | Phi | Lift | n RR | n |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for edge in top_edges:
        lines.append(
            f"| {edge['drug_a']} | {edge['drug_b']} | {safe_text(edge['phi'])} | "
            f"{safe_text(edge['rr_lift'])} | {edge['n_rr']} | {edge['n_both_known']} |"
        )
    lines += [
        "",
        "## Audit Summary",
        "",
        "| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |",
        "| --- | --- | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in summary:
        lines.append(
            f"| {row['site']} | {row['organism']} | {row['drug']} | {safe_text(row['raw_auc'])} | "
            f"{safe_text(row['stratum_centered_auc'])} | {safe_text(row['matched_retention'])} | "
            f"{row['adequacy_label']} | {row['interpretation_category']} |"
        )
    path.write_text("\n".join(lines) + "\n")


def write_minimal_network_svg(path: Path, edges: Sequence[dict]) -> None:
    all_edges = [edge for edge in edges if edge["site"] == "ALL"]
    drugs = []
    for edge in all_edges:
        for drug in (edge["drug_a"], edge["drug_b"]):
            if drug not in drugs:
                drugs.append(drug)
    if not drugs:
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="400" height="120"><text x="20" y="40">No network edges</text></svg>')
        return
    width = 760
    height = 520
    cx = width / 2
    cy = height / 2 + 20
    radius = 170
    pos = {}
    for idx, drug in enumerate(drugs):
        angle = -math.pi / 2 + 2 * math.pi * idx / len(drugs)
        pos[drug] = (cx + radius * math.cos(angle), cy + radius * math.sin(angle))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="24" y="34" font-family="Arial" font-size="19" font-weight="700">Cross-resistance network</text>',
    ]
    for edge in all_edges:
        phi = float(edge["phi"])
        if abs(phi) < 0.15:
            continue
        x1, y1 = pos[edge["drug_a"]]
        x2, y2 = pos[edge["drug_b"]]
        color = "#3b5bdb" if phi >= 0 else "#d9480f"
        width_px = 1.0 + 7.0 * min(abs(phi), 1.0)
        parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width_px:.2f}" stroke-opacity="0.58"/>')
    for drug, (x, y) in pos.items():
        label = drug if len(drug) <= 18 else drug[:17] + "..."
        parts.append(f'<circle cx="{x}" cy="{y}" r="30" fill="#f8f9fa" stroke="#222"/>')
        parts.append(f'<text x="{x}" y="{y + 4}" text-anchor="middle" font-family="Arial" font-size="11">{label}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    background_drugs = [part.strip() for part in args.background_drugs.split(",") if part.strip()]

    rows = read_long_predictions(
        args.predictions_csv,
        id_col=args.id_col,
        site_col=args.site_col,
        year_col=args.year_col,
        organism_col=args.organism_col,
        drug_col=args.drug_col,
        label_col=args.label_col,
        prob_col=args.prob_col,
        background_signature_col=args.background_signature_col,
        model_name=args.model_name,
    )
    rows = add_background_signatures(rows, background_drugs or None)
    summary, retained = compute_background_matched_summary(
        rows,
        min_pos_per_stratum=args.min_pos_per_stratum,
        min_neg_per_stratum=args.min_neg_per_stratum,
        match_year=args.match_year,
        bootstrap_n=args.bootstrap_n,
        permutation_n=args.permutation_n,
        random_seed=args.random_seed,
        adequacy_min_n_matched=args.adequacy_min_n_matched,
        adequacy_min_retention=args.adequacy_min_retention,
        adequacy_min_pairwise=args.adequacy_min_pairwise,
    )
    edges, prevalence = build_cross_resistance_edges(rows, args.min_edge_n)
    annotated = add_ecology_annotations(summary, edges)

    write_csv(output_dir / "normalized_predictions.csv", rows)
    write_csv(output_dir / "background_matched_audit_summary.csv", summary)
    write_csv(output_dir / "background_matched_retained_rows.csv", retained)
    write_csv(output_dir / "cross_resistance_edges.csv", edges)
    write_csv(output_dir / "cross_resistance_prevalence.csv", prevalence)
    write_csv(output_dir / "background_audit_with_resistance_ecology.csv", annotated)
    write_markdown_table(
        output_dir / "background_matched_audit_summary.md",
        "Background-Matched Audit Summary",
        summary,
        [
            "site",
            "organism",
            "drug",
            "raw_auc",
            "matched_auc",
            "stratum_centered_auc",
            "matched_retention",
            "adequacy_label",
            "interpretation_category",
        ],
    )

    sensitivity_rows = []
    for threshold in parse_thresholds(args.sensitivity_thresholds):
        sens_summary, _ = compute_background_matched_summary(
            rows,
            min_pos_per_stratum=threshold,
            min_neg_per_stratum=threshold,
            match_year=args.match_year,
            bootstrap_n=0,
            permutation_n=0,
            random_seed=args.random_seed,
            adequacy_min_n_matched=args.adequacy_min_n_matched,
            adequacy_min_retention=args.adequacy_min_retention,
            adequacy_min_pairwise=args.adequacy_min_pairwise,
        )
        for row in sens_summary:
            sensitivity_rows.append({"sensitivity_threshold": threshold, **row})
    write_csv(output_dir / "background_matched_sensitivity.csv", sensitivity_rows)
    write_report(output_dir / "background_audit_report.md", summary, edges)
    write_minimal_network_svg(output_dir / "cross_resistance_network.svg", edges)

    print(f"Read {len(rows)} normalized prediction rows.")
    print(f"Wrote audit outputs to {output_dir}")
    print(f"Summary rows: {len(summary)}")
    print(f"Cross-resistance edges: {len(edges)}")


if __name__ == "__main__":
    main()
