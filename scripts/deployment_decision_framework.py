#!/usr/bin/env python3
"""Convert audit and calibration results into deployment-facing decisions."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Sequence


DEFAULT_AUDIT = "outputs/final_framework_outputs/table_1_primary_background_matched_audit.csv"
DEFAULT_CALIBRATION = "outputs/analysis_outputs/calibration_analysis/calibration_summary.csv"
DEFAULT_OUTPUT_DIR = "outputs/analysis_outputs/deployment_decision_framework"

RULE_FIELDS = ["scenario", "decision_category", "interpretation", "recommended_action"]
READINESS_FIELDS = [
    "site",
    "pair",
    "drug",
    "raw_auc",
    "background_centered_auc",
    "matched_retention_pct",
    "n_matched",
    "valid_strata",
    "adequacy",
    "calibration_label",
    "brier",
    "expected_calibration_error",
    "decision_category",
    "recommended_action",
]

DRUG_SHORT_TO_LONG = {
    "Cipro": "Ciprofloxacin",
    "Norflox": "Norfloxacin",
    "Amox-Clav": "Amoxicillin-Clavulanic acid",
    "CRO": "Ceftriaxone",
    "CAZ": "Ceftazidime",
    "FEP": "Cefepime",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build deployment decision tables from audit and calibration outputs.")
    parser.add_argument("--audit-csv", default=DEFAULT_AUDIT)
    parser.add_argument("--calibration-csv", default=DEFAULT_CALIBRATION)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def parse_float(value: object) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def parse_auc(text: object) -> float | None:
    match = re.search(r"[-+]?\d*\.?\d+", str(text or ""))
    return parse_float(match.group(0)) if match else None


def drug_from_pair(pair: str) -> str:
    if "/" not in str(pair):
        return str(pair).strip()
    return str(pair).rsplit("/", 1)[-1].strip()


def calibration_lookup(calibration_rows: Sequence[dict]) -> dict[tuple[str, str], dict]:
    lookup = {}
    for row in calibration_rows:
        drug = str(row.get("drug", ""))
        lookup[(str(row.get("site", "")), drug)] = row
        short = next((key for key, value in DRUG_SHORT_TO_LONG.items() if value == drug), drug)
        lookup[(str(row.get("site", "")), short)] = row
    return lookup


def decision_category(
    *,
    raw_auc: float | None,
    centered_auc: float | None,
    adequacy: str,
    ece: float | None = None,
    brier: float | None = None,
) -> str:
    raw = raw_auc if raw_auc is not None and math.isfinite(raw_auc) else math.nan
    centered = centered_auc if centered_auc is not None and math.isfinite(centered_auc) else math.nan
    adequate = str(adequacy) == "interpretable"
    poor_calibration = (
        (ece is not None and math.isfinite(ece) and ece >= 0.15)
        or (brier is not None and math.isfinite(brier) and brier >= 0.25)
    )

    if not adequate or not math.isfinite(centered):
        return "insufficient_matched_evidence"
    if math.isfinite(raw) and raw < 0.60 and centered < 0.60:
        return "not_deployment_ready"
    if math.isfinite(raw) and raw >= 0.70 and centered >= 0.60:
        if poor_calibration:
            return "ranking_only_recalibrate_before_clinical_use"
        return "candidate_for_controlled_deployment"
    if math.isfinite(raw) and raw >= 0.70 and centered < 0.60:
        return "background_dependent_retrain_locally"
    if centered >= 0.60:
        return "promising_but_needs_local_validation"
    return "not_deployment_ready"


def action_for_category(category: str) -> str:
    return {
        "candidate_for_controlled_deployment": "Proceed only with local calibration check, locked thresholds, and ongoing drift monitoring.",
        "ranking_only_recalibrate_before_clinical_use": "Use only as a ranking signal until calibration, thresholds, and clinical operating points are locally reset.",
        "background_dependent_retrain_locally": "Do not deploy as portable focal-drug model; validate locally and retrain or recalibrate on current ecology.",
        "insufficient_matched_evidence": "Collect more matched support or additional external labels before making a deployment claim.",
        "promising_but_needs_local_validation": "Treat as a candidate requiring site-specific validation and calibration before deployment.",
        "not_deployment_ready": "Do not deploy; use for research triage or redesign the model/evaluation target.",
    }[category]


def decision_rules() -> list[dict]:
    return [
        {
            "scenario": "Raw high + background-centered high + acceptable calibration",
            "decision_category": "candidate_for_controlled_deployment",
            "interpretation": "Signal survives observable co-resistance background and probabilities are not obviously miscalibrated.",
            "recommended_action": action_for_category("candidate_for_controlled_deployment"),
        },
        {
            "scenario": "Raw high + background-centered high + poor calibration",
            "decision_category": "ranking_only_recalibrate_before_clinical_use",
            "interpretation": "Ranking may be useful, but probability thresholds are not ready for clinical decisions.",
            "recommended_action": action_for_category("ranking_only_recalibrate_before_clinical_use"),
        },
        {
            "scenario": "Raw high + background-centered low",
            "decision_category": "background_dependent_retrain_locally",
            "interpretation": "Raw performance likely depends on local resistant-population or co-resistance background.",
            "recommended_action": action_for_category("background_dependent_retrain_locally"),
        },
        {
            "scenario": "Raw high + underpowered matching",
            "decision_category": "insufficient_matched_evidence",
            "interpretation": "The audit cannot distinguish focal signal from background with the available matched strata.",
            "recommended_action": action_for_category("insufficient_matched_evidence"),
        },
        {
            "scenario": "Raw low + background-centered low",
            "decision_category": "not_deployment_ready",
            "interpretation": "Neither ordinary external performance nor background-controlled performance supports deployment.",
            "recommended_action": action_for_category("not_deployment_ready"),
        },
    ]


def build_readiness_rows(audit_rows: Sequence[dict], calibration_rows: Sequence[dict]) -> list[dict]:
    cal_lookup = calibration_lookup(calibration_rows)
    output = []
    for row in audit_rows:
        site = str(row.get("site", ""))
        pair = str(row.get("pair", ""))
        drug = drug_from_pair(pair)
        cal = cal_lookup.get((site, drug), {})
        raw_auc = parse_auc(row.get("raw_auc_95ci"))
        centered_auc = parse_auc(row.get("stratum_centered_auc_95ci"))
        ece = parse_float(cal.get("expected_calibration_error", ""))
        brier = parse_float(cal.get("brier", ""))
        category = decision_category(
            raw_auc=raw_auc,
            centered_auc=centered_auc,
            adequacy=str(row.get("adequacy", "")),
            ece=ece,
            brier=brier,
        )
        output.append(
            {
                "site": site,
                "pair": pair,
                "drug": drug,
                "raw_auc": raw_auc,
                "background_centered_auc": centered_auc,
                "matched_retention_pct": parse_float(row.get("matched_retention_pct", "")),
                "n_matched": parse_float(row.get("n_matched", "")),
                "valid_strata": parse_float(row.get("valid_strata", "")),
                "adequacy": row.get("adequacy", ""),
                "calibration_label": cal.get("calibration_label", ""),
                "brier": brier,
                "expected_calibration_error": ece,
                "decision_category": category,
                "recommended_action": action_for_category(category),
            }
        )
    return output


def format_value(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.6g}"
    if value is None:
        return ""
    return value


def write_csv(path: Path, rows: Sequence[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: format_value(row.get(field, "")) for field in fields})


def markdown_table(rows: Sequence[dict], fields: Sequence[str]) -> str:
    if not rows:
        return "_No rows._"
    rendered = [{field: str(format_value(row.get(field, ""))) for field in fields} for row in rows]
    widths = {field: max(len(field), *(len(row[field]) for row in rendered)) for field in fields}
    lines = [
        "| " + " | ".join(field.ljust(widths[field]) for field in fields) + " |",
        "| " + " | ".join("-" * widths[field] for field in fields) + " |",
    ]
    for row in rendered:
        lines.append("| " + " | ".join(row[field].ljust(widths[field]) for field in fields) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    with Path(args.audit_csv).open() as f:
        audit_rows = list(csv.DictReader(f))
    with Path(args.calibration_csv).open() as f:
        calibration_rows = list(csv.DictReader(f))

    rules = decision_rules()
    readiness = build_readiness_rows(audit_rows, calibration_rows)
    write_csv(output_dir / "deployment_decision_rules.csv", rules, RULE_FIELDS)
    write_csv(output_dir / "deployment_readiness_by_pair.csv", readiness, READINESS_FIELDS)
    md = [
        "# Deployment Decision Framework",
        "",
        "Action table for interpreting raw AUC, background-matched AUC, support, and calibration together.",
        "",
        "## Rules",
        "",
        markdown_table(rules, RULE_FIELDS),
        "",
        "## Pair-Level Readiness",
        "",
        markdown_table(readiness, READINESS_FIELDS),
        "",
    ]
    (output_dir / "deployment_decision_framework.md").write_text("\n".join(md))
    print(f"Wrote deployment decision framework to {output_dir}")


if __name__ == "__main__":
    main()
