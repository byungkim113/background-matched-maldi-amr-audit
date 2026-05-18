#!/usr/bin/env python3
"""Build a harmonized model-class audit matrix from background-audit summaries."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "outputs" / "analysis_outputs" / "model_class_matrix"

FIELDS = [
    "status",
    "model_class",
    "model_variant",
    "scope",
    "organism",
    "drug",
    "site",
    "raw_auc",
    "centered_auc",
    "pairwise_accuracy",
    "matched_retention",
    "n_total",
    "n_r",
    "n_matched",
    "n_matched_r",
    "valid_strata",
    "adequacy_label",
    "interpretation",
    "source_path",
    "notes",
]

DRUG_ALIASES = {
    "Amox-Clav": "Amoxicillin-Clavulanic acid",
    "Cipro": "Ciprofloxacin",
    "Norflox": "Norfloxacin",
    "CRO": "Ceftriaxone",
    "CAZ": "Ceftazidime",
    "FEP": "Cefepime",
}


@dataclass(frozen=True)
class ModelSpec:
    model_class: str
    model_variant: str
    organism: str
    drug: str
    summary_path: Path
    scope: str
    format: str = "audit_summary"
    notes: str = ""


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path)


def parse_float(value: object) -> float:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return math.nan
    return number if math.isfinite(number) else math.nan


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.6f}"
    text = str(value).strip()
    if text.upper() in {"NA", "NAN", "NONE"}:
        return ""
    number = parse_float(text)
    if math.isfinite(number) and any(ch in text for ch in ".eE"):
        return f"{number:.6f}"
    return text


def first_auc_token(value: str) -> str:
    text = str(value or "").strip()
    if not text or text.upper() == "NA":
        return ""
    return fmt(text.split()[0])


def missing_row(spec: ModelSpec, status: str = "missing") -> dict:
    return {
        "status": status,
        "model_class": spec.model_class,
        "model_variant": spec.model_variant,
        "scope": spec.scope,
        "organism": spec.organism,
        "drug": spec.drug,
        "site": "",
        "source_path": rel(spec.summary_path),
        "notes": spec.notes or "Run/export predictions, then audit with run_background_audit_framework.py.",
    }


def rows_from_audit_summary(spec: ModelSpec) -> list[dict]:
    if not spec.summary_path.exists():
        return [missing_row(spec)]
    rows = []
    with spec.summary_path.open(newline="") as f:
        for record in csv.DictReader(f):
            organism = record.get("organism", spec.organism)
            drug = record.get("drug", spec.drug)
            if organism != spec.organism or drug != spec.drug:
                continue
            rows.append(
                {
                    "status": "complete",
                    "model_class": spec.model_class,
                    "model_variant": spec.model_variant,
                    "scope": spec.scope,
                    "organism": organism,
                    "drug": drug,
                    "site": record.get("site", ""),
                    "raw_auc": fmt(record.get("raw_auc", "")),
                    "centered_auc": fmt(record.get("stratum_centered_auc", "")),
                    "pairwise_accuracy": fmt(record.get("pairwise_accuracy", "")),
                    "matched_retention": fmt(record.get("matched_retention", "")),
                    "n_total": fmt(record.get("n_total", "")),
                    "n_r": fmt(record.get("n_r", "")),
                    "n_matched": fmt(record.get("n_matched", "")),
                    "n_matched_r": fmt(record.get("n_matched_r", "")),
                    "valid_strata": fmt(record.get("n_valid_strata", "")),
                    "adequacy_label": record.get("adequacy_label", ""),
                    "interpretation": record.get("interpretation_category", ""),
                    "source_path": rel(spec.summary_path),
                    "notes": spec.notes,
                }
            )
    return rows or [missing_row(spec, status="no_matching_rows")]


def parse_pair(pair: str) -> tuple[str, str]:
    parts = [part.strip() for part in str(pair).split("/", 1)]
    if len(parts) != 2:
        return "", str(pair)
    organism, drug = parts
    organism = organism.replace("E. coli", "Escherichia coli")
    return organism, DRUG_ALIASES.get(drug, drug)


def rows_from_figure_table(spec: ModelSpec) -> list[dict]:
    if not spec.summary_path.exists():
        return [missing_row(spec)]
    rows = []
    with spec.summary_path.open(newline="") as f:
        for record in csv.DictReader(f):
            organism, drug = parse_pair(record.get("pair", ""))
            if organism != spec.organism or drug != spec.drug:
                continue
            retention = parse_float(record.get("matched_retention_pct", ""))
            rows.append(
                {
                    "status": "complete",
                    "model_class": spec.model_class,
                    "model_variant": spec.model_variant,
                    "scope": spec.scope,
                    "organism": organism,
                    "drug": drug,
                    "site": record.get("site", ""),
                    "raw_auc": first_auc_token(record.get("raw_auc_95ci", "")),
                    "centered_auc": first_auc_token(record.get("stratum_centered_auc_95ci", "")),
                    "pairwise_accuracy": "",
                    "matched_retention": fmt(retention / 100.0 if math.isfinite(retention) else math.nan),
                    "n_total": "",
                    "n_r": "",
                    "n_matched": fmt(record.get("n_matched", "")),
                    "n_matched_r": "",
                    "valid_strata": fmt(record.get("valid_strata", "")),
                    "adequacy_label": record.get("adequacy", ""),
                    "interpretation": record.get("interpretation", ""),
                    "source_path": rel(spec.summary_path),
                    "notes": spec.notes,
                }
            )
    return rows or [missing_row(spec, status="no_matching_rows")]


def read_spec_rows(spec: ModelSpec) -> list[dict]:
    if spec.format == "figure_table":
        return rows_from_figure_table(spec)
    return rows_from_audit_summary(spec)


def default_specs() -> list[ModelSpec]:
    return [
        ModelSpec(
            "CNN/Mega",
            "multi-task CNN",
            "Escherichia coli",
            "Ciprofloxacin",
            ROOT / "outputs/analysis_outputs/background_matched_transfer_audit_figure_table.csv",
            "E. coli expanded panel",
            format="figure_table",
        ),
        ModelSpec(
            "CNN/Mega",
            "multi-task CNN",
            "Escherichia coli",
            "Amoxicillin-Clavulanic acid",
            ROOT / "outputs/analysis_outputs/background_matched_transfer_audit_figure_table.csv",
            "E. coli expanded panel",
            format="figure_table",
        ),
        ModelSpec(
            "LightGBM",
            "multi-task",
            "Escherichia coli",
            "Ciprofloxacin",
            ROOT / "outputs/analysis_outputs/ecoli_lgbm_multi_background_audit/background_matched_audit_summary.csv",
            "E. coli expanded panel",
        ),
        ModelSpec(
            "LightGBM",
            "multi-task",
            "Escherichia coli",
            "Amoxicillin-Clavulanic acid",
            ROOT / "outputs/analysis_outputs/ecoli_lgbm_multi_background_audit/background_matched_audit_summary.csv",
            "E. coli expanded panel",
        ),
        ModelSpec(
            "LightGBM",
            "single-task",
            "Escherichia coli",
            "Ciprofloxacin",
            ROOT / "outputs/analysis_outputs/ecoli_lgbm_single_background_audit/background_matched_audit_summary.csv",
            "E. coli expanded panel",
        ),
        ModelSpec(
            "LightGBM",
            "single-task",
            "Escherichia coli",
            "Amoxicillin-Clavulanic acid",
            ROOT / "outputs/analysis_outputs/ecoli_lgbm_single_background_audit/background_matched_audit_summary.csv",
            "E. coli expanded panel",
        ),
        ModelSpec(
            "CNN/Mega",
            "multi-task CNN",
            "Staphylococcus aureus",
            "Oxacillin",
            ROOT / "outputs/analysis_outputs/saureus_panel_oxa_background_audit/background_matched_audit_summary.csv",
            "S. aureus/Oxacillin panel",
        ),
        ModelSpec(
            "LightGBM",
            "multi-task",
            "Staphylococcus aureus",
            "Oxacillin",
            ROOT / "outputs/analysis_outputs/saureus_lgbm_multi_oxa_background_audit/background_matched_audit_summary.csv",
            "S. aureus/Oxacillin panel",
        ),
        ModelSpec(
            "LightGBM",
            "single-task",
            "Staphylococcus aureus",
            "Oxacillin",
            ROOT / "outputs/analysis_outputs/saureus_lgbm_single_oxa_background_audit/background_matched_audit_summary.csv",
            "S. aureus/Oxacillin panel",
        ),
        ModelSpec(
            "Weis LR",
            "official published workflow",
            "Escherichia coli",
            "Ceftriaxone",
            ROOT / "outputs/analysis_outputs/weis_lr_official_panel_parity/audit_summaries/ecoli_ceftriaxone_background_matched_audit_summary.csv",
            "Weis official panel only",
            notes="Published-workflow compatibility evidence; not the expanded E. coli panel.",
        ),
        ModelSpec(
            "Weis LR",
            "official published workflow",
            "Klebsiella pneumoniae",
            "Ceftriaxone",
            ROOT / "outputs/analysis_outputs/weis_lr_official_panel_parity/audit_summaries/klebsiella_ceftriaxone_background_matched_audit_summary.csv",
            "Weis official panel only",
            notes="Published-workflow compatibility evidence.",
        ),
        ModelSpec(
            "Weis LR",
            "official published workflow",
            "Staphylococcus aureus",
            "Oxacillin",
            ROOT / "outputs/analysis_outputs/weis_lr_official_panel_parity/audit_summaries/saureus_oxacillin_background_matched_audit_summary.csv",
            "Weis official panel only",
            notes="Separate from the CNN/Mega S. aureus panel; weak DRIAMS-C raw signal.",
        ),
    ]


def build_matrix_rows(specs: Sequence[ModelSpec] | None = None) -> list[dict]:
    output = []
    for spec in specs or default_specs():
        output.extend(read_spec_rows(spec))
    return sorted(
        ({field: row.get(field, "") for field in FIELDS} for row in output),
        key=lambda row: (
            row["scope"],
            row["organism"],
            row["drug"],
            row["model_class"],
            row["model_variant"],
            row["site"],
        ),
    )


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: Sequence[dict]) -> str:
    if not rows:
        return "_No rows._"
    widths = {field: max(len(field), *(len(str(row.get(field, ""))) for row in rows)) for field in FIELDS}
    lines = [
        "| " + " | ".join(field.ljust(widths[field]) for field in FIELDS) + " |",
        "| " + " | ".join("-" * widths[field] for field in FIELDS) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")).ljust(widths[field]) for field in FIELDS) + " |")
    return "\n".join(lines)


def write_missing_commands(path: Path) -> None:
    path.write_text(
        """# Missing LGBM Model-Class Cells

Run this on Kaggle or another machine with the DRIAMS data mounted. The pipeline
exports isolate-level LGBM predictions, runs the model-agnostic background audit,
and rebuilds the matrix.

```bash
python scripts/run_model_class_matrix_pipeline.py \\
  --data-root /kaggle/input/datasets/drscarlat/driams \\
  --ecoli-run-dir /kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30 \\
  --saureus-run-dir /kaggle/working/runs/exp_saureus_panel_oxa_background_mae30
```

Use `--dry-run` first to print every command without running it.

The equivalent manual commands are below.

## E. coli LGBM single-task

```bash
python scripts/export_lgbm_predictions_for_audit.py \\
  --data-root /kaggle/input/datasets/drscarlat/driams \\
  --pair-profile ecoli_mechanism6 \\
  --run-dir /kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30 \\
  --variants single \\
  --output-dir /kaggle/working/lgbm_exports/ecoli

python run_background_audit_framework.py \\
  --predictions-csv /kaggle/working/lgbm_exports/ecoli/lgbm_single_predictions_long.csv \\
  --background-signature-col background_signature \\
  --model-name LGBM-single-ecoli6 \\
  --output-dir /kaggle/working/ecoli_lgbm_single_background_audit
```

## S. aureus/Oxacillin LGBM single-task and multi-task

```bash
python scripts/export_lgbm_predictions_for_audit.py \\
  --data-root /kaggle/input/datasets/drscarlat/driams \\
  --pair-profile saureus_panel \\
  --run-dir /kaggle/working/runs/exp_saureus_panel_oxa_background_mae30 \\
  --variants single,multi \\
  --train-if-missing \\
  --output-dir /kaggle/working/lgbm_exports/saureus_oxa

python run_background_audit_framework.py \\
  --predictions-csv /kaggle/working/lgbm_exports/saureus_oxa/lgbm_multi_predictions_long.csv \\
  --background-signature-col background_signature \\
  --model-name LGBM-multi-saureus-oxa \\
  --output-dir /kaggle/working/saureus_lgbm_multi_oxa_background_audit

python run_background_audit_framework.py \\
  --predictions-csv /kaggle/working/lgbm_exports/saureus_oxa/lgbm_single_predictions_long.csv \\
  --background-signature-col background_signature \\
  --model-name LGBM-single-saureus-oxa \\
  --output-dir /kaggle/working/saureus_lgbm_single_oxa_background_audit
```
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the model-class background-audit matrix.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = build_matrix_rows()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "model_class_matrix.csv", rows)
    (args.output_dir / "model_class_matrix.md").write_text(
        "# Model-Class Background-Audit Matrix\n\n" + markdown_table(rows) + "\n",
        encoding="utf-8",
    )
    write_missing_commands(args.output_dir / "run_missing_lgbm_commands.md")
    print(f"Wrote {args.output_dir / 'model_class_matrix.csv'}")
    print(f"Wrote {args.output_dir / 'model_class_matrix.md'}")


if __name__ == "__main__":
    main()
