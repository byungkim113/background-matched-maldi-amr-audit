#!/usr/bin/env python3
"""Run the missing LGBM audits and rebuild the model-class matrix.

Edit the DEFAULT_* paths below for your Kaggle session, or pass equivalent CLI
flags. The script assumes it is run from this repository, with DRIAMS mounted.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


# ─────────────────────────────────────────────────────────────────────────────
# Edit these defaults in Kaggle if your paths differ.
# ─────────────────────────────────────────────────────────────────────────────
DEFAULT_DATA_ROOT = Path("/kaggle/input/datasets/drscarlat/driams")
DEFAULT_ECOLI_RUN_DIR = Path("/kaggle/working/runs/exp_ecoli_mechanism6_drugid_mae30")
DEFAULT_SAUREUS_RUN_DIR = Path("/kaggle/working/runs/exp_saureus_panel_oxa_background_mae30")
DEFAULT_OUTPUT_ROOT = Path("outputs/analysis_outputs")
DEFAULT_BOOTSTRAP_N = 500
DEFAULT_PERMUTATION_N = 500
DEFAULT_TRAIN_IF_MISSING = True


@dataclass(frozen=True)
class PipelineConfig:
    repo_root: Path
    data_root: Path
    ecoli_run_dir: Path
    saureus_run_dir: Path
    output_root: Path
    bootstrap_n: int = DEFAULT_BOOTSTRAP_N
    permutation_n: int = DEFAULT_PERMUTATION_N
    train_if_missing: bool = DEFAULT_TRAIN_IF_MISSING


@dataclass(frozen=True)
class PipelineStep:
    name: str
    command: list[str]
    expected_output: Path | None = None


def under_repo(repo_root: Path, path: Path) -> Path:
    return path if path.is_absolute() else repo_root / path


def script_path(config: PipelineConfig, relative: str) -> str:
    return str(config.repo_root / relative)


def output_path(config: PipelineConfig, relative: str) -> Path:
    return under_repo(config.repo_root, config.output_root) / relative


def audit_command(
    config: PipelineConfig,
    *,
    predictions_csv: Path,
    model_name: str,
    output_dir: Path,
) -> list[str]:
    return [
        sys.executable,
        script_path(config, "run_background_audit_framework.py"),
        "--predictions-csv",
        str(predictions_csv),
        "--background-signature-col",
        "background_signature",
        "--model-name",
        model_name,
        "--output-dir",
        str(output_dir),
        "--bootstrap-n",
        str(config.bootstrap_n),
        "--permutation-n",
        str(config.permutation_n),
    ]


def export_lgbm_command(
    config: PipelineConfig,
    *,
    pair_profile: str,
    run_dir: Path,
    variants: str,
    output_dir: Path,
) -> list[str]:
    command = [
        sys.executable,
        script_path(config, "scripts/export_lgbm_predictions_for_audit.py"),
        "--data-root",
        str(config.data_root),
        "--pair-profile",
        pair_profile,
        "--run-dir",
        str(run_dir),
        "--variants",
        variants,
        "--output-dir",
        str(output_dir),
    ]
    if config.train_if_missing:
        command.append("--train-if-missing")
    return command


def build_pipeline_steps(config: PipelineConfig) -> list[PipelineStep]:
    ecoli_export_dir = output_path(config, "lgbm_prediction_exports/ecoli")
    saureus_export_dir = output_path(config, "lgbm_prediction_exports/saureus_oxa")

    ecoli_single_audit_dir = output_path(config, "ecoli_lgbm_single_background_audit")
    saureus_multi_audit_dir = output_path(config, "saureus_lgbm_multi_oxa_background_audit")
    saureus_single_audit_dir = output_path(config, "saureus_lgbm_single_oxa_background_audit")
    matrix_dir = output_path(config, "model_class_matrix")

    return [
        PipelineStep(
            "export_ecoli_lgbm_single",
            export_lgbm_command(
                config,
                pair_profile="ecoli_mechanism6",
                run_dir=config.ecoli_run_dir,
                variants="single",
                output_dir=ecoli_export_dir,
            ),
            ecoli_export_dir / "lgbm_single_predictions_long.csv",
        ),
        PipelineStep(
            "audit_ecoli_lgbm_single",
            audit_command(
                config,
                predictions_csv=ecoli_export_dir / "lgbm_single_predictions_long.csv",
                model_name="LGBM-single-ecoli6",
                output_dir=ecoli_single_audit_dir,
            ),
            ecoli_single_audit_dir / "background_matched_audit_summary.csv",
        ),
        PipelineStep(
            "export_saureus_lgbm_single_multi",
            export_lgbm_command(
                config,
                pair_profile="saureus_panel",
                run_dir=config.saureus_run_dir,
                variants="single,multi",
                output_dir=saureus_export_dir,
            ),
            saureus_export_dir / "lgbm_multi_predictions_long.csv",
        ),
        PipelineStep(
            "audit_saureus_lgbm_multi",
            audit_command(
                config,
                predictions_csv=saureus_export_dir / "lgbm_multi_predictions_long.csv",
                model_name="LGBM-multi-saureus-oxa",
                output_dir=saureus_multi_audit_dir,
            ),
            saureus_multi_audit_dir / "background_matched_audit_summary.csv",
        ),
        PipelineStep(
            "audit_saureus_lgbm_single",
            audit_command(
                config,
                predictions_csv=saureus_export_dir / "lgbm_single_predictions_long.csv",
                model_name="LGBM-single-saureus-oxa",
                output_dir=saureus_single_audit_dir,
            ),
            saureus_single_audit_dir / "background_matched_audit_summary.csv",
        ),
        PipelineStep(
            "rebuild_model_class_matrix",
            [
                sys.executable,
                script_path(config, "scripts/build_model_class_matrix.py"),
                "--output-dir",
                str(matrix_dir),
            ],
            matrix_dir / "model_class_matrix.csv",
        ),
    ]


def validate_config(config: PipelineConfig) -> None:
    required_files = [
        config.repo_root / "run_background_audit_framework.py",
        config.repo_root / "scripts/export_lgbm_predictions_for_audit.py",
        config.repo_root / "scripts/build_model_class_matrix.py",
    ]
    missing_files = [path for path in required_files if not path.exists()]
    if missing_files:
        raise FileNotFoundError("Missing repository files:\n" + "\n".join(str(path) for path in missing_files))
    if not config.data_root.exists():
        raise FileNotFoundError(f"DRIAMS data root does not exist: {config.data_root}")


def run_steps(config: PipelineConfig, steps: Sequence[PipelineStep], *, dry_run: bool, skip_existing: bool) -> None:
    for index, step in enumerate(steps, start=1):
        if skip_existing and step.expected_output and step.expected_output.exists():
            print(f"[{index}/{len(steps)}] SKIP {step.name}: {step.expected_output} exists", flush=True)
            continue
        print(f"\n[{index}/{len(steps)}] {step.name}", flush=True)
        print(shlex.join(step.command), flush=True)
        if dry_run:
            continue
        subprocess.run(step.command, cwd=config.repo_root, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LGBM exports/audits and rebuild model-class matrix.")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--ecoli-run-dir", type=Path, default=DEFAULT_ECOLI_RUN_DIR)
    parser.add_argument("--saureus-run-dir", type=Path, default=DEFAULT_SAUREUS_RUN_DIR)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP_N)
    parser.add_argument("--permutation-n", type=int, default=DEFAULT_PERMUTATION_N)
    parser.add_argument("--no-train-if-missing", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = PipelineConfig(
        repo_root=args.repo_root.resolve(),
        data_root=args.data_root,
        ecoli_run_dir=args.ecoli_run_dir,
        saureus_run_dir=args.saureus_run_dir,
        output_root=args.output_root,
        bootstrap_n=args.bootstrap_n,
        permutation_n=args.permutation_n,
        train_if_missing=not args.no_train_if_missing,
    )
    steps = build_pipeline_steps(config)
    if not args.dry_run:
        validate_config(config)
    run_steps(config, steps, dry_run=args.dry_run, skip_existing=args.skip_existing)
    print("\nPipeline finished.", flush=True)


if __name__ == "__main__":
    main()
