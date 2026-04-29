#!/usr/bin/env python3
"""Run the clinical4 Mega_Model experiment.

This is a thin reproducibility wrapper around Mega_Model.py. It intentionally
does not duplicate model code; it only pins the pair profile and common paper
defaults so the command used for the clinical4 result is explicit.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MEGA = ROOT / "Mega_Model.py"
DEFAULT_OUTPUT = ROOT / "runs"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run Mega_Model clinical4 training/evaluation.")
    p.add_argument("--mega-model", type=Path, default=DEFAULT_MEGA)
    p.add_argument("--data-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--experiment", default="exp_clinical4_mae30")
    p.add_argument("--mae-epochs", type=int, default=30)
    p.add_argument("--early-stop", choices=["macro", "primary"], default="macro")
    p.add_argument("--seed-policy", choices=["all", "threshold", "topk"], default="all")
    p.add_argument("--top-k-seeds", type=int, default=5)
    p.add_argument("--drug-conditioning", default="task_id")
    p.add_argument("--with-random-cv", action="store_true")
    p.add_argument("--with-ablation", action="store_true")
    p.add_argument("--with-leave-one-drug-out", action="store_true")
    p.add_argument("--no-lgbm", action="store_true")
    p.add_argument("--no-saliency", action="store_true")
    p.add_argument("--no-bn-adapt", action="store_true")
    p.add_argument("--no-prevalence-shift", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    cmd = [
        sys.executable,
        str(args.mega_model),
        "--pair-profile", "clinical4",
        "--experiment", args.experiment,
        "--data-root", str(args.data_root),
        "--output-dir", str(args.output_dir),
        "--mae-epochs", str(args.mae_epochs),
        "--early-stop", args.early_stop,
        "--seed-policy", args.seed_policy,
        "--top-k-seeds", str(args.top_k_seeds),
        "--drug-conditioning", args.drug_conditioning,
    ]
    if args.with_random_cv:
        cmd.append("--with-random-cv")
    if args.with_ablation:
        cmd.append("--with-ablation")
    if args.with_leave_one_drug_out:
        cmd.append("--with-leave-one-drug-out")
    if args.no_lgbm:
        cmd.append("--no-lgbm")
    if args.no_saliency:
        cmd.append("--no-saliency")
    if args.no_bn_adapt:
        cmd.append("--no-bn-adapt")
    if args.no_prevalence_shift:
        cmd.append("--no-prevalence-shift")

    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
