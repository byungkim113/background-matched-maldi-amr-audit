#!/usr/bin/env python3
"""Run the model-agnostic Background-Matched MALDI-AMR Audit."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AUDIT = ROOT / "run_background_audit_framework.py"
DEFAULT_OUTPUT = ROOT / "analysis_outputs" / "background_audit"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run background-matched audit on a long prediction CSV.")
    p.add_argument("--audit-script", type=Path, default=DEFAULT_AUDIT)
    p.add_argument("--predictions-csv", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--model-name", default="")
    p.add_argument("--min-pos-per-stratum", type=int, default=3)
    p.add_argument("--min-neg-per-stratum", type=int, default=3)
    p.add_argument("--bootstrap-n", type=int, default=500)
    p.add_argument("--permutation-n", type=int, default=500)
    p.add_argument("--match-year", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    cmd = [
        sys.executable,
        str(args.audit_script),
        "--predictions-csv", str(args.predictions_csv),
        "--output-dir", str(args.output_dir),
        "--min-pos-per-stratum", str(args.min_pos_per_stratum),
        "--min-neg-per-stratum", str(args.min_neg_per_stratum),
        "--bootstrap-n", str(args.bootstrap_n),
        "--permutation-n", str(args.permutation_n),
    ]
    if args.model_name:
        cmd += ["--model-name", args.model_name]
    if args.match_year:
        cmd.append("--match-year")
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
