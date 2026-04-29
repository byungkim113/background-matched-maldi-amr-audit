#!/usr/bin/env python3
"""Run the public Basel/Cuenod UPEC WGS-linked MALDI support analyses."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "data_manifests" / "upec_master_metadata.tsv"
DEFAULT_WGS_OUT = ROOT / "outputs" / "analysis_outputs" / "upec_wgs_validation_outputs"
DEFAULT_PROTEOMIC_OUT = ROOT / "outputs" / "analysis_outputs" / "updated_proteomic_overlap_outputs"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run UPEC WGS-linked MALDI and proteomic-overlap analyses.")
    p.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    p.add_argument("--median-peaks", type=Path, required=True)
    p.add_argument("--wgs-output-dir", type=Path, default=DEFAULT_WGS_OUT)
    p.add_argument("--proteomic-output-dir", type=Path, default=DEFAULT_PROTEOMIC_OUT)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--skip-proteomic", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    validation_script = ROOT / "scripts" / "upec_wgs_validation_analysis.py"
    proteomic_script = ROOT / "scripts" / "updated_proteomic_overlap_analysis.py"

    cmd = [
        sys.executable, str(validation_script),
        "--metadata", str(args.metadata),
        "--median-peaks", str(args.median_peaks),
        "--output-dir", str(args.wgs_output_dir),
        "--folds", str(args.folds),
    ]
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    if not args.skip_proteomic:
        cmd = [
            sys.executable, str(proteomic_script),
            "--metadata", str(args.metadata),
            "--median-peaks", str(args.median_peaks),
            "--output-dir", str(args.proteomic_output_dir),
            "--permutations", str(args.permutations),
        ]
        print("Running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
