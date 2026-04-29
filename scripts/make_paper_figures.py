#!/usr/bin/env python3
"""Build final framework tables and paper-facing figures."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUILDER = ROOT / "scripts" / "make_final_framework_tables_figures.py"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build final paper tables and figures.")
    p.add_argument("--builder", type=Path, default=DEFAULT_BUILDER)
    return p


def main() -> None:
    args = build_parser().parse_args()
    cmd = [sys.executable, str(args.builder)]
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
