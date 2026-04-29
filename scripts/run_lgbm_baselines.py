#!/usr/bin/env python3
"""Run Mega_Model's LightGBM baselines without CNN retraining.

Mega_Model.py can train LGBM baselines during the full CNN run, but this wrapper
calls the LGBM helper functions directly so a reviewer can reproduce the
classical-ML comparison faster.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MEGA = ROOT / "Mega_Model.py"
DEFAULT_OUTPUT = ROOT / "runs"


def load_mega(path: Path):
    spec = importlib.util.spec_from_file_location("Mega_Model_wrapped", path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load Mega_Model module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["Mega_Model_wrapped"] = module
    spec.loader.exec_module(module)
    return module


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run Mega_Model LightGBM baselines only.")
    p.add_argument("--mega-model", type=Path, default=DEFAULT_MEGA)
    p.add_argument("--data-root", type=Path, required=True)
    p.add_argument("--pair-profile", default="ecoli_mechanism6")
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--experiment", default="lgbm_ecoli_mechanism6")
    p.add_argument("--with-random-cv", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    import pandas as pd

    mega = load_mega(args.mega_model)
    mega.init_config(args.pair_profile)

    out_dir = args.output_dir / args.experiment
    metrics_dir = out_dir / "metrics"
    models_dir = out_dir / "models"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading active pairs for {args.pair_profile}", flush=True)
    active_pairs = mega.screen_active_pairs(str(args.data_root), args.pair_profile)
    all_source = mega.load_all_organisms(str(args.data_root), mega.TRAIN_SITE, active_pairs)
    train_s, val_s, test_s = mega.make_split(all_source, active_pairs)

    print("Training LightGBM single-task baselines", flush=True)
    lgbm_single = mega.train_lgbm_singletask(train_s, val_s, active_pairs)
    print("Training LightGBM multi-task baseline", flush=True)
    lgbm_multi = mega.train_lgbm_multitask(train_s, val_s)

    rows = []
    if test_s:
        rows.extend(mega.evaluate_lgbm_site(lgbm_single, lgbm_multi, test_s, active_pairs))
        for row in rows:
            if row.get("site") == "test":
                row["site"] = f"A-{mega.TEST_YEAR}"

    for site in mega.TEST_SITES:
        if not (args.data_root / site).exists():
            continue
        site_rows = mega.load_all_organisms(str(args.data_root), site, active_pairs)
        site_samples = [(p, l, o) for p, l, o, _ in site_rows]
        if not site_samples:
            continue
        site_eval = mega.evaluate_lgbm_site(lgbm_single, lgbm_multi, site_samples, active_pairs)
        for row in site_eval:
            row["site"] = site
        rows.extend(site_eval)

    pd.DataFrame(rows).to_csv(metrics_dir / "lgbm_results.csv", index=False)
    print(f"Wrote {metrics_dir / 'lgbm_results.csv'}", flush=True)

    if args.with_random_cv:
        print("Running random-CV inflation diagnostic", flush=True)
        random_rows = mega.run_random_cv_inflation_analysis(all_source, rows, active_pairs)
        pd.DataFrame(random_rows).to_csv(metrics_dir / "temporal_vs_random_cv.csv", index=False)
        print(f"Wrote {metrics_dir / 'temporal_vs_random_cv.csv'}", flush=True)


if __name__ == "__main__":
    main()
