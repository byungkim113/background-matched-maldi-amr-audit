#!/usr/bin/env python3
"""Export Mega_Model checkpoints to a model-agnostic long prediction CSV.

The Background-Matched MALDI-AMR Audit needs isolate-level prediction rows, not
only aggregated AUC tables. This script loads a completed Mega_Model run and
writes one row per isolate/drug prediction:

    isolate_id,site,year,organism,drug,label,prob,model_name

That CSV can then be passed to scripts/run_background_audit.py.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MEGA = ROOT / "Mega_Model.py"


def load_mega(path: Path):
    spec = importlib.util.spec_from_file_location("Mega_Model_export", path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(f"Could not load Mega_Model module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["Mega_Model_export"] = module
    spec.loader.exec_module(module)
    return module


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Export Mega_Model predictions for background audit.")
    p.add_argument("--mega-model", type=Path, default=DEFAULT_MEGA)
    p.add_argument("--run-dir", type=Path, required=True, help="Completed run dir containing config.json/models.")
    p.add_argument("--data-root", type=Path, default=None, help="Override data root from config.json.")
    p.add_argument("--output-csv", type=Path, default=None)
    p.add_argument("--model-name", default="Mega-CNN")
    p.add_argument("--tta-passes", type=int, default=None)
    p.add_argument("--no-bn-adapt", action="store_true")
    return p


def sample_to_row(sample, prob: float, mega, active_lookup: dict[int, tuple[str, str]], site: str, model_name: str):
    path_text, label, org_id = sample
    path = Path(path_text)
    organism, drug = active_lookup[int(org_id)]
    return {
        "model_name": model_name,
        "site": site,
        "year": path.parent.name,
        "isolate_id": path.stem,
        "organism": organism,
        "drug": drug,
        "label": int(label),
        "prob": float(prob),
    }


def main() -> None:
    args = build_parser().parse_args()
    import pandas as pd
    import torch

    config_path = args.run_dir / "config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Missing config.json: {config_path}")
    config = json.loads(config_path.read_text())

    mega = load_mega(args.mega_model)
    pair_profile = config.get("pair_profile", "ecoli_mechanism6")
    mega.init_config(pair_profile)

    data_root = Path(args.data_root or config.get("data_root") or mega.DATA_ROOT)
    temperature = float(config["temperature"])
    drug_conditioning = config.get("drug_conditioning", "task_id")
    tta_passes = int(args.tta_passes or config.get("tta_passes", mega.TTA_PASSES))
    use_bn_adapt = bool(config.get("use_bn_adapt", False)) and not args.no_bn_adapt
    active_pairs = [tuple(row) for row in config.get("active_pairs", [])]
    if not active_pairs:
        active_pairs = mega.screen_active_pairs(str(data_root), pair_profile)
    active_lookup = {int(oid): (org, drug) for oid, org, drug in active_pairs}

    selected = config.get("selected_seed_indices")
    if selected is None:
        selected = []
    checkpoint_dir = Path(config.get("ckpt_dir") or args.run_dir / "models")
    checkpoints = []
    if selected:
        checkpoints = [checkpoint_dir / f"maldi_amr_seed{int(seed)}.pt" for seed in selected]
    else:
        checkpoints = sorted(checkpoint_dir.glob("maldi_amr_seed*.pt"))
    checkpoints = [path for path in checkpoints if path.exists()]
    if not checkpoints:
        raise FileNotFoundError(f"No maldi_amr_seed*.pt checkpoints found in {checkpoint_dir}")

    models = []
    for checkpoint in checkpoints:
        state = torch.load(checkpoint, map_location=mega.DEVICE)
        if any(key.startswith("module.") for key in state):
            state = {key.removeprefix("module."): value for key, value in state.items()}
        model = mega.create_maldi_model(
            n_sites=mega.N_SITES,
            n_organisms=mega.N_ORGANISMS,
            drug_conditioning=drug_conditioning,
        ).to(mega.DEVICE)
        model.load_state_dict(state)
        model.eval()
        models.append(model)
        print(f"Loaded {checkpoint}", flush=True)

    output_csv = args.output_csv or (args.run_dir / "metrics" / "mega_predictions_long.csv")
    rows = []

    all_source = mega.load_all_organisms(str(data_root), mega.TRAIN_SITE, active_pairs)
    _, _, test_s = mega.make_split(all_source, active_pairs)
    eval_sites = {}
    if test_s:
        eval_sites[f"A-{mega.TEST_YEAR}"] = (test_s, False)
    for site in mega.TEST_SITES:
        if not (data_root / site).exists():
            continue
        site_rows = mega.load_all_organisms(str(data_root), site, active_pairs)
        site_samples = [(p, l, o) for p, l, o, _ in site_rows]
        if site_samples:
            eval_sites[site] = (site_samples, True)

    for site, (samples, is_external) in eval_sites.items():
        print(f"Predicting {site}: n={len(samples)}", flush=True)
        with mega.adapted_batchnorm(models, samples, use_adapt=is_external and use_bn_adapt):
            probs, labels, orgs = mega.ensemble_predict(models, samples, temperature, tta_passes)
        for sample, prob in zip(samples, probs):
            rows.append(sample_to_row(sample, float(prob), mega, active_lookup, site, args.model_name))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_csv, index=False)
    print(f"Wrote {output_csv} ({len(rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
