"""
Background-matched contrastive analysis for Mega_Model checkpoints.

Run this on Kaggle after an experiment has finished. It does not fit models.
It loads saved MALDI checkpoints, exports per-spectrum prediction scores, and
asks whether focal-drug prediction survives after matching on co-resistance
background excluding the focal drug.

Typical Kaggle usage:
  python background_matched_contrastive_kaggle.py \
    --mega-model-path /kaggle/working/Mega_Model.py \
    --data-root /kaggle/input/datasets/drscarlat/driams

The script auto-detects completed runs mounted under paths such as:
  /kaggle/input/newruns/runs/exp_ecoli_mechanism6_drugid_mae30
  /kaggle/input/datasets/bfdf121/newruns/runs/exp_ecoli_mechanism6_drugid_mae30

Outputs under /kaggle/working/background_matched_contrastive/ on Kaggle:
  background_matched_predictions.csv
  background_matched_retained_rows.csv
  background_matched_contrastive_summary.csv
  background_matched_contrastive_summary.md
  background_matched_sensitivity.csv
  background_matched_sensitivity.md
  background_matched_lgbm_single_predictions.csv
  background_matched_lgbm_multi_predictions.csv
"""

import argparse
import importlib.util
import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset


LABEL_MAP = {"S": 0, "R": 1, "I": -1, "-": -1}
LABEL_CHAR = {0: "S", 1: "R", -1: "U"}
DEFAULT_EXPERIMENT = "exp_ecoli_mechanism6_drugid_mae30"
DEFAULT_KAGGLE_OUTPUT = Path("/kaggle/working/background_matched_contrastive")
DEFAULT_BOOTSTRAP_N = 500
DEFAULT_PERMUTATION_N = 500
DEFAULT_STAT_SEED = 20260427


def load_mega_module(path):
    module_path = resolve_mega_model_path(path)
    spec = importlib.util.spec_from_file_location("Mega_Model_eval_only", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_config(run_dir):
    config_path = Path(run_dir) / "config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"config.json not found in {run_dir}")
    return json.loads(config_path.read_text())


def _unique_paths(paths):
    seen = set()
    unique = []
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def discover_mega_model_paths(input_root):
    root = Path(input_root)
    if not root.exists():
        return []
    discovered = []
    for top_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        discovered.extend(sorted(top_dir.glob("Mega_Model.py")))
        discovered.extend(sorted(top_dir.glob("*/Mega_Model.py")))
        if _is_probably_raw_data_mount(top_dir):
            continue
        discovered.extend(sorted(top_dir.rglob("Mega_Model.py")))
    return _unique_paths(discovered)


def candidate_mega_model_paths(requested_path):
    requested = Path(requested_path)
    candidates = [
        requested,
        Path.cwd() / requested,
        Path("/kaggle/working/Mega_Model.py"),
        Path("/kaggle/input/Mega_Model.py"),
    ]
    input_root = Path("/kaggle/input")
    if input_root.exists():
        candidates.extend(discover_mega_model_paths(input_root))
    return _unique_paths(candidates)


def resolve_mega_model_path(path):
    checked = candidate_mega_model_paths(path)
    for candidate in checked:
        if candidate.exists():
            if candidate != Path(path):
                print(f"Using Mega_Model.py: {candidate}")
            return candidate

    checked_text = "\n  ".join(str(path) for path in checked)
    raise FileNotFoundError(
        f"Mega_Model.py not found. Checked:\n  {checked_text}\n\n"
        f"{len(checked)} candidate Mega_Model.py paths checked.\n"
        "Please upload or copy Mega_Model.py to /kaggle/working/Mega_Model.py "
        "or pass --mega-model-path /path/to/Mega_Model.py."
    )


def _is_probably_raw_data_mount(path):
    name = path.name.lower()
    if name in {"datasets", "driams", "drscarlat"}:
        return True
    return (path / "drscarlat").exists() or (path / "DRIAMS-A").exists()


def discover_run_dirs(input_root, experiment_name=DEFAULT_EXPERIMENT):
    root = Path(input_root)
    if not root.exists():
        return []

    discovered = []
    top_dirs = [p for p in sorted(root.iterdir()) if p.is_dir()]
    for top_dir in top_dirs:
        discovered.extend(sorted(top_dir.glob(f"runs/{experiment_name}")))
        discovered.extend(sorted(top_dir.glob(f"*/runs/{experiment_name}")))
        discovered.extend(sorted(top_dir.glob(f"*/newruns/runs/{experiment_name}")))
        discovered.extend(sorted(top_dir.glob(f"*/*/runs/{experiment_name}")))
        if _is_probably_raw_data_mount(top_dir):
            continue
        for config_path in sorted(top_dir.rglob("config.json")):
            run_dir = config_path.parent
            if run_dir.name == experiment_name or experiment_name in str(run_dir):
                discovered.append(run_dir)

    # Last-resort fallback: useful when Kaggle nests the uploaded run oddly.
    # This can be slower, so it only runs after the small mounts were checked.
    if not discovered:
        for config_path in sorted(root.rglob("config.json")):
            run_dir = config_path.parent
            if run_dir.name == experiment_name or experiment_name in str(run_dir):
                discovered.append(run_dir)
    return _unique_paths(discovered)


def candidate_run_dirs(experiment_name=DEFAULT_EXPERIMENT):
    candidates = [
        Path("/kaggle/working/runs") / experiment_name,
        Path("/kaggle/input/newruns/runs") / experiment_name,
        Path("/kaggle/input/datasets/bfdf121/newruns/runs") / experiment_name,
    ]
    input_root = Path("/kaggle/input")
    if input_root.exists():
        candidates.extend(sorted(input_root.glob(f"*/runs/{experiment_name}")))
        candidates.extend(sorted(input_root.glob(f"datasets/*/newruns/runs/{experiment_name}")))
        candidates.extend(discover_run_dirs(input_root, experiment_name))
    return _unique_paths(candidates)


def resolve_run_dir(run_dir, experiment_name=DEFAULT_EXPERIMENT):
    requested = Path(run_dir)
    checked = [requested]
    if (requested / "config.json").exists():
        return requested

    for candidate in candidate_run_dirs(experiment_name):
        if candidate in checked:
            continue
        checked.append(candidate)
        if (candidate / "config.json").exists():
            print(f"Using run directory: {candidate}")
            return candidate

    checked_text = "\n  ".join(str(path) for path in checked)
    raise FileNotFoundError(
        f"config.json not found. Checked:\n  {checked_text}\n\n"
        f"{len(checked)} candidate run directories discovered.\n"
        "If your run is mounted in Kaggle, pass --run-dir "
        "/kaggle/input/<dataset-name>/runs/<experiment-name>."
    )


def default_output_dir(run_dir):
    if Path("/kaggle/working").exists():
        return DEFAULT_KAGGLE_OUTPUT
    return Path(run_dir) / "background_matched_contrastive"


def looks_like_driams_root(path):
    root = Path(path)
    return (
        (root / "DRIAMS-A" / "binned_6000").exists()
        and (root / "DRIAMS-A" / "id").exists()
    )


def discover_data_roots(input_root):
    root = Path(input_root)
    if not root.exists():
        return []
    discovered = []
    for candidate in [
        root / "datasets" / "drscarlat" / "driams",
        root / "drscarlat" / "driams",
        root / "driams",
    ]:
        if looks_like_driams_root(candidate):
            discovered.append(candidate)
    for binned_dir in sorted(root.rglob("DRIAMS-A/binned_6000")):
        site_dir = binned_dir.parent
        data_root = site_dir.parent
        if looks_like_driams_root(data_root):
            discovered.append(data_root)
    return _unique_paths(discovered)


def candidate_data_roots(requested_root):
    requested = Path(requested_root)
    candidates = [
        requested,
        Path("/kaggle/input/datasets/drscarlat/driams"),
        Path("/kaggle/input/drscarlat/driams"),
        Path("/kaggle/input/driams"),
    ]
    input_root = Path("/kaggle/input")
    if input_root.exists():
        candidates.extend(discover_data_roots(input_root))
    return _unique_paths(candidates)


def resolve_data_root(data_root):
    checked = candidate_data_roots(data_root)
    for candidate in checked:
        if looks_like_driams_root(candidate):
            if candidate != Path(data_root):
                print(f"Using DRIAMS data root: {candidate}")
            return candidate
    checked_text = "\n  ".join(str(path) for path in checked)
    raise FileNotFoundError(
        f"DRIAMS data root not found. Checked:\n  {checked_text}\n\n"
        f"{len(checked)} candidate DRIAMS data roots checked. "
        "Pass --data-root /path/containing/DRIAMS-A."
    )


def normalize_name(text):
    return "".join(ch.lower() for ch in str(text) if ch.isalnum())


def find_column(df, target):
    if target in df.columns:
        return target
    target_norm = normalize_name(target)
    for col in df.columns:
        if normalize_name(col) == target_norm:
            return col
    for col in df.columns:
        col_norm = normalize_name(col)
        if target_norm in col_norm or col_norm in target_norm:
            return col
    return None


def find_id_column(df):
    for col in df.columns:
        if col.lower() in ("sample_id", "code", "uuid", "id"):
            return col
    return df.columns[0]


def find_species_column(df):
    for col in df.columns:
        lower = col.lower()
        if "species" in lower or "organism" in lower:
            return col
    return None


def read_site_clean_table(data_root, site):
    id_root = Path(data_root) / site / "id"
    frames = []
    for year_dir in sorted(id_root.iterdir()):
        clean_path = year_dir / f"{year_dir.name}_clean.csv"
        if clean_path.exists():
            frame = pd.read_csv(clean_path, low_memory=False)
            frame["__source_year"] = str(year_dir.name)
            frames.append(frame)
    if not frames:
        raise FileNotFoundError(f"No *_clean.csv files found under {id_root}")
    return pd.concat(frames, ignore_index=True)


def value_to_label(value):
    return LABEL_MAP.get(str(value).strip(), -1)


def make_background_signature(label_by_drug, focal_drug, background_drugs,
                              exclude_focal=True):
    parts = []
    known = 0
    resistant = 0
    for drug in background_drugs:
        if exclude_focal and drug == focal_drug:
            continue
        label = int(label_by_drug.get(drug, -1))
        if label in (0, 1):
            known += 1
            resistant += int(label == 1)
        parts.append(f"{drug}={LABEL_CHAR.get(label, 'U')}")
    return "|".join(parts), known, resistant


def resolve_eval_site(site_label, train_site, test_year):
    if site_label == f"A-{test_year}":
        return train_site, test_year
    return site_label, None


def build_prediction_rows(mega, data_root, site_label, active_pairs,
                          panel_drugs_by_organism, train_site, test_year,
                          match_year=False):
    site, year_filter = resolve_eval_site(site_label, train_site, test_year)
    spectrum_index = mega.build_spectrum_index(data_root, site)
    if spectrum_index is None:
        print(f"{site_label}: diagnostics spectrum_files=0 reason=binned_6000_missing site={site}")
        return []
    _, stem_to_paths = spectrum_index
    df = read_site_clean_table(data_root, site)
    id_col = find_id_column(df)
    species_col = find_species_column(df)
    rows = []

    active_by_organism = {}
    for org_id, organism, drug in active_pairs:
        active_by_organism.setdefault(organism, []).append((org_id, drug))

    for organism, org_drugs in active_by_organism.items():
        org_df = df
        if species_col is not None:
            org_df = org_df[
                org_df[species_col].astype(str).str.contains(
                    organism, case=False, na=False)
            ]
        if year_filter is not None:
            org_df = org_df[org_df["__source_year"].astype(str) == str(year_filter)]
        if org_df.empty:
            print(
                f"{site_label}: diagnostics organism={organism} "
                f"metadata_rows=0 species_col={species_col} year_filter={year_filter}"
            )
            continue

        background_drugs = panel_drugs_by_organism[organism]
        drug_cols = {drug: find_column(org_df, drug) for drug in background_drugs}
        rows_with_spectrum = 0
        valid_label_counts = {drug: 0 for drug in background_drugs}
        for _, record in org_df.iterrows():
            uid = str(record[id_col]).strip()
            paths = stem_to_paths.get(uid, [])
            if not paths:
                continue
            rows_with_spectrum += 1
            label_by_drug = {}
            for drug, col in drug_cols.items():
                label_by_drug[drug] = value_to_label(record[col]) if col else -1
                valid_label_counts[drug] += int(label_by_drug[drug] in (0, 1))
            for org_id, focal_drug in org_drugs:
                focal_label = int(label_by_drug.get(focal_drug, -1))
                if focal_label not in (0, 1):
                    continue
                signature, known_count, resistant_count = make_background_signature(
                    label_by_drug, focal_drug, background_drugs, exclude_focal=True)
                for path in paths:
                    rows.append(dict(
                        site=site_label,
                        raw_site=site,
                        year=str(record["__source_year"]),
                        uid=uid,
                        sample_path=str(path),
                        organism=organism,
                        drug=focal_drug,
                        org_id=int(org_id),
                        label=focal_label,
                        background_signature=signature,
                        background_known_count=known_count,
                        background_resistant_count=resistant_count,
                        match_year=bool(match_year),
                    ))
        print(
            f"{site_label}: diagnostics organism={organism} "
            f"spectrum_files={len(stem_to_paths)} metadata_rows={len(org_df)} "
            f"rows_with_spectrum={rows_with_spectrum} id_col={id_col} "
            f"species_col={species_col} drug_columns={drug_cols} "
            f"valid_label_counts={valid_label_counts}"
        )
    return rows


class PredictionDataset(Dataset):
    def __init__(self, rows, mega, use_augment=False):
        self.rows = rows
        self.mega = mega
        self.use_augment = use_augment

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, idx):
        row = self.rows[idx]
        x = self.mega.load_spectrum(row["sample_path"])
        if self.use_augment:
            x = self.mega.augment(x)
        x = torch.from_numpy(x).unsqueeze(0)
        org = torch.tensor(int(row["org_id"]), dtype=torch.long)
        return x, org


def load_checkpoint_models(mega, run_dir, config, ckpt_dir=None):
    ckpt_root = Path(ckpt_dir) if ckpt_dir else Path(run_dir) / "models"
    if not ckpt_root.exists():
        raise FileNotFoundError(f"Checkpoint directory not found: {ckpt_root}")
    seed_indices = config.get("selected_seed_indices") or list(range(config.get("n_seeds", 5)))
    drug_conditioning = config.get("drug_conditioning", "task_id")
    models = []
    for seed in seed_indices:
        ckpt = ckpt_root / f"maldi_amr_seed{int(seed)}.pt"
        if not ckpt.exists():
            raise FileNotFoundError(f"Missing checkpoint: {ckpt}")
        state = torch.load(ckpt, map_location=mega.DEVICE)
        if any(str(k).startswith("module.") for k in state):
            state = {str(k).removeprefix("module."): v for k, v in state.items()}
        model = mega.create_maldi_model(
            n_sites=mega.N_SITES,
            n_organisms=mega.N_ORGANISMS,
            drug_conditioning=drug_conditioning,
        ).to(mega.DEVICE)
        model.load_state_dict(state)
        model.eval()
        models.append(model)
    if not models:
        raise RuntimeError("No checkpoint models were loaded.")
    return models


@torch.no_grad()
def score_rows(mega, rows, models, temperature, tta_passes=5, batch_size=256):
    if not rows:
        return rows
    all_model_probs = []
    n_passes = max(1, int(tta_passes))
    for model in models:
        pass_probs = []
        for pass_idx in range(n_passes):
            dataset = PredictionDataset(
                rows, mega, use_augment=(n_passes > 1))
            loader = DataLoader(
                dataset, batch_size=batch_size, shuffle=False,
                num_workers=0, pin_memory=torch.cuda.is_available())
            probs = []
            for x, org in loader:
                logits = model(x.to(mega.DEVICE), org.to(mega.DEVICE)) / float(temperature)
                probs.append(torch.sigmoid(logits).detach().cpu().numpy())
            pass_probs.append(np.concatenate(probs))
        all_model_probs.append(np.mean(pass_probs, axis=0))
    mean_probs = np.mean(all_model_probs, axis=0)
    for row, prob in zip(rows, mean_probs):
        row["prob"] = float(prob)
    return rows


def load_lgbm_models(run_dir, model_dir=None, score_family="lgbm_multi"):
    model_root = Path(model_dir) if model_dir else Path(run_dir) / "models"
    if not model_root.exists():
        raise FileNotFoundError(f"LightGBM model directory not found: {model_root}")

    if score_family == "lgbm_multi":
        model_path = model_root / "lgbm_multi.pkl"
        if not model_path.exists():
            raise FileNotFoundError(f"Missing LightGBM multi-task model: {model_path}")
        with open(model_path, "rb") as f:
            return pickle.load(f)

    if score_family == "lgbm_single":
        models = {}
        for path in sorted(model_root.glob("lgbm_single_org*.pkl")):
            stem = path.stem
            org_id = int(stem.replace("lgbm_single_org", ""))
            with open(path, "rb") as f:
                models[org_id] = pickle.load(f)
        if not models:
            raise FileNotFoundError(f"No lgbm_single_org*.pkl files found in {model_root}")
        return models

    raise ValueError(f"Unsupported LightGBM score family: {score_family}")


def score_lgbm_rows(mega, rows, lgbm_models, score_family="lgbm_multi",
                    batch_size=4096):
    if not rows:
        return rows

    frame = pd.DataFrame(rows).copy()
    frame["prob"] = float("nan")
    unique_paths = frame["sample_path"].drop_duplicates().tolist()
    path_to_pos = {path: i for i, path in enumerate(unique_paths)}
    spectra = []
    for path in unique_paths:
        spectra.append(mega.load_spectrum(path))
    x_unique = np.stack(spectra).astype(np.float32)
    row_positions = frame["sample_path"].map(path_to_pos).to_numpy(dtype=int)

    for org_id, group in frame.groupby("org_id", sort=False):
        row_idx = group.index.to_numpy()
        x = x_unique[row_positions[row_idx]]
        if score_family == "lgbm_single":
            org_id_int = int(org_id)
            if org_id_int not in lgbm_models:
                raise KeyError(f"No single-task LightGBM model for org_id={org_id_int}")
            model = lgbm_models[org_id_int]
            probs = model.predict(x)
        elif score_family == "lgbm_multi":
            model = lgbm_models
            org_feature = np.full((len(x), 1), float(org_id), dtype=np.float32)
            probs = model.predict(np.hstack([x, org_feature]))
        else:
            raise ValueError(f"Unsupported score family: {score_family}")
        frame.loc[row_idx, "prob"] = np.asarray(probs, dtype=float)

    if frame["prob"].isna().any():
        raise RuntimeError("LightGBM scoring left some rows without probabilities.")
    return frame.to_dict("records")


def safe_auc(labels, scores):
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    if len(labels) == 0 or labels.sum() == 0 or labels.sum() == len(labels):
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    ranks = np.empty(len(scores), dtype=float)
    i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and sorted_scores[j + 1] == sorted_scores[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def safe_aupr(labels, scores):
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    n_pos = int(labels.sum())
    if n_pos == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    tp = np.cumsum(y)
    precision = tp / (np.arange(len(y)) + 1)
    return float(precision[y == 1].sum() / n_pos)


def pairwise_accuracy_for_strata(df, stratum_cols):
    wins = 0.0
    total = 0
    for _, group in df.groupby(stratum_cols, dropna=False):
        pos = group[group["label"] == 1]["prob"].to_numpy(dtype=float)
        neg = group[group["label"] == 0]["prob"].to_numpy(dtype=float)
        if len(pos) == 0 or len(neg) == 0:
            continue
        diff = pos[:, None] - neg[None, :]
        wins += float((diff > 0).sum()) + 0.5 * float((diff == 0).sum())
        total += int(diff.size)
    return wins / total if total else float("nan"), total


def bootstrap_metric_ci(labels, scores, metric_fn, n_boot=DEFAULT_BOOTSTRAP_N,
                        seed=DEFAULT_STAT_SEED, ci=0.95):
    if n_boot <= 0:
        return float("nan"), float("nan")
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)
    mask = np.isfinite(scores)
    labels = labels[mask]
    scores = scores[mask]
    if len(labels) == 0 or labels.sum() == 0:
        return float("nan"), float("nan")

    rng = np.random.default_rng(int(seed))
    values = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, len(labels), len(labels))
        value = metric_fn(labels[idx], scores[idx])
        if not math.isnan(value):
            values.append(value)
    if not values:
        return float("nan"), float("nan")

    alpha = (1.0 - float(ci)) / 2.0
    low, high = np.quantile(np.asarray(values, dtype=float), [alpha, 1.0 - alpha])
    return float(low), float(high)


def permutation_null_within_strata(matched_df, stratum_cols, score_col="centered_prob",
                                  observed=None, n_perm=DEFAULT_PERMUTATION_N,
                                  seed=DEFAULT_STAT_SEED):
    if n_perm <= 0 or matched_df.empty:
        return float("nan"), float("nan"), float("nan")

    frame = matched_df.reset_index(drop=True)
    labels = frame["label"].to_numpy(dtype=int)
    scores = frame[score_col].to_numpy(dtype=float)
    observed_auc = safe_auc(labels, scores) if observed is None else float(observed)
    if math.isnan(observed_auc):
        return float("nan"), float("nan"), float("nan")

    group_indices = [
        group.index.to_numpy(dtype=int)
        for _, group in frame.groupby(stratum_cols, dropna=False)
    ]
    if not group_indices:
        return float("nan"), float("nan"), float("nan")

    rng = np.random.default_rng(int(seed))
    null_values = []
    for _ in range(int(n_perm)):
        shuffled = labels.copy()
        for idx in group_indices:
            shuffled[idx] = rng.permutation(shuffled[idx])
        value = safe_auc(shuffled, scores)
        if not math.isnan(value):
            null_values.append(value)
    if not null_values:
        return float("nan"), float("nan"), float("nan")

    null = np.asarray(null_values, dtype=float)
    p_value = (1.0 + float((null >= observed_auc).sum())) / (1.0 + float(len(null)))
    return float(p_value), float(null.mean()), float(null.std(ddof=1) if len(null) > 1 else 0.0)


def assign_adequacy_label(n_matched, matched_retention, n_valid_strata,
                          pairwise_comparisons, min_n_matched=100,
                          min_retention=0.10, min_pairwise=100):
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


def compute_background_matched_summary(pred_df, min_pos_per_stratum=3,
                                       min_neg_per_stratum=3,
                                       match_year=False,
                                       bootstrap_n=DEFAULT_BOOTSTRAP_N,
                                       permutation_n=DEFAULT_PERMUTATION_N,
                                       random_seed=DEFAULT_STAT_SEED,
                                       adequacy_min_n_matched=100,
                                       adequacy_min_retention=0.10,
                                       adequacy_min_pairwise=100):
    rows = []
    retained_frames = []
    stratum_cols = ["background_signature"]
    if match_year:
        stratum_cols.append("year")

    for group_idx, ((site, organism, drug), pair_df) in enumerate(
            pred_df.groupby(["site", "organism", "drug"])):
        raw_auc = safe_auc(pair_df["label"], pair_df["prob"])
        raw_aupr = safe_aupr(pair_df["label"], pair_df["prob"])
        seed_base = int(random_seed) + int(group_idx) * 100003
        raw_auc_ci_low, raw_auc_ci_high = bootstrap_metric_ci(
            pair_df["label"], pair_df["prob"], safe_auc,
            n_boot=bootstrap_n, seed=seed_base + 1)

        valid_indices = []
        for key, group in pair_df.groupby(stratum_cols, dropna=False):
            n_pos = int((group["label"] == 1).sum())
            n_neg = int((group["label"] == 0).sum())
            if n_pos >= min_pos_per_stratum and n_neg >= min_neg_per_stratum:
                valid_indices.extend(group.index.tolist())

        if valid_indices:
            matched = pair_df.loc[valid_indices].copy()
            matched["centered_prob"] = matched["prob"] - matched.groupby(
                stratum_cols, dropna=False)["prob"].transform("mean")
            matched_auc = safe_auc(matched["label"], matched["prob"])
            matched_aupr = safe_aupr(matched["label"], matched["prob"])
            stratum_centered_auc = safe_auc(matched["label"], matched["centered_prob"])
            matched_auc_ci_low, matched_auc_ci_high = bootstrap_metric_ci(
                matched["label"], matched["prob"], safe_auc,
                n_boot=bootstrap_n, seed=seed_base + 2)
            centered_auc_ci_low, centered_auc_ci_high = bootstrap_metric_ci(
                matched["label"], matched["centered_prob"], safe_auc,
                n_boot=bootstrap_n, seed=seed_base + 3)
            permutation_p, permutation_null_mean, permutation_null_std = (
                permutation_null_within_strata(
                    matched, stratum_cols, score_col="centered_prob",
                    observed=stratum_centered_auc, n_perm=permutation_n,
                    seed=seed_base + 4))
            pair_acc, pair_n = pairwise_accuracy_for_strata(matched, stratum_cols)
            matched_retention = len(matched) / len(pair_df)
            matched["matched_valid_stratum"] = True
            n_valid_strata = matched.groupby(stratum_cols, dropna=False).ngroups
            retained_frames.append(matched)
        else:
            matched = pd.DataFrame()
            matched_auc = float("nan")
            matched_aupr = float("nan")
            stratum_centered_auc = float("nan")
            matched_auc_ci_low = float("nan")
            matched_auc_ci_high = float("nan")
            centered_auc_ci_low = float("nan")
            centered_auc_ci_high = float("nan")
            permutation_p = float("nan")
            permutation_null_mean = float("nan")
            permutation_null_std = float("nan")
            pair_acc = float("nan")
            pair_n = 0
            matched_retention = 0.0
            n_valid_strata = 0

        adequacy_label = assign_adequacy_label(
            n_matched=len(matched),
            matched_retention=matched_retention,
            n_valid_strata=n_valid_strata,
            pairwise_comparisons=pair_n,
            min_n_matched=adequacy_min_n_matched,
            min_retention=adequacy_min_retention,
            min_pairwise=adequacy_min_pairwise,
        )

        rows.append(dict(
            site=site,
            organism=organism,
            drug=drug,
            raw_auc=raw_auc,
            raw_auc_ci_low=raw_auc_ci_low,
            raw_auc_ci_high=raw_auc_ci_high,
            raw_aupr=raw_aupr,
            matched_auc=matched_auc,
            matched_auc_ci_low=matched_auc_ci_low,
            matched_auc_ci_high=matched_auc_ci_high,
            matched_aupr=matched_aupr,
            stratum_centered_auc=stratum_centered_auc,
            stratum_centered_auc_ci_low=centered_auc_ci_low,
            stratum_centered_auc_ci_high=centered_auc_ci_high,
            pairwise_accuracy=pair_acc,
            pairwise_comparisons=pair_n,
            permutation_p=permutation_p,
            permutation_null_mean=permutation_null_mean,
            permutation_null_std=permutation_null_std,
            matched_retention=matched_retention,
            adequacy_label=adequacy_label,
            n_total=len(pair_df),
            n_r=int(pair_df["label"].sum()),
            n_matched=len(matched),
            n_matched_r=int(matched["label"].sum()) if not matched.empty else 0,
            n_valid_strata=n_valid_strata,
            min_pos_per_stratum=min_pos_per_stratum,
            min_neg_per_stratum=min_neg_per_stratum,
            match_year=bool(match_year),
        ))

    retained = pd.concat(retained_frames, ignore_index=True) if retained_frames else pd.DataFrame()
    return pd.DataFrame(rows), retained


def parse_sensitivity_thresholds(text):
    if text is None or str(text).strip() == "":
        return []
    thresholds = []
    for part in str(text).split(","):
        part = part.strip()
        if not part:
            continue
        value = int(part)
        if value <= 0:
            raise ValueError("--sensitivity-thresholds values must be positive integers")
        thresholds.append(value)
    return sorted(set(thresholds))


def run_sensitivity_summaries(pred_df, thresholds, match_year=False,
                              random_seed=DEFAULT_STAT_SEED,
                              adequacy_min_n_matched=100,
                              adequacy_min_retention=0.10,
                              adequacy_min_pairwise=100):
    frames = []
    for threshold in thresholds:
        summary, _ = compute_background_matched_summary(
            pred_df,
            min_pos_per_stratum=threshold,
            min_neg_per_stratum=threshold,
            match_year=match_year,
            bootstrap_n=0,
            permutation_n=0,
            random_seed=random_seed + threshold * 1009,
            adequacy_min_n_matched=adequacy_min_n_matched,
            adequacy_min_retention=adequacy_min_retention,
            adequacy_min_pairwise=adequacy_min_pairwise,
        )
        summary.insert(0, "sensitivity_threshold", threshold)
        frames.append(summary)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def output_names_for_score_family(score_family):
    if score_family == "cnn":
        return dict(
            predictions="background_matched_predictions.csv",
            retained="background_matched_retained_rows.csv",
            summary="background_matched_contrastive_summary.csv",
            summary_md="background_matched_contrastive_summary.md",
            sensitivity="background_matched_sensitivity.csv",
            sensitivity_md="background_matched_sensitivity.md",
        )
    return dict(
        predictions=f"background_matched_{score_family}_predictions.csv",
        retained=f"background_matched_{score_family}_retained_rows.csv",
        summary=f"background_matched_{score_family}_summary.csv",
        summary_md=f"background_matched_{score_family}_summary.md",
        sensitivity=f"background_matched_{score_family}_sensitivity.csv",
        sensitivity_md=f"background_matched_{score_family}_sensitivity.md",
    )


def markdown_table(df):
    if df.empty:
        return "_No rows._"
    frame = df.copy()
    for col in frame.columns:
        if pd.api.types.is_float_dtype(frame[col]):
            frame[col] = frame[col].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
        else:
            frame[col] = frame[col].map(lambda x: "" if pd.isna(x) else str(x))
    widths = {col: max(len(col), *(len(v) for v in frame[col].tolist())) for col in frame.columns}
    header = "| " + " | ".join(col.ljust(widths[col]) for col in frame.columns) + " |"
    sep = "| " + " | ".join("-" * widths[col] for col in frame.columns) + " |"
    lines = [header, sep]
    for row in frame.to_dict("records"):
        lines.append("| " + " | ".join(row[col].ljust(widths[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


def write_markdown(summary_df, output_path):
    lines = [
        "# Background-Matched Contrastive Analysis\n\n",
        "This evaluation loads saved checkpoints and asks whether focal-drug prediction ",
        "survives after retaining only co-resistance background strata that contain both ",
        "resistant and susceptible isolates for the focal drug. The background signature ",
        "excludes the focal drug itself.\n\n",
        "Interpretation: if matched/stratum-centered AUC collapses relative to raw AUC, ",
        "the original model likely relied on population or co-resistance background. ",
        "If it remains high, there may be focal-drug-associated spectral signal.\n\n",
        markdown_table(summary_df),
        "\n",
    ]
    Path(output_path).write_text("".join(lines))


def main():
    parser = argparse.ArgumentParser(description="Evaluation-only background-matched contrastive analysis")
    parser.add_argument("--mega-model-path", default="Mega_Model.py")
    parser.add_argument("--run-dir", default=f"/kaggle/working/runs/{DEFAULT_EXPERIMENT}")
    parser.add_argument("--ckpt-dir", default=None)
    parser.add_argument("--data-root", default="/kaggle/input/datasets/drscarlat/driams")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--eval-sites", default="A-2018,DRIAMS-B,DRIAMS-C,DRIAMS-D")
    parser.add_argument("--pair-profile", default=None)
    parser.add_argument(
        "--score-family", default="cnn",
        choices=["cnn", "lgbm_single", "lgbm_multi"],
        help="Which saved model family to score before running the matched audit.")
    parser.add_argument("--lgbm-model-dir", default=None)
    parser.add_argument(
        "--row-template-csv", default=None,
        help="Optional existing background_matched_predictions.csv to reuse rows/background strata but rescore probabilities.")
    parser.add_argument("--tta-passes", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--min-pos-per-stratum", type=int, default=3)
    parser.add_argument("--min-neg-per-stratum", type=int, default=3)
    parser.add_argument("--predictions-csv", default=None)
    parser.add_argument("--force-rescore", action="store_true")
    parser.add_argument("--match-year", action="store_true")
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP_N)
    parser.add_argument("--permutation-n", type=int, default=DEFAULT_PERMUTATION_N)
    parser.add_argument("--stat-seed", type=int, default=DEFAULT_STAT_SEED)
    parser.add_argument("--sensitivity-thresholds", default="2,3,5")
    parser.add_argument("--no-sensitivity", action="store_true")
    parser.add_argument("--adequacy-min-n-matched", type=int, default=100)
    parser.add_argument("--adequacy-min-retention", type=float, default=0.10)
    parser.add_argument("--adequacy-min-pairwise", type=int, default=100)
    args, unknown = parser.parse_known_args()
    if unknown:
        print(f"Ignoring notebook/kernel arguments: {' '.join(unknown)}")

    run_dir = resolve_run_dir(args.run_dir, DEFAULT_EXPERIMENT)
    output_dir = Path(args.output_dir) if args.output_dir else default_output_dir(run_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    config = read_config(run_dir)
    pair_profile = args.pair_profile or config.get("pair_profile", "ecoli_mechanism6")
    mega = load_mega_module(args.mega_model_path)
    mega.init_config(pair_profile)
    data_root = resolve_data_root(args.data_root)
    active_pairs = [tuple(row) for row in config.get("active_pairs", [])]
    if not active_pairs:
        active_pairs = [(i, org, drug) for i, (org, drug) in enumerate(mega.ORGANISM_DRUG_PAIRS)]

    panel_drugs_by_organism = {}
    for _, organism, drug in active_pairs:
        panel_drugs_by_organism.setdefault(organism, []).append(drug)

    output_names = output_names_for_score_family(args.score_family)
    default_predictions_csv = output_dir / output_names["predictions"]
    predictions_csv = Path(args.predictions_csv) if args.predictions_csv else default_predictions_csv
    if predictions_csv.exists() and not args.force_rescore:
        print(f"Reusing scored predictions: {predictions_csv}")
        pred_df = pd.read_csv(predictions_csv)
    else:
        if args.row_template_csv:
            row_template = Path(args.row_template_csv)
            print(f"Reusing row/background template: {row_template}")
            template_df = pd.read_csv(row_template)
            if "prob" in template_df.columns:
                template_df = template_df.drop(columns=["prob"])
            all_rows = template_df.to_dict("records")
        else:
            eval_sites = [s.strip() for s in args.eval_sites.split(",") if s.strip()]
            all_rows = []
            for site_label in eval_sites:
                site_rows = build_prediction_rows(
                    mega, data_root, site_label, active_pairs,
                    panel_drugs_by_organism, mega.TRAIN_SITE, mega.TEST_YEAR,
                    match_year=args.match_year,
                )
                print(f"{site_label}: built {len(site_rows)} focal prediction rows")
                all_rows.extend(site_rows)
        if not all_rows:
            raise RuntimeError("No prediction rows were built. Check DRIAMS paths and labels.")

        if args.score_family == "cnn":
            models = load_checkpoint_models(mega, run_dir, config, ckpt_dir=args.ckpt_dir)
            scored_rows = score_rows(
                mega, all_rows, models,
                temperature=float(config.get("temperature", 1.0)),
                tta_passes=args.tta_passes,
                batch_size=args.batch_size,
            )
        else:
            models = load_lgbm_models(
                run_dir, model_dir=args.lgbm_model_dir,
                score_family=args.score_family)
            scored_rows = score_lgbm_rows(
                mega, all_rows, models,
                score_family=args.score_family,
                batch_size=args.batch_size,
            )
        pred_df = pd.DataFrame(scored_rows)
        pred_df.to_csv(default_predictions_csv, index=False)

    summary_df, retained_df = compute_background_matched_summary(
        pred_df,
        min_pos_per_stratum=args.min_pos_per_stratum,
        min_neg_per_stratum=args.min_neg_per_stratum,
        match_year=args.match_year,
        bootstrap_n=args.bootstrap_n,
        permutation_n=args.permutation_n,
        random_seed=args.stat_seed,
        adequacy_min_n_matched=args.adequacy_min_n_matched,
        adequacy_min_retention=args.adequacy_min_retention,
        adequacy_min_pairwise=args.adequacy_min_pairwise,
    )
    summary_df.to_csv(output_dir / output_names["summary"], index=False)
    retained_df.to_csv(output_dir / output_names["retained"], index=False)
    write_markdown(summary_df, output_dir / output_names["summary_md"])

    wrote = [
        output_names["predictions"],
        output_names["retained"],
        output_names["summary"],
        output_names["summary_md"],
    ]
    if not args.no_sensitivity:
        thresholds = parse_sensitivity_thresholds(args.sensitivity_thresholds)
        sensitivity_df = run_sensitivity_summaries(
            pred_df,
            thresholds,
            match_year=args.match_year,
            random_seed=args.stat_seed,
            adequacy_min_n_matched=args.adequacy_min_n_matched,
            adequacy_min_retention=args.adequacy_min_retention,
            adequacy_min_pairwise=args.adequacy_min_pairwise,
        )
        sensitivity_df.to_csv(output_dir / output_names["sensitivity"], index=False)
        write_markdown(sensitivity_df, output_dir / output_names["sensitivity_md"])
        wrote.extend([
            output_names["sensitivity"],
            output_names["sensitivity_md"],
        ])

    print("Wrote:")
    for name in wrote:
        print(f"  {output_dir / name}")


if __name__ == "__main__":
    main()
