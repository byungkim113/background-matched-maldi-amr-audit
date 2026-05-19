#!/usr/bin/env python3
"""Run the missing LGBM audits and rebuild the model-class matrix.

This runner is intentionally Kaggle-friendly. It can be run as a normal script
from the repository, or pasted into a Kaggle notebook cell. When paths are not
given, it tries to find:

* this repository under /kaggle/working, cloning it if needed;
* the DRIAMS data root by locating a DRIAMS-A directory;
* the E. coli and S. aureus run folders by their experiment names.

Use --dry-run first. It prints the resolved paths and all commands without
running the expensive exports/audits.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


REPO_NAME = "background-matched-maldi-amr-audit"
REPO_URL = "https://github.com/byungkim113/background-matched-maldi-amr-audit.git"
REPO_BRANCH = "weis-lr-results"
ECOLI_RUN_NAME = "exp_ecoli_mechanism6_drugid_mae30"
SAUREUS_RUN_NAME = "exp_saureus_panel_oxa_background_mae30"
DEFAULT_BOOTSTRAP_N = 500
DEFAULT_PERMUTATION_N = 500
DEFAULT_TRAIN_IF_MISSING = True
REPO_MARKERS = (
    "run_background_audit_framework.py",
    "scripts/export_lgbm_predictions_for_audit.py",
    "scripts/build_model_class_matrix.py",
)
SAUREUS_PANEL_PROFILE = (
    ("Staphylococcus aureus", "Oxacillin"),
    ("Staphylococcus aureus", "Penicillin"),
    ("Staphylococcus aureus", "Ciprofloxacin"),
    ("Staphylococcus aureus", "Erythromycin"),
    ("Staphylococcus aureus", "Clindamycin"),
    ("Staphylococcus aureus", "Gentamicin"),
    ("Staphylococcus aureus", "Fusidic acid"),
)


@dataclass(frozen=True)
class PipelineConfig:
    repo_root: Path
    data_root: Path
    ecoli_run_dir: Path
    saureus_run_dir: Path
    output_root: Path
    mega_model: Path
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


def unique_paths(paths: Sequence[Path]) -> list[Path]:
    seen: set[str] = set()
    unique = []
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def default_output_root() -> Path:
    kaggle_working = Path("/kaggle/working")
    if kaggle_working.exists():
        return kaggle_working / "model_class_matrix_outputs"
    return Path("outputs/analysis_outputs")


def is_repo_root(path: Path) -> bool:
    return all((path / marker).exists() for marker in REPO_MARKERS)


def default_repo_root() -> Path:
    file_name = globals().get("__file__")
    candidates = []
    if file_name:
        candidates.append(Path(file_name).resolve().parents[1])
    cwd = Path.cwd().resolve()
    candidates.extend([cwd, cwd / REPO_NAME, Path("/kaggle/working") / REPO_NAME])
    for candidate in unique_paths(candidates):
        if is_repo_root(candidate):
            return candidate.resolve()
    return cwd


def default_search_roots(extra: Sequence[Path] | None = None) -> list[Path]:
    roots: list[Path] = []
    if extra is not None:
        roots.extend(extra)
    else:
        env_roots = os.environ.get("MALDI_AUDIT_SEARCH_ROOTS")
        if env_roots:
            roots.extend(Path(item) for item in env_roots.split(os.pathsep) if item)
        cwd = Path.cwd().resolve()
        # Limit parent traversal to avoid including filesystem roots like / or C:\
        roots.extend([cwd, *list(cwd.parents)[:3], Path("/kaggle/working"), Path("/kaggle/input")])
    return [path for path in unique_paths(roots) if path.exists()]


def shallow_child_dirs(path: Path) -> list[Path]:
    try:
        return [child for child in path.iterdir() if child.is_dir()]
    except OSError:
        return []


def repo_candidates(
    explicit: Path | None,
    search_roots: Sequence[Path],
    *,
    include_runtime_candidates: bool = True,
) -> list[Path]:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.extend([explicit, explicit / REPO_NAME])
    if include_runtime_candidates:
        file_name = globals().get("__file__")
        if file_name:
            candidates.append(Path(file_name).resolve().parents[1])
        cwd = Path.cwd().resolve()
        candidates.extend([cwd, cwd / REPO_NAME, Path("/kaggle/working") / REPO_NAME])
    for root in search_roots:
        candidates.extend([root, root / REPO_NAME])
        candidates.extend(shallow_child_dirs(root))
    return unique_paths(candidates)


def clone_destination(explicit: Path | None) -> Path | None:
    if explicit is not None:
        return explicit if explicit.name == REPO_NAME else explicit / REPO_NAME
    kaggle_working = Path("/kaggle/working")
    if kaggle_working.exists():
        return kaggle_working / REPO_NAME
    return None


def clone_repo(destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    print(f"Repository not found. Cloning {REPO_BRANCH} into {destination}...", flush=True)
    subprocess.run(
        ["git", "clone", "-b", REPO_BRANCH, REPO_URL, str(destination)],
        check=True,
    )


def discover_repo_root(
    *,
    explicit: Path | None,
    search_roots: Sequence[Path] | None = None,
    auto_clone: bool = True,
) -> Path:
    include_runtime = explicit is None and search_roots is None
    roots = default_search_roots(search_roots)
    for candidate in repo_candidates(explicit, roots, include_runtime_candidates=include_runtime):
        if is_repo_root(candidate):
            return candidate.resolve()

    destination = clone_destination(explicit)
    if auto_clone and destination is not None:
        clone_repo(destination)
        if is_repo_root(destination):
            return destination.resolve()

    searched = "\n".join(
        str(path) for path in repo_candidates(explicit, roots, include_runtime_candidates=include_runtime)[:30]
    )
    raise FileNotFoundError(
        "Could not find the repository root with required files. "
        "Run from the cloned repo or allow auto-clone.\nSearched:\n" + searched
    )


def discover_driams_root(search_roots: Sequence[Path], explicit: Path | None = None) -> Path:
    if explicit is not None:
        return explicit.resolve()

    candidates = [Path("/kaggle/input/datasets/drscarlat/driams")]
    for root in default_search_roots(search_roots):
        candidates.extend([root, root / "DRIAMS", root / "driams"])
    for candidate in unique_paths(candidates):
        if (candidate / "DRIAMS-A").exists():
            return candidate.resolve()

    # Bounded 3-level search — avoids scanning filesystem roots via rglob
    for root in default_search_roots(search_roots):
        for pattern in ("DRIAMS-A", "*/DRIAMS-A", "*/*/DRIAMS-A"):
            for site_dir in sorted(root.glob(pattern)):
                if site_dir.is_dir():
                    return site_dir.parent.resolve()

    raise FileNotFoundError("Could not discover DRIAMS data root. Pass --data-root explicitly.")


def discover_run_dir(run_name: str, search_roots: Sequence[Path], explicit: Path | None = None) -> Path:
    if explicit is not None:
        return explicit.resolve()

    default_working = Path("/kaggle/working/runs") / run_name
    if default_working.exists():
        return default_working.resolve()

    for root in default_search_roots(search_roots):
        direct_candidates = [root / run_name, root / "runs" / run_name]
        for candidate in direct_candidates:
            if candidate.exists():
                return candidate.resolve()
        for candidate in sorted(root.rglob(run_name)):
            if candidate.is_dir():
                return candidate.resolve()

    raise FileNotFoundError(f"Could not discover run directory {run_name}. Pass the matching flag explicitly.")


def saureus_panel_compat_block() -> str:
    rows = ",\n    ".join(repr(pair) for pair in SAUREUS_PANEL_PROFILE)
    return (
        "\n# Added by run_model_class_matrix_pipeline.py for Sa/Oxa audit compatibility.\n"
        f"PAIR_PROFILES.setdefault(\"saureus_panel\", [\n    {rows},\n])\n\n"
    )


def ensure_mega_model_compat(repo_root: Path, output_root: Path) -> Path:
    source = repo_root / "Mega_Model.py"
    if not source.exists():
        raise FileNotFoundError(f"Mega_Model.py not found: {source}")

    text = source.read_text()
    if '"saureus_panel"' in text or "'saureus_panel'" in text:
        return source.resolve()

    marker = 'RUN14_OVERLAP_PAIRS = list(PAIR_PROFILES["run14"])'
    if marker not in text:
        raise ValueError(f"Could not patch Sa/Oxa profile into {source}; marker not found.")

    compat_dir = output_root / "_compat"
    compat_dir.mkdir(parents=True, exist_ok=True)
    compat_path = compat_dir / "Mega_Model_with_saureus_panel.py"
    compat_text = text.replace(marker, saureus_panel_compat_block() + marker, 1)
    compat_path.write_text(compat_text)
    print(f"Using compatibility Mega_Model with saureus_panel profile: {compat_path}", flush=True)
    return compat_path.resolve()


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
        "--mega-model",
        str(config.mega_model),
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
    if not config.ecoli_run_dir.exists():
        raise FileNotFoundError(f"E. coli run directory does not exist: {config.ecoli_run_dir}")
    if not config.saureus_run_dir.exists():
        raise FileNotFoundError(f"S. aureus run directory does not exist: {config.saureus_run_dir}")
    if not config.mega_model.exists():
        raise FileNotFoundError(f"Mega_Model file does not exist: {config.mega_model}")


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


def split_notebook_kernel_args(args: Sequence[str]) -> tuple[list[str], list[str]]:
    ignored: list[str] = []
    remaining: list[str] = []
    index = 0
    while index < len(args):
        arg = args[index]
        if arg == "-f" and index + 1 < len(args):
            ignored.extend([arg, args[index + 1]])
            index += 2
            continue
        if arg.startswith("--HistoryManager.") or arg.startswith("--IPKernelApp."):
            ignored.append(arg)
            index += 1
            continue
        remaining.append(arg)
        index += 1
    return ignored, remaining


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run LGBM exports/audits and rebuild model-class matrix.")
    parser.add_argument("--repo-root", type=Path, default=None, help="Repository root. Auto-discovered if omitted.")
    parser.add_argument("--data-root", type=Path, default=None, help="DRIAMS root containing DRIAMS-A/B/C/D.")
    parser.add_argument("--ecoli-run-dir", type=Path, default=None, help=f"Run directory for {ECOLI_RUN_NAME}.")
    parser.add_argument("--saureus-run-dir", type=Path, default=None, help=f"Run directory for {SAUREUS_RUN_NAME}.")
    parser.add_argument("--output-root", type=Path, default=default_output_root())
    parser.add_argument("--bootstrap-n", type=int, default=DEFAULT_BOOTSTRAP_N)
    parser.add_argument("--permutation-n", type=int, default=DEFAULT_PERMUTATION_N)
    parser.add_argument("--no-train-if-missing", action="store_true")
    parser.add_argument("--no-auto-clone", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parsed, unknown = parser.parse_known_args(argv)
    ignored, remaining = split_notebook_kernel_args(unknown)
    if ignored:
        print(f"Ignoring notebook/kernel arguments: {shlex.join(ignored)}", file=sys.stderr, flush=True)
    if remaining:
        parser.error(f"unrecognized arguments: {shlex.join(remaining)}")
    return parsed


def resolve_config(args: argparse.Namespace, search_roots: Sequence[Path] | None = None) -> PipelineConfig:
    roots = default_search_roots(search_roots)
    repo_root = discover_repo_root(
        explicit=args.repo_root,
        search_roots=roots,
        auto_clone=not args.no_auto_clone,
    )
    roots = default_search_roots([repo_root, repo_root.parent, *roots])
    data_root = discover_driams_root(roots, explicit=args.data_root)
    ecoli_run_dir = discover_run_dir(ECOLI_RUN_NAME, roots, explicit=args.ecoli_run_dir)
    saureus_run_dir = discover_run_dir(SAUREUS_RUN_NAME, roots, explicit=args.saureus_run_dir)
    resolved_output_root = under_repo(repo_root, args.output_root)
    mega_model = ensure_mega_model_compat(repo_root, resolved_output_root)
    return PipelineConfig(
        repo_root=repo_root,
        data_root=data_root,
        ecoli_run_dir=ecoli_run_dir,
        saureus_run_dir=saureus_run_dir,
        output_root=args.output_root,
        mega_model=mega_model,
        bootstrap_n=args.bootstrap_n,
        permutation_n=args.permutation_n,
        train_if_missing=not args.no_train_if_missing,
    )


def print_config(config: PipelineConfig, *, dry_run: bool) -> None:
    mode = "dry run" if dry_run else "real run"
    print(f"Resolved pipeline paths ({mode}):", flush=True)
    print(f"  repo_root       = {config.repo_root}", flush=True)
    print(f"  data_root       = {config.data_root}", flush=True)
    print(f"  ecoli_run_dir   = {config.ecoli_run_dir}", flush=True)
    print(f"  saureus_run_dir = {config.saureus_run_dir}", flush=True)
    print(f"  mega_model      = {config.mega_model}", flush=True)
    print(f"  output_root     = {under_repo(config.repo_root, config.output_root)}", flush=True)
    print(f"  bootstrap_n     = {config.bootstrap_n}", flush=True)
    print(f"  permutation_n   = {config.permutation_n}", flush=True)
    print(f"  train_missing   = {config.train_if_missing}", flush=True)


def main() -> None:
    args = parse_args()
    config = resolve_config(args)
    print_config(config, dry_run=args.dry_run)
    steps = build_pipeline_steps(config)
    if not args.dry_run:
        validate_config(config)
    run_steps(config, steps, dry_run=args.dry_run, skip_existing=args.skip_existing)
    print("\nPipeline finished.", flush=True)


if __name__ == "__main__":
    main()
