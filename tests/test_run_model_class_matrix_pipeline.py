import pathlib
import tempfile
from types import SimpleNamespace
import unittest

from scripts import run_model_class_matrix_pipeline as pipeline


def make_repo_root(path: pathlib.Path) -> pathlib.Path:
    (path / "scripts").mkdir(parents=True)
    (path / "run_background_audit_framework.py").write_text("")
    (path / "scripts" / "export_lgbm_predictions_for_audit.py").write_text("")
    (path / "scripts" / "build_model_class_matrix.py").write_text("")
    (path / "Mega_Model.py").write_text(
        'PAIR_PROFILES = {"run14": []}\nRUN14_OVERLAP_PAIRS = list(PAIR_PROFILES["run14"])\n'
    )
    return path


class RunModelClassMatrixPipelineTests(unittest.TestCase):
    def test_default_repo_root_falls_back_to_cwd_without_file_global(self):
        old_file = pipeline.__dict__.pop("__file__", None)
        try:
            self.assertEqual(pipeline.default_repo_root(), pathlib.Path.cwd().resolve())
        finally:
            if old_file is not None:
                pipeline.__dict__["__file__"] = old_file

    def test_parse_args_ignores_notebook_kernel_arguments(self):
        args = pipeline.parse_args(
            [
                "--dry-run",
                "--repo-root",
                "/repo",
                "-f",
                "/tmp/tmp31x05qkj.json",
                "--HistoryManager.hist_file=:memory:",
            ]
        )

        self.assertTrue(args.dry_run)
        self.assertEqual(args.repo_root, pathlib.Path("/repo"))

    def test_discovers_repo_root_nested_under_kaggle_working(self):
        with tempfile.TemporaryDirectory() as tmp:
            working = pathlib.Path(tmp) / "kaggle" / "working"
            repo = make_repo_root(working / "background-matched-maldi-amr-audit")

            self.assertEqual(
                pipeline.discover_repo_root(explicit=working, search_roots=[working], auto_clone=False),
                repo.resolve(),
            )

    def test_discovers_driams_root_from_site_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            data_root = pathlib.Path(tmp) / "input" / "driams-dataset" / "DRIAMS"
            (data_root / "DRIAMS-A").mkdir(parents=True)

            self.assertEqual(pipeline.discover_driams_root([pathlib.Path(tmp) / "input"]), data_root.resolve())

    def test_resolves_config_from_discovered_kaggle_like_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            repo = make_repo_root(root / "working" / "background-matched-maldi-amr-audit")
            data_root = root / "input" / "driams-dataset" / "DRIAMS"
            (data_root / "DRIAMS-A").mkdir(parents=True)
            ecoli_run = root / "input" / "NewRuns" / "runs" / "exp_ecoli_mechanism6_drugid_mae30"
            saureus_run = root / "input" / "SA-OxaPanelData" / "runs" / "exp_saureus_panel_oxa_background_mae30"
            ecoli_run.mkdir(parents=True)
            saureus_run.mkdir(parents=True)

            args = SimpleNamespace(
                repo_root=None,
                data_root=None,
                ecoli_run_dir=None,
                saureus_run_dir=None,
                output_root=root / "working" / "model_class_matrix_outputs",
                bootstrap_n=7,
                permutation_n=9,
                no_train_if_missing=False,
                no_auto_clone=True,
            )

            config = pipeline.resolve_config(args, search_roots=[root / "working", root / "input"])

            self.assertEqual(config.repo_root, repo.resolve())
            self.assertEqual(config.data_root, data_root.resolve())
            self.assertEqual(config.ecoli_run_dir, ecoli_run.resolve())
            self.assertEqual(config.saureus_run_dir, saureus_run.resolve())
            self.assertTrue(config.mega_model.exists())
            self.assertIn("saureus_panel", config.mega_model.read_text())
            self.assertEqual(config.bootstrap_n, 7)
            self.assertEqual(config.permutation_n, 9)

    def test_lgbm_export_commands_use_compatible_mega_model_path(self):
        config = pipeline.PipelineConfig(
            repo_root=pathlib.Path("/repo"),
            data_root=pathlib.Path("/data/driams"),
            ecoli_run_dir=pathlib.Path("/runs/ecoli"),
            saureus_run_dir=pathlib.Path("/runs/saureus"),
            output_root=pathlib.Path("/repo/outputs/analysis_outputs"),
            mega_model=pathlib.Path("/repo/outputs/_compat/Mega_Model_with_saureus_panel.py"),
            bootstrap_n=11,
            permutation_n=13,
            train_if_missing=True,
        )

        command_text = "\n".join(" ".join(step.command) for step in pipeline.build_pipeline_steps(config))

        expected = f"--mega-model {pathlib.Path('/repo/outputs/_compat/Mega_Model_with_saureus_panel.py')}"
        self.assertIn(expected, command_text)

    def test_builds_expected_lgbm_export_audit_and_matrix_steps(self):
        config = pipeline.PipelineConfig(
            repo_root=pathlib.Path("/repo"),
            data_root=pathlib.Path("/data/driams"),
            ecoli_run_dir=pathlib.Path("/runs/ecoli"),
            saureus_run_dir=pathlib.Path("/runs/saureus"),
            output_root=pathlib.Path("/repo/outputs/analysis_outputs"),
            mega_model=pathlib.Path("/repo/Mega_Model.py"),
            bootstrap_n=11,
            permutation_n=13,
            train_if_missing=True,
        )

        steps = pipeline.build_pipeline_steps(config)
        names = [step.name for step in steps]

        self.assertEqual(
            names,
            [
                "export_ecoli_lgbm_single",
                "audit_ecoli_lgbm_single",
                "export_saureus_lgbm_single_multi",
                "audit_saureus_lgbm_multi",
                "audit_saureus_lgbm_single",
                "rebuild_model_class_matrix",
            ],
        )
        command_text = "\n".join(" ".join(step.command) for step in steps)
        self.assertIn("--pair-profile ecoli_mechanism6", command_text)
        self.assertIn("--pair-profile saureus_panel", command_text)
        self.assertIn("--variants single,multi", command_text)
        self.assertIn("--train-if-missing", command_text)
        self.assertIn("--bootstrap-n 11", command_text)
        self.assertIn("--permutation-n 13", command_text)
        self.assertIn("ecoli_lgbm_single_background_audit", command_text)
        self.assertIn("saureus_lgbm_multi_oxa_background_audit", command_text)
        self.assertIn("saureus_lgbm_single_oxa_background_audit", command_text)


if __name__ == "__main__":
    unittest.main()
