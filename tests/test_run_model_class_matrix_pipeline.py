import pathlib
import unittest

from scripts import run_model_class_matrix_pipeline as pipeline


class RunModelClassMatrixPipelineTests(unittest.TestCase):
    def test_builds_expected_lgbm_export_audit_and_matrix_steps(self):
        config = pipeline.PipelineConfig(
            repo_root=pathlib.Path("/repo"),
            data_root=pathlib.Path("/data/driams"),
            ecoli_run_dir=pathlib.Path("/runs/ecoli"),
            saureus_run_dir=pathlib.Path("/runs/saureus"),
            output_root=pathlib.Path("/repo/outputs/analysis_outputs"),
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
