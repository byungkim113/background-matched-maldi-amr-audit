import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import deployment_decision_framework as deploy


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "deployment_decision_framework.py"


class DeploymentDecisionFrameworkTests(unittest.TestCase):
    def test_decision_categories_cover_matched_high_matched_low_and_underpowered(self):
        self.assertEqual(
            deploy.decision_category(raw_auc=0.82, centered_auc=0.70, adequacy="interpretable", ece=0.04),
            "candidate_for_controlled_deployment",
        )
        self.assertEqual(
            deploy.decision_category(raw_auc=0.82, centered_auc=0.53, adequacy="interpretable", ece=0.04),
            "background_dependent_retrain_locally",
        )
        self.assertEqual(
            deploy.decision_category(raw_auc=0.82, centered_auc=None, adequacy="caution_low_n_matched", ece=0.04),
            "insufficient_matched_evidence",
        )
        self.assertEqual(
            deploy.decision_category(raw_auc=0.82, centered_auc=0.70, adequacy="interpretable", ece=0.16),
            "ranking_only_recalibrate_before_clinical_use",
        )

    def test_build_readiness_rows_joins_audit_and_calibration_metrics(self):
        audit_rows = [
            {
                "pair": "E. coli / Cipro",
                "site": "A-2018",
                "raw_auc_95ci": "0.823 (0.799-0.849)",
                "stratum_centered_auc_95ci": "0.703 (0.664-0.747)",
                "matched_retention_pct": "61.9",
                "n_matched": "836",
                "valid_strata": "9",
                "adequacy": "interpretable",
                "interpretation": "Survives background matching",
            },
            {
                "pair": "E. coli / Amox-Clav",
                "site": "DRIAMS-C",
                "raw_auc_95ci": "0.535 (0.494-0.577)",
                "stratum_centered_auc_95ci": "0.497 (0.449-0.547)",
                "matched_retention_pct": "89.2",
                "n_matched": "805",
                "valid_strata": "10",
                "adequacy": "interpretable",
                "interpretation": "Collapses to chance after matching",
            },
        ]
        calibration_rows = [
            {"site": "A-2018", "drug": "Ciprofloxacin", "expected_calibration_error": "0.04", "brier": "0.12", "calibration_label": "well_calibrated"},
            {"site": "DRIAMS-C", "drug": "Amoxicillin-Clavulanic acid", "expected_calibration_error": "0.14", "brier": "0.28", "calibration_label": "poorly_calibrated"},
        ]

        rows = deploy.build_readiness_rows(audit_rows, calibration_rows)

        self.assertEqual(rows[0]["decision_category"], "candidate_for_controlled_deployment")
        self.assertEqual(rows[0]["recommended_action"], "Proceed only with local calibration check, locked thresholds, and ongoing drift monitoring.")
        self.assertEqual(rows[1]["decision_category"], "not_deployment_ready")

    def test_cli_writes_rule_and_readiness_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            audit_csv = tmp / "audit.csv"
            calibration_csv = tmp / "calibration.csv"
            output_dir = tmp / "deployment"
            audit_rows = [
                {
                    "pair": "E. coli / Cipro",
                    "site": "A-2018",
                    "raw_auc_95ci": "0.823 (0.799-0.849)",
                    "stratum_centered_auc_95ci": "0.703 (0.664-0.747)",
                    "matched_retention_pct": "61.9",
                    "n_matched": "836",
                    "valid_strata": "9",
                    "adequacy": "interpretable",
                    "interpretation": "Survives background matching",
                }
            ]
            calibration_rows = [
                {"site": "A-2018", "drug": "Ciprofloxacin", "expected_calibration_error": "0.04", "brier": "0.12", "calibration_label": "well_calibrated"}
            ]
            with audit_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()))
                writer.writeheader()
                writer.writerows(audit_rows)
            with calibration_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(calibration_rows[0].keys()))
                writer.writeheader()
                writer.writerows(calibration_rows)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--audit-csv",
                    str(audit_csv),
                    "--calibration-csv",
                    str(calibration_csv),
                    "--output-dir",
                    str(output_dir),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "deployment_decision_rules.csv").exists())
            self.assertTrue((output_dir / "deployment_readiness_by_pair.csv").exists())
            self.assertTrue((output_dir / "deployment_decision_framework.md").exists())


if __name__ == "__main__":
    unittest.main()
