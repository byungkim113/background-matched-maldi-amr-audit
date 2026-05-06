import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import temporal_reliability_audit as temporal


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "temporal_reliability_audit.py"


class TemporalReliabilityAuditTests(unittest.TestCase):
    def test_single_period_is_flagged_as_insufficient_for_reliability_duration(self):
        rows = [
            {"site": "A-2018", "year": "2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.1", "background_resistant_count": "0"},
            {"site": "A-2018", "year": "2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.9", "background_resistant_count": "1"},
            {"site": "A-2018", "year": "2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.2", "background_resistant_count": "0"},
            {"site": "A-2018", "year": "2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.8", "background_resistant_count": "1"},
        ]

        out = temporal.build_temporal_rows(rows, min_n=4, min_pos=2, min_neg=2)

        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["reliability_status"], "insufficient_periods")
        self.assertEqual(out[0]["recommended_action"], "collect_future_periods")

    def test_later_period_with_auc_drop_or_ecology_shift_flags_recalibration_review(self):
        rows = [
            {"site": "DRIAMS-X", "year": "2020", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.1", "background_resistant_count": "0"},
            {"site": "DRIAMS-X", "year": "2020", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.2", "background_resistant_count": "0"},
            {"site": "DRIAMS-X", "year": "2020", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.8", "background_resistant_count": "1"},
            {"site": "DRIAMS-X", "year": "2020", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.9", "background_resistant_count": "1"},
            {"site": "DRIAMS-X", "year": "2021", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.9", "background_resistant_count": "3"},
            {"site": "DRIAMS-X", "year": "2021", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.8", "background_resistant_count": "3"},
            {"site": "DRIAMS-X", "year": "2021", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.2", "background_resistant_count": "4"},
            {"site": "DRIAMS-X", "year": "2021", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.1", "background_resistant_count": "4"},
        ]

        out = temporal.build_temporal_rows(
            rows,
            min_n=4,
            min_pos=2,
            min_neg=2,
            auc_drop_alert=0.2,
            burden_shift_alert=1.0,
        )
        later = next(row for row in out if row["period"] == "2021")

        self.assertEqual(later["reliability_status"], "drift_alert")
        self.assertEqual(later["recommended_action"], "recalibration_or_retraining_review")
        self.assertLessEqual(later["auc_delta_from_baseline"], -0.2)
        self.assertGreaterEqual(later["mean_background_resistant_delta"], 1.0)

    def test_cli_writes_outputs(self):
        rows = [
            {"site": "A", "year": "2020", "organism": "E", "drug": "D", "label": "0", "prob": "0.1", "background_resistant_count": "0"},
            {"site": "A", "year": "2020", "organism": "E", "drug": "D", "label": "1", "prob": "0.9", "background_resistant_count": "1"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "predictions.csv"
            output_dir = tmp / "temporal"
            with input_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                writer.writeheader()
                writer.writerows(rows)

            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--predictions-csv", str(input_csv), "--output-dir", str(output_dir), "--min-n", "2", "--min-pos", "1", "--min-neg", "1"],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "temporal_reliability.csv").exists())
            self.assertTrue((output_dir / "temporal_reliability.md").exists())


if __name__ == "__main__":
    unittest.main()
