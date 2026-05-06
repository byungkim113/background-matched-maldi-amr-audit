import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import calibration_analysis as cal


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "calibration_analysis.py"


class CalibrationAnalysisTests(unittest.TestCase):
    def rows(self):
        return [
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.10"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.20"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.80"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.90"},
        ]

    def test_grouped_calibration_reports_brier_ece_and_threshold_metrics(self):
        rows = cal.build_calibration_rows(self.rows(), n_bins=2, threshold=0.5)

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["n"], 4)
        self.assertEqual(row["n_r"], 2)
        self.assertAlmostEqual(row["auc"], 1.0)
        self.assertAlmostEqual(row["brier"], 0.025, places=6)
        self.assertAlmostEqual(row["expected_calibration_error"], 0.15, places=6)
        self.assertAlmostEqual(row["sensitivity"], 1.0)
        self.assertAlmostEqual(row["specificity"], 1.0)
        self.assertEqual(row["calibration_label"], "well_calibrated")

    def test_cli_writes_csv_and_markdown_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "predictions.csv"
            output_dir = tmp / "calibration"
            with input_csv.open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(self.rows()[0].keys()))
                writer.writeheader()
                writer.writerows(self.rows())

            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--predictions-csv",
                    str(input_csv),
                    "--output-dir",
                    str(output_dir),
                    "--n-bins",
                    "2",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "calibration_summary.csv").exists())
            self.assertTrue((output_dir / "calibration_summary.md").exists())


if __name__ == "__main__":
    unittest.main()
