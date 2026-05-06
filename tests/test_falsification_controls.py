import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import falsification_controls as falsify


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "falsification_controls.py"


class FalsificationControlsTests(unittest.TestCase):
    def rows(self):
        return [
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.1", "background_signature": "B=0", "background_resistant_count": "0"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.9", "background_signature": "B=0", "background_resistant_count": "0"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "0", "prob": "0.2", "background_signature": "B=1", "background_resistant_count": "1"},
            {"site": "A-2018", "organism": "Escherichia coli", "drug": "Ciprofloxacin", "label": "1", "prob": "0.8", "background_signature": "B=1", "background_resistant_count": "1"},
        ]

    def test_controls_report_score_exceeds_background_only_and_shuffled_null(self):
        rows = falsify.build_falsification_rows(self.rows(), n_permutations=50, min_n=4, min_pos=2, min_neg=2, seed=7)

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertAlmostEqual(row["observed_auc"], 1.0)
        self.assertAlmostEqual(row["background_burden_auc"], 0.5)
        self.assertGreater(row["observed_minus_burden_auc"], 0.4)
        self.assertGreater(row["observed_minus_shuffle_null_auc"], 0.1)
        self.assertEqual(row["control_interpretation"], "focal_score_exceeds_controls")

    def test_cli_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "predictions.csv"
            output_dir = tmp / "falsification"
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
                    "--permutations",
                    "10",
                    "--min-n",
                    "4",
                    "--min-pos",
                    "2",
                    "--min-neg",
                    "2",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "falsification_controls.csv").exists())
            self.assertTrue((output_dir / "falsification_controls.md").exists())


if __name__ == "__main__":
    unittest.main()
