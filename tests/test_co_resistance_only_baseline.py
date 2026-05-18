import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import co_resistance_only_baseline as baseline


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "co_resistance_only_baseline.py"


class CoResistanceOnlyBaselineTests(unittest.TestCase):
    def rows(self):
        rows = []
        for idx, (background, labels) in enumerate(
            [
                ("A=R|B=S", [1, 1, 1, 0]),
                ("A=S|B=S", [1, 0, 0, 0]),
            ],
            start=1,
        ):
            for rep, label in enumerate(labels, start=1):
                rows.append(
                    {
                        "site": "A-2018",
                        "organism": "Escherichia coli",
                        "drug": "Ciprofloxacin",
                        "label": str(label),
                        "prob": str(0.9 if label == 1 else 0.1),
                        "background_signature": background,
                        "background_resistant_count": str(background.count("=R")),
                        "uid": f"i{idx}_{rep}",
                    }
                )
        return rows

    def test_exact_background_leave_one_out_baseline_uses_non_focal_pattern(self):
        rows = baseline.build_baseline_rows(self.rows(), min_n=8, min_pos=2, min_neg=2, smoothing_alpha=1.0)

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertAlmostEqual(row["observed_auc"], 1.0)
        self.assertAlmostEqual(row["background_burden_auc"], 0.75)
        self.assertAlmostEqual(row["exact_background_auc"], 0.5625)
        self.assertAlmostEqual(row["observed_minus_exact_background_auc"], 0.4375)
        self.assertEqual(row["adequacy_label"], "interpretable")

    def test_cli_writes_csv_and_markdown_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "predictions.csv"
            output_dir = tmp / "baseline"
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
                    "--min-n",
                    "8",
                    "--min-pos",
                    "2",
                    "--min-neg",
                    "2",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "co_resistance_only_baseline.csv").exists())
            self.assertTrue((output_dir / "co_resistance_only_baseline.md").exists())


if __name__ == "__main__":
    unittest.main()
