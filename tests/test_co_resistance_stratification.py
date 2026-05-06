import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest

from scripts import co_resistance_stratification as strat


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "co_resistance_stratification.py"


class CoResistanceStratificationTests(unittest.TestCase):
    def rows(self):
        return [
            {
                "site": "A-2018",
                "year": "2018",
                "uid": "i1",
                "organism": "Escherichia coli",
                "drug": "Ciprofloxacin",
                "label": "0",
                "prob": "0.10",
                "background_signature": "Norfloxacin=S|Ceftriaxone=S|Cefepime=S",
            },
            {
                "site": "A-2018",
                "year": "2018",
                "uid": "i2",
                "organism": "Escherichia coli",
                "drug": "Ciprofloxacin",
                "label": "1",
                "prob": "0.90",
                "background_signature": "Norfloxacin=S|Ceftriaxone=S|Cefepime=S",
            },
            {
                "site": "A-2018",
                "year": "2018",
                "uid": "i3",
                "organism": "Escherichia coli",
                "drug": "Ciprofloxacin",
                "label": "0",
                "prob": "0.80",
                "background_signature": "Norfloxacin=R|Ceftriaxone=R|Cefepime=S",
            },
            {
                "site": "A-2018",
                "year": "2018",
                "uid": "i4",
                "organism": "Escherichia coli",
                "drug": "Ciprofloxacin",
                "label": "1",
                "prob": "0.20",
                "background_signature": "Norfloxacin=R|Ceftriaxone=R|Cefepime=S",
            },
        ]

    def test_build_stratification_rows_reports_burden_and_exact_background_auc(self):
        rows = strat.build_stratification_rows(self.rows(), min_n=2, min_pos=1, min_neg=1)

        burden_rows = [r for r in rows if r["stratum_type"] == "background_resistance_burden"]
        exact_rows = [r for r in rows if r["stratum_type"] == "exact_background_signature"]

        burden_zero = next(r for r in burden_rows if r["stratum_value"] == "0")
        burden_two = next(r for r in burden_rows if r["stratum_value"] == "2")

        self.assertEqual(burden_zero["n"], 2)
        self.assertEqual(burden_zero["n_r"], 1)
        self.assertEqual(burden_zero["auc"], 1.0)
        self.assertEqual(burden_two["auc"], 0.0)
        self.assertEqual(len(exact_rows), 2)

    def test_cli_writes_csv_and_markdown_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "predictions.csv"
            output_dir = tmp / "stratification"
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
                    "2",
                    "--min-pos",
                    "1",
                    "--min-neg",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "co_resistance_stratification.csv").exists())
            self.assertTrue((output_dir / "co_resistance_stratification.md").exists())


if __name__ == "__main__":
    unittest.main()
