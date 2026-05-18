import math
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from scripts import marisma_end_to_end_kaggle as marisma


class MarismaEndToEndKaggleTests(unittest.TestCase):
    def test_bruker_axis_matches_validated_marisma_example(self):
        params = {
            "DELAY": 19626.0,
            "DW": 2.0,
            "ML1": 5419123.44049644,
            "ML2": 437.772115091131,
            "ML3": -0.012041034464694,
            "TD": 21409.0,
        }

        mz = marisma.bruker_mz_axis(params)

        self.assertEqual(mz.shape, (21409,))
        self.assertAlmostEqual(float(mz[0]), 2000.2688, places=3)
        self.assertAlmostEqual(float(mz[-1]), 21004.2601, places=3)
        self.assertTrue(np.all(np.diff(mz) > 0))

    def test_vectorize_spectrum_returns_finite_6000_bin_vector(self):
        mz = np.linspace(2000.0, 21000.0, 21409)
        intensity = np.zeros_like(mz, dtype=np.float32)
        intensity[1000:1010] = 100
        intensity[9000:9020] = 250

        vector = marisma.vectorize_spectrum(mz, intensity)

        self.assertEqual(vector.shape, (6000,))
        self.assertTrue(np.isfinite(vector).all())
        self.assertAlmostEqual(float(vector.mean()), 0.0, places=5)
        self.assertAlmostEqual(float(vector.std()), 1.0, places=5)

    def test_resolve_spot_path_strips_marisma_prefix(self):
        root = Path("/tmp/kaggle/input/marisma/MARISMa")
        resolved = marisma.resolve_spot_path(
            root,
            "/MARISMa/2024/Escherichia/Coli/b23ff733/0_A10",
        )

        self.assertEqual(
            resolved,
            root / "2024" / "Escherichia" / "Coli" / "b23ff733" / "0_A10",
        )

    def test_read_spot_uses_little_endian_int32_and_metadata_axis(self):
        params = {
            "DELAY": 19626.0,
            "DW": 2.0,
            "ML1": 5419123.44049644,
            "ML2": 437.772115091131,
            "ML3": -0.012041034464694,
            "TD": 8.0,
        }
        acqu_text = "\n".join(f"##${key}= {value}" for key, value in params.items())

        with tempfile.TemporaryDirectory() as tmp:
            spot = Path(tmp) / "b23ff733" / "0_A10"
            analysis = spot / "1" / "1SLin"
            processed = analysis / "pdata" / "1"
            processed.mkdir(parents=True)
            (analysis / "acqu").write_text(acqu_text)
            np.arange(8, dtype="<i4").tofile(processed / "1r")

            mz, intensity = marisma.read_bruker_spot(spot)

        self.assertEqual(mz.shape, (8,))
        self.assertEqual(intensity.tolist(), list(range(8)))
        self.assertTrue(math.isclose(float(mz[0]), 2000.2688, rel_tol=0, abs_tol=1e-3))

    def test_aggregate_predictions_to_isolate_drug_averages_spot_rows_and_reports_duplicates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prediction_csv = root / "predictions.csv"
            output_csv = root / "isolate_level.csv"
            report_json = root / "report.json"
            prediction_csv.write_text(
                "\n".join(
                    [
                        "model_name,site,year,isolate_id,spot_id,organism,drug,label,prob",
                        "Mega,MARISMa,2024,iso1,spot1,Escherichia coli,Ciprofloxacin,1,0.8",
                        "Mega,MARISMa,2024,iso1,spot2,Escherichia coli,Ciprofloxacin,1,0.6",
                        "Mega,MARISMa,2024,iso1,spot2,Escherichia coli,Ciprofloxacin,1,0.6",
                        "Mega,MARISMa,2024,iso1,spot1,Escherichia coli,Cefepime,0,0.2",
                    ]
                )
                + "\n"
            )

            aggregated = marisma.aggregate_predictions_to_isolate_drug(
                prediction_csv,
                output_csv,
                report_json,
            )

            self.assertEqual(len(aggregated), 2)
            cipro = aggregated[aggregated["drug"].eq("Ciprofloxacin")].iloc[0]
            self.assertAlmostEqual(float(cipro["prob"]), (0.8 + 0.6 + 0.6) / 3.0)
            self.assertEqual(int(cipro["n_prediction_rows"]), 3)
            self.assertEqual(int(cipro["n_unique_spots"]), 2)

            report = json.loads(report_json.read_text())
            self.assertEqual(report["input_rows"], 4)
            self.assertEqual(report["isolate_drug_rows"], 2)
            self.assertEqual(report["duplicate_site_year_isolate_drug_rows"], 2)
            self.assertEqual(report["exact_duplicate_rows"], 1)

    def test_aggregate_predictions_excludes_conflicting_label_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prediction_csv = root / "predictions.csv"
            output_csv = root / "isolate_level.csv"
            report_json = root / "report.json"
            prediction_csv.write_text(
                "\n".join(
                    [
                        "site,year,isolate_id,spot_id,organism,drug,label,prob",
                        "MARISMa,2024,iso1,spot1,Escherichia coli,Cefepime,0,0.8",
                        "MARISMa,2024,iso1,spot1,Escherichia coli,Cefepime,1,0.8",
                        "MARISMa,2024,iso2,spot1,Escherichia coli,Cefepime,1,0.7",
                    ]
                )
                + "\n"
            )

            aggregated = marisma.aggregate_predictions_to_isolate_drug(
                prediction_csv,
                output_csv,
                report_json,
            )

            self.assertEqual(len(aggregated), 1)
            self.assertEqual(aggregated.iloc[0]["isolate_id"], "iso2")
            report = json.loads(report_json.read_text())
            self.assertEqual(report["conflicting_isolate_drug_groups_excluded"], 1)
            self.assertTrue(output_csv.with_name("marisma_conflicting_isolate_drug_rows.csv").exists())


if __name__ == "__main__":
    unittest.main()
