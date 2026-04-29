import csv
import importlib.util
import pathlib
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
FRAMEWORK_PATH = ROOT / "run_background_audit_framework.py"


def load_framework():
    spec = importlib.util.spec_from_file_location("run_background_audit_framework", FRAMEWORK_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class BackgroundAuditFrameworkTests(unittest.TestCase):
    def write_csv(self, path, rows):
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def example_rows(self):
        rows = []
        isolates = [
            ("A1", "A", "2018", "Escherichia coli", {"Cipro": "S", "Norflox": "S", "Amox": "S"}, {"Cipro": 0.10, "Norflox": 0.10, "Amox": 0.20}),
            ("A2", "A", "2018", "Escherichia coli", {"Cipro": "R", "Norflox": "S", "Amox": "S"}, {"Cipro": 0.70, "Norflox": 0.20, "Amox": 0.30}),
            ("A3", "A", "2018", "Escherichia coli", {"Cipro": "S", "Norflox": "R", "Amox": "R"}, {"Cipro": 0.40, "Norflox": 0.80, "Amox": 0.70}),
            ("A4", "A", "2018", "Escherichia coli", {"Cipro": "R", "Norflox": "R", "Amox": "R"}, {"Cipro": 0.90, "Norflox": 0.85, "Amox": 0.60}),
            ("B1", "B", "2018", "Escherichia coli", {"Cipro": "S", "Norflox": "S", "Amox": "S"}, {"Cipro": 0.20, "Norflox": 0.10, "Amox": 0.10}),
            ("B2", "B", "2018", "Escherichia coli", {"Cipro": "R", "Norflox": "S", "Amox": "S"}, {"Cipro": 0.65, "Norflox": 0.20, "Amox": 0.20}),
            ("B3", "B", "2018", "Escherichia coli", {"Cipro": "S", "Norflox": "R", "Amox": "R"}, {"Cipro": 0.30, "Norflox": 0.75, "Amox": 0.80}),
            ("B4", "B", "2018", "Escherichia coli", {"Cipro": "R", "Norflox": "R", "Amox": "R"}, {"Cipro": 0.80, "Norflox": 0.90, "Amox": 0.70}),
        ]
        for isolate_id, site, year, organism, labels, probs in isolates:
            for drug in ["Cipro", "Norflox", "Amox"]:
                rows.append(
                    {
                        "sample_id": isolate_id,
                        "hospital": site,
                        "collection_year": year,
                        "species": organism,
                        "antibiotic": drug,
                        "phenotype": labels[drug],
                        "score": probs[drug],
                    }
                )
        return rows

    def test_framework_normalizes_custom_prediction_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "weis_like_predictions.csv"
            self.write_csv(input_csv, self.example_rows())

            fw = load_framework()
            rows = fw.read_long_predictions(
                input_csv,
                id_col="sample_id",
                site_col="hospital",
                year_col="collection_year",
                organism_col="species",
                drug_col="antibiotic",
                label_col="phenotype",
                prob_col="score",
            )

            self.assertEqual(len(rows), 24)
            self.assertEqual(rows[0]["uid"], "A1")
            self.assertEqual(rows[0]["site"], "A")
            self.assertEqual(rows[0]["drug"], "Cipro")
            self.assertEqual(rows[0]["label"], 0)
            self.assertEqual(rows[1]["label"], 0)
            self.assertIsInstance(rows[0]["prob"], float)

    def test_framework_accepts_precomputed_background_signature(self):
        fw = load_framework()
        records = [
            {
                "sample_id": "A1",
                "hospital": "A",
                "collection_year": "2018",
                "species": "Escherichia coli",
                "antibiotic": "Cipro",
                "phenotype": "S",
                "score": "0.2",
                "bg": "Norflox=S|Amox=R|Ceftriaxone=U",
            },
            {
                "sample_id": "A2",
                "hospital": "A",
                "collection_year": "2018",
                "species": "Escherichia coli",
                "antibiotic": "Cipro",
                "phenotype": "R",
                "score": "0.8",
                "bg": "Norflox=S|Amox=R|Ceftriaxone=U",
            },
        ]
        rows = fw.read_rows_from_records(
            records,
            id_col="sample_id",
            site_col="hospital",
            year_col="collection_year",
            organism_col="species",
            drug_col="antibiotic",
            label_col="phenotype",
            prob_col="score",
            background_signature_col="bg",
        )
        enriched = fw.add_background_signatures(rows)
        self.assertEqual(enriched[0]["background_signature"], "Norflox=S|Amox=R|Ceftriaxone=U")
        self.assertEqual(enriched[0]["background_known_count"], 2)
        self.assertEqual(enriched[0]["background_resistant_count"], 1)

    def test_framework_builds_background_signatures_excluding_focal_drug(self):
        fw = load_framework()
        rows = fw.read_rows_from_records(
            self.example_rows(),
            id_col="sample_id",
            site_col="hospital",
            year_col="collection_year",
            organism_col="species",
            drug_col="antibiotic",
            label_col="phenotype",
            prob_col="score",
        )
        enriched = fw.add_background_signatures(rows)
        cipro_a1 = next(r for r in enriched if r["uid"] == "A1" and r["drug"] == "Cipro")

        self.assertIn("Norflox=S", cipro_a1["background_signature"])
        self.assertIn("Amox=S", cipro_a1["background_signature"])
        self.assertNotIn("Cipro", cipro_a1["background_signature"])
        self.assertEqual(cipro_a1["background_known_count"], 2)

    def test_cli_writes_audit_and_ecology_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            input_csv = tmp / "weis_like_predictions.csv"
            output_dir = tmp / "audit"
            self.write_csv(input_csv, self.example_rows())

            result = subprocess.run(
                [
                    sys.executable,
                    str(FRAMEWORK_PATH),
                    "--predictions-csv",
                    str(input_csv),
                    "--id-col",
                    "sample_id",
                    "--site-col",
                    "hospital",
                    "--year-col",
                    "collection_year",
                    "--organism-col",
                    "species",
                    "--drug-col",
                    "antibiotic",
                    "--label-col",
                    "phenotype",
                    "--prob-col",
                    "score",
                    "--output-dir",
                    str(output_dir),
                    "--min-pos-per-stratum",
                    "1",
                    "--min-neg-per-stratum",
                    "1",
                    "--bootstrap-n",
                    "25",
                    "--permutation-n",
                    "25",
                    "--adequacy-min-n-matched",
                    "1",
                    "--adequacy-min-pairwise",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output_dir / "normalized_predictions.csv").exists())
            self.assertTrue((output_dir / "background_matched_audit_summary.csv").exists())
            self.assertTrue((output_dir / "background_matched_retained_rows.csv").exists())
            self.assertTrue((output_dir / "cross_resistance_edges.csv").exists())
            self.assertTrue((output_dir / "background_audit_with_resistance_ecology.csv").exists())
            self.assertTrue((output_dir / "background_audit_report.md").exists())

            with (output_dir / "background_matched_audit_summary.csv").open(newline="") as f:
                summary = list(csv.DictReader(f))
            cipro = next(r for r in summary if r["site"] == "A" and r["drug"] == "Cipro")
            self.assertEqual(cipro["adequacy_label"], "interpretable")
            self.assertIn("stratum_centered_auc", cipro)


if __name__ == "__main__":
    unittest.main()
