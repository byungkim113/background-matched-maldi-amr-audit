import csv
import pathlib
import tempfile
import unittest

from scripts import build_model_class_matrix as matrix


class BuildModelClassMatrixTests(unittest.TestCase):
    def write_summary(self, path: pathlib.Path) -> None:
        fields = [
            "site",
            "organism",
            "drug",
            "raw_auc",
            "stratum_centered_auc",
            "pairwise_accuracy",
            "matched_retention",
            "adequacy_label",
            "interpretation_category",
            "n_total",
            "n_r",
            "n_matched",
            "n_matched_r",
            "n_valid_strata",
            "min_pos_per_stratum",
            "min_neg_per_stratum",
        ]
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerow(
                {
                    "site": "A-2018",
                    "organism": "Staphylococcus aureus",
                    "drug": "Oxacillin",
                    "raw_auc": "0.83",
                    "stratum_centered_auc": "0.76",
                    "pairwise_accuracy": "0.77",
                    "matched_retention": "0.60",
                    "adequacy_label": "interpretable",
                    "interpretation_category": "focal_signal_retained",
                    "n_total": "100",
                    "n_r": "20",
                    "n_matched": "60",
                    "n_matched_r": "15",
                    "n_valid_strata": "4",
                    "min_pos_per_stratum": "3",
                    "min_neg_per_stratum": "3",
                }
            )

    def test_builds_completed_and_missing_rows_from_model_specs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            summary_path = tmp / "summary.csv"
            self.write_summary(summary_path)

            specs = [
                matrix.ModelSpec(
                    model_class="CNN/Mega",
                    model_variant="multi-task CNN",
                    organism="Staphylococcus aureus",
                    drug="Oxacillin",
                    summary_path=summary_path,
                    scope="Sa/Oxa panel",
                ),
                matrix.ModelSpec(
                    model_class="LightGBM",
                    model_variant="single-task",
                    organism="Staphylococcus aureus",
                    drug="Oxacillin",
                    summary_path=tmp / "missing.csv",
                    scope="Sa/Oxa panel",
                ),
            ]

            rows = matrix.build_matrix_rows(specs)

        self.assertEqual(len(rows), 2)
        complete, missing = rows
        self.assertEqual(complete["status"], "complete")
        self.assertEqual(complete["model_class"], "CNN/Mega")
        self.assertEqual(complete["site"], "A-2018")
        self.assertEqual(complete["centered_auc"], "0.760000")
        self.assertEqual(complete["valid_strata"], "4")
        self.assertEqual(missing["status"], "missing")
        self.assertEqual(missing["model_variant"], "single-task")

    def test_completed_matrix_writes_completion_note_instead_of_missing_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "run_missing_lgbm_commands.md"
            matrix.write_missing_commands(path, [{"status": "complete"}])

            text = path.read_text()

        self.assertIn("LGBM Model-Class Cells Completed", text)
        self.assertNotIn("Missing LGBM Model-Class Cells", text)


if __name__ == "__main__":
    unittest.main()
