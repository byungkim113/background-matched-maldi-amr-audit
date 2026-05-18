import unittest

from scripts import export_lgbm_predictions_for_audit as exporter


class ExportLgbmPredictionsForAuditTests(unittest.TestCase):
    def test_prediction_rows_include_model_variant_and_background_signature(self):
        active_lookup = {
            0: ("Escherichia coli", "Ciprofloxacin"),
            1: ("Escherichia coli", "Ceftriaxone"),
        }
        rows = [
            exporter.sample_to_prediction_row(
                ("/data/DRIAMS-A/binned_6000/2018/iso1_MALDI1.txt", 1, 0),
                prob=0.9,
                active_lookup=active_lookup,
                site="A-2018",
                model_name="LGBM-single",
            ),
            exporter.sample_to_prediction_row(
                ("/data/DRIAMS-A/binned_6000/2018/iso1_MALDI1.txt", 0, 1),
                prob=0.2,
                active_lookup=active_lookup,
                site="A-2018",
                model_name="LGBM-single",
            ),
        ]

        enriched = exporter.add_background_signatures(rows, active_lookup)

        self.assertEqual(enriched[0]["model_name"], "LGBM-single")
        self.assertEqual(enriched[0]["isolate_id"], "iso1_MALDI1")
        self.assertEqual(enriched[0]["year"], "2018")
        self.assertEqual(enriched[0]["background_signature"], "Ceftriaxone=S")
        self.assertEqual(enriched[1]["background_signature"], "Ciprofloxacin=R")
        self.assertEqual(enriched[1]["background_resistant_count"], 1)


if __name__ == "__main__":
    unittest.main()
