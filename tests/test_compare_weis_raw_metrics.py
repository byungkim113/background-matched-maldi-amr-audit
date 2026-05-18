import importlib.util
import pathlib
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compare_weis_raw_metrics.py"


def load_module():
    spec = importlib.util.spec_from_file_location("compare_weis_raw_metrics", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CompareWeisRawMetricsTests(unittest.TestCase):
    def test_compare_matches_reference_by_model_site_pair_and_seed(self):
        module = load_module()
        raw = [
            {
                "model": "lr",
                "train_site": "DRIAMS-A",
                "test_site": "DRIAMS-C",
                "species": "Escherichia coli",
                "drug": "Ciprofloxacin",
                "seed": 35,
                "auroc": 0.75,
                "auprc": 0.6,
                "accuracy": 0.8,
            }
        ]
        reference = {
            ("lr", "DRIAMS-A", "DRIAMS-C", "Escherichia coli", "Ciprofloxacin", "35"): {
                "auroc": 0.7500001,
                "auprc": 0.6,
                "accuracy": 0.8,
                "_reference_path": "ref.json",
            }
        }

        rows = module.compare(raw, reference, tolerance=1e-5)

        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0]["reference_found"])
        self.assertTrue(rows[0]["within_tolerance"])

    def test_load_reference_results_skips_unrelated_json(self):
        module = load_module()
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "bad.json").write_text('{"not": "a result"}')
            (root / "good.json").write_text(
                '{"model":"lr","train_site":"DRIAMS-A","test_site":"DRIAMS-B",'
                '"species":"Escherichia coli","antibiotic":"Ceftriaxone","seed":35,'
                '"auroc":0.7,"auprc":0.5,"accuracy":0.8}'
            )

            index = module.load_reference_results(root)

        self.assertEqual(list(index), [("lr", "DRIAMS-A", "DRIAMS-B", "Escherichia coli", "Ceftriaxone", "35")])


if __name__ == "__main__":
    unittest.main()
