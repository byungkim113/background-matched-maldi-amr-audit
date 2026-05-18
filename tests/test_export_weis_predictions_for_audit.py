import importlib.util
import json
import pathlib
import sys
import tempfile
import types
import unittest

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[1]
EXPORTER_PATH = ROOT / "scripts" / "export_weis_predictions_for_audit.py"


def load_exporter():
    if "sklearn.metrics" not in sys.modules:
        sklearn = types.ModuleType("sklearn")
        metrics = types.ModuleType("sklearn.metrics")

        def _unused_metric(*args, **kwargs):
            raise AssertionError("metric functions are not needed for exporter unit tests")

        metrics.accuracy_score = _unused_metric
        metrics.average_precision_score = _unused_metric
        metrics.roc_auc_score = _unused_metric
        sklearn.metrics = metrics
        sys.modules.setdefault("sklearn", sklearn)
        sys.modules.setdefault("sklearn.metrics", metrics)

    spec = importlib.util.spec_from_file_location("export_weis_predictions_for_audit", EXPORTER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FakeIndex:
    def __init__(self, labels):
        self.labels = list(labels)

    def get_indexer(self, labels):
        lookup = {label: idx for idx, label in enumerate(self.labels)}
        return np.asarray([lookup.get(label, -1) for label in labels], dtype=int)

    def astype(self, dtype):
        return np.asarray([dtype(label) for label in self.labels], dtype=object)


class FakeMeta:
    def __init__(self, labels):
        self.index = FakeIndex(labels)
        self.labels = list(labels)

    def __len__(self):
        return len(self.labels)


class ExportWeisPredictionsTests(unittest.TestCase):
    def test_all_external_row_policy_keeps_every_eligible_external_row(self):
        exporter = load_exporter()
        meta = FakeMeta(["a", "b", "c", "d"])

        def split_fn(*args, **kwargs):
            raise AssertionError("stratification should not be called for all-row external scoring")

        selected = exporter.select_external_indices(
            meta,
            split_fn,
            drug="Ciprofloxacin",
            seed=35,
            policy="all",
            name="DRIAMS-B/Ciprofloxacin/test",
        )

        np.testing.assert_array_equal(selected, np.arange(4))

    def test_stratified_external_row_policy_preserves_legacy_split_behavior(self):
        exporter = load_exporter()
        meta = FakeMeta(["iso-a", "iso-b", "iso-c", "iso-d"])

        def split_fn(frame, antibiotic, random_state):
            self.assertEqual(antibiotic, "Ciprofloxacin")
            self.assertEqual(random_state, 35)
            return ["iso-a"], ["iso-b", "iso-d"]

        selected = exporter.select_external_indices(
            meta,
            split_fn,
            drug="Ciprofloxacin",
            seed=35,
            policy="stratified",
            name="DRIAMS-B/Ciprofloxacin/test",
        )

        np.testing.assert_array_equal(selected, np.asarray([1, 3]))

    def test_custom_panel_uses_requested_species_and_drugs(self):
        exporter = load_exporter()
        args = type(
            "Args",
            (),
            {
                "panel": "custom",
                "species": "Escherichia coli",
                "drugs": "Ciprofloxacin,Cefepime",
            },
        )()

        self.assertEqual(
            exporter.resolve_species_drug_panels(args),
            [("Escherichia coli", ["Ciprofloxacin", "Cefepime"])],
        )

    def test_weis_core_panel_matches_original_repository_pair_list(self):
        exporter = load_exporter()
        args = type(
            "Args",
            (),
            {
                "panel": "weis-core",
                "species": "ignored",
                "drugs": "ignored",
            },
        )()

        panels = dict(exporter.resolve_species_drug_panels(args))
        self.assertIn("Escherichia coli", panels)
        self.assertIn("Staphylococcus aureus", panels)
        self.assertIn("Klebsiella pneumoniae", panels)
        self.assertIn("Piperacillin-Tazobactam", panels["Escherichia coli"])
        self.assertIn("Amoxicillin-Clavulanic acid", panels["Klebsiella pneumoniae"])

    def test_progress_helpers_report_pair_fraction_and_status_files(self):
        exporter = load_exporter()
        species_panels = [
            ("Escherichia coli", ["Ciprofloxacin", "Cefepime"]),
            ("Staphylococcus aureus", ["Oxacillin"]),
        ]

        total = exporter.planned_pair_count(["DRIAMS-B", "DRIAMS-C"], species_panels)
        self.assertEqual(total, 6)
        self.assertEqual(
            exporter.format_pair_progress(3, total, "DRIAMS-C", "Escherichia coli", "Cefepime"),
            "[pair 3/6 50.0%] DRIAMS-C | Escherichia coli / Cefepime",
        )

        with tempfile.TemporaryDirectory() as tmp:
            out = pathlib.Path(tmp)
            exporter.write_status(
                out,
                stage="training",
                pair_index=3,
                total_pairs=6,
                test_site="DRIAMS-C",
                species="Escherichia coli",
                drug="Cefepime",
                rows_written=123,
            )

            status = json.loads((out / "STATUS.json").read_text())
            self.assertEqual(status["stage"], "training")
            self.assertEqual(status["pair_index"], 3)
            self.assertEqual(status["total_pairs"], 6)
            self.assertEqual(status["percent_pairs_complete"], 50.0)
            self.assertEqual(status["rows_written"], 123)
            self.assertIn("training", (out / "CURRENT_STAGE.txt").read_text())

    def test_git_metadata_handles_non_git_source_directory(self):
        exporter = load_exporter()
        with tempfile.TemporaryDirectory() as tmp:
            metadata = exporter.git_metadata(pathlib.Path(tmp))
        self.assertFalse(metadata["source_is_git_checkout"])
        self.assertIsNone(metadata["source_commit"])
        self.assertIsNone(metadata["source_branch"])


if __name__ == "__main__":
    unittest.main()
