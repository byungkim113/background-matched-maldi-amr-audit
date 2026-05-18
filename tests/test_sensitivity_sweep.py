"""Tests for scripts/sensitivity_sweep.py."""
from __future__ import annotations

import csv
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import run_background_audit_framework as fw
from sensitivity_sweep import (
    _stratum_centered_auc,
    _summarize,
    run_sweep,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_rows(pairs: list[dict]) -> list[dict]:
    """Build minimal prediction rows with background signatures already set."""
    rows = []
    for p in pairs:
        rows.append(dict(
            uid=p["uid"], site=p.get("site", "Site-A"),
            year=p.get("year", "2021"), organism=p.get("organism", "Escherichia coli"),
            drug=p["drug"], label=p["label"], prob=p["prob"],
            model_name="test",
        ))
    return rows


def _with_sigs(rows: list[dict]) -> list[dict]:
    return fw.add_background_signatures(rows)


# ---------------------------------------------------------------------------
# Unit tests for _stratum_centered_auc
# ---------------------------------------------------------------------------

class TestStratumCenteredAuc:

    def _perfect_focal_rows(self) -> list[dict]:
        """
        Two strata with identical background labels.
        Within each stratum, model perfectly discriminates by focal label.
        Centered AUC should stay high.
        """
        rows = []
        for strat in ["S", "R"]:            # background stratum value (Norfloxacin status)
            for focal_r in [0, 1, 0, 1]:   # 2R + 2S per stratum
                uid = f"{strat}_{focal_r}_{len(rows)}"
                rows.append(dict(
                    uid=uid, site="S", year="2021",
                    organism="Escherichia coli", drug="Ciprofloxacin",
                    label=focal_r,
                    prob=0.85 if focal_r else 0.20,
                    model_name="t",
                    background_signature=f"Norfloxacin={strat}",
                ))
        return rows

    def _background_driven_rows(self) -> list[dict]:
        """
        Model scores correlate with background (Norfloxacin=R → high prob),
        not with focal label.  Centered AUC should be near 0.50.
        """
        rows = []
        for strat, prob_high in [("R", 0.80), ("S", 0.25)]:
            for focal_r in [0, 1, 0, 1]:
                uid = f"{strat}_{focal_r}_{len(rows)}"
                rows.append(dict(
                    uid=uid, site="S", year="2021",
                    organism="Escherichia coli", drug="Ceftriaxone",
                    label=focal_r,
                    prob=prob_high,         # same prob regardless of focal label
                    model_name="t",
                    background_signature=f"Norfloxacin={strat}",
                ))
        return rows

    def test_focal_signal_survives(self):
        rows = self._perfect_focal_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        assert result["centered_auc"] > 0.90, result["centered_auc"]
        assert result["n_valid_strata"] == 2
        assert result["matched_retention"] == 1.0

    def test_background_driven_collapses(self):
        rows = self._background_driven_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        # After centering, the background-driven signal should be near chance
        assert result["centered_auc"] < 0.60, result["centered_auc"]
        # raw AUC can still be > 0.5 because label is confounded with background
        assert not math.isnan(result["raw_auc"])

    def test_threshold_too_strict_returns_nan(self):
        rows = self._perfect_focal_rows()   # 2 pos + 2 neg per stratum
        result = _stratum_centered_auc(rows, min_t=3)
        assert math.isnan(result["centered_auc"])
        assert result["n_valid_strata"] == 0
        assert result["n_matched"] == 0

    def test_pairwise_accuracy_perfect(self):
        rows = self._perfect_focal_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        assert result["pairwise_accuracy"] == 1.0

    def test_pairwise_accuracy_background_driven(self):
        rows = self._background_driven_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        # With same prob within each stratum, 50% pairwise accuracy
        assert result["pairwise_accuracy"] == pytest.approx(0.5, abs=0.01)

    def test_matched_retention_fraction(self):
        rows = self._perfect_focal_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        assert result["matched_retention"] == 1.0

    def test_raw_auc_computed_for_all_rows(self):
        rows = self._perfect_focal_rows()
        result = _stratum_centered_auc(rows, min_t=2)
        assert math.isfinite(result["raw_auc"])


class TestRunSweep:

    def _make_sweep_rows(self) -> list[dict]:
        """Two isolates per (site, drug) combination with known background signatures."""
        base = []
        for site in ["Site-A", "Site-B"]:
            for strat in ["X=R|Y=S", "X=S|Y=S"]:
                for focal_r in [0, 1, 0, 1]:
                    base.append(dict(
                        uid=f"{site}_{strat}_{focal_r}_{len(base)}",
                        site=site, year="2021",
                        organism="Escherichia coli", drug="Ciprofloxacin",
                        label=focal_r,
                        prob=0.80 if focal_r else 0.25,
                        model_name="t",
                        background_signature=strat,
                    ))
        return base

    def test_sweep_produces_correct_threshold_count(self):
        rows = self._make_sweep_rows()
        thresholds = [2, 3]
        detail, summary = run_sweep(rows, thresholds)
        assert len(summary) == 2
        assert {r["min_stratum"] for r in summary} == {2, 3}

    def test_detail_covers_all_pairs(self):
        rows = self._make_sweep_rows()
        detail, _ = run_sweep(rows, [2])
        pair_keys = {(r["site"], r["organism"], r["drug"]) for r in detail}
        assert ("Site-A", "Escherichia coli", "Ciprofloxacin") in pair_keys
        assert ("Site-B", "Escherichia coli", "Ciprofloxacin") in pair_keys

    def test_summary_adequate_count_decreases_with_threshold(self):
        rows = self._make_sweep_rows()    # 2 pos + 2 neg per stratum
        _, summary = run_sweep(rows, [2, 3])
        s2 = next(r for r in summary if r["min_stratum"] == 2)
        s3 = next(r for r in summary if r["min_stratum"] == 3)
        # At threshold=3 no stratum qualifies (only 2 of each)
        assert int(s2["n_adequate"]) >= int(s3["n_adequate"])

    def test_summary_fields_present(self):
        rows = self._make_sweep_rows()
        _, summary = run_sweep(rows, [2])
        row = summary[0]
        for field in ["min_stratum", "n_drug_site_pairs", "n_adequate",
                      "mean_matched_retention", "macro_raw_auc",
                      "macro_centered_auc", "mean_delta"]:
            assert field in row, f"missing field: {field}"


class TestSummarize:

    def _detail(self, raw, cen, ret):
        return dict(min_stratum=3, site="S", organism="O", drug="D",
                    raw_auc=f"{raw:.4f}", matched_auc="", centered_auc=f"{cen:.4f}",
                    pairwise_accuracy="", pairwise_comparisons=10,
                    matched_retention=f"{ret:.4f}",
                    n_total=100, n_r=20, n_matched=50, n_matched_r=10, n_valid_strata=3)

    def test_mean_delta_sign(self):
        rows = [self._detail(0.75, 0.62, 0.8), self._detail(0.70, 0.65, 0.9)]
        s = _summarize(3, rows)
        expected_delta = (0.75 - 0.62 + 0.70 - 0.65) / 2
        assert float(s["mean_delta"]) == pytest.approx(expected_delta, abs=1e-4)

    def test_empty_centered_counted_as_inadequate(self):
        rows = [
            dict(min_stratum=3, site="S", organism="O", drug="D",
                 raw_auc="0.70", matched_auc="", centered_auc="",
                 pairwise_accuracy="", pairwise_comparisons=0,
                 matched_retention="0.00",
                 n_total=10, n_r=5, n_matched=0, n_matched_r=0, n_valid_strata=0),
        ]
        s = _summarize(3, rows)
        assert s["n_adequate"] == 0
        assert s["macro_centered_auc"] == ""


class TestCLI:
    """Integration test: generate example predictions, run sweep via main()."""

    def test_cli_produces_output_files(self, tmp_path):
        import subprocess
        example_csv = ROOT / "example_predictions.csv"
        if not example_csv.exists():
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "generate_example_data.py")],
                check=True, capture_output=True,
            )
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "sensitivity_sweep.py"),
             "--predictions-csv", str(example_csv),
             "--output-dir", str(tmp_path),
             "--thresholds", "2,3"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "sensitivity_detail.csv").exists()
        assert (tmp_path / "sensitivity_summary.csv").exists()
        assert (tmp_path / "sensitivity_summary.md").exists()

    def test_summary_has_two_threshold_rows(self, tmp_path):
        import subprocess
        example_csv = ROOT / "example_predictions.csv"
        if not example_csv.exists():
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "generate_example_data.py")],
                check=True, capture_output=True,
            )
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "sensitivity_sweep.py"),
             "--predictions-csv", str(example_csv),
             "--output-dir", str(tmp_path),
             "--thresholds", "2,5"],
            check=True, capture_output=True,
        )
        with (tmp_path / "sensitivity_summary.csv").open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 2
        thresholds = {int(r["min_stratum"]) for r in rows}
        assert thresholds == {2, 5}

    def test_site_robustness_flag(self, tmp_path):
        import subprocess
        example_csv = ROOT / "example_predictions.csv"
        if not example_csv.exists():
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "generate_example_data.py")],
                check=True, capture_output=True,
            )
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "sensitivity_sweep.py"),
             "--predictions-csv", str(example_csv),
             "--output-dir", str(tmp_path),
             "--thresholds", "3",
             "--site-robustness"],
            check=True, capture_output=True,
        )
        assert (tmp_path / "sensitivity_site_robustness.csv").exists()
        with (tmp_path / "sensitivity_site_robustness.csv").open() as f:
            rows = list(csv.DictReader(f))
        sites = {r["site"] for r in rows}
        assert len(sites) >= 2, "Expected at least 2 sites in robustness table"


import pytest
