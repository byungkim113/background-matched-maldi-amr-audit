# Deployment Reliability Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reproducible supplementary tools for deployment decisions, calibration, temporal reliability monitoring, and falsification controls for the MALDI-AMR background-matched audit.

**Architecture:** Use standalone standard-library analysis scripts that consume the existing long prediction CSV and background audit tables. The final framework builder imports the generated CSV outputs and exports them as numbered final tables with short Markdown summaries.

**Tech Stack:** Python standard library for new analysis scripts and unit tests; existing pandas/Pillow final-framework generator for final table export.

---

### Task 1: Calibration Metrics

**Files:**
- Create: `scripts/calibration_analysis.py`
- Create: `tests/test_calibration_analysis.py`
- Output: `outputs/analysis_outputs/calibration_analysis/calibration_summary.csv`

- [ ] Write tests for Brier score, expected calibration error, threshold metrics, and CLI output.
- [ ] Verify tests fail before the script exists.
- [ ] Implement grouped calibration analysis over `site + organism + drug`.
- [ ] Verify tests pass.

### Task 2: Temporal Reliability Monitor

**Files:**
- Create: `scripts/temporal_reliability_audit.py`
- Create: `tests/test_temporal_reliability_audit.py`
- Output: `outputs/analysis_outputs/temporal_reliability_audit/temporal_reliability.csv`

- [ ] Write tests for multi-period drift flags and single-period insufficiency flags.
- [ ] Verify tests fail before the script exists.
- [ ] Implement grouped period metrics and baseline-to-period drift calculations.
- [ ] Verify tests pass.

### Task 3: Falsification Controls

**Files:**
- Create: `scripts/falsification_controls.py`
- Create: `tests/test_falsification_controls.py`
- Output: `outputs/analysis_outputs/falsification_controls/falsification_controls.csv`

- [ ] Write tests for within-background label shuffle, burden-only AUC, and score-vs-burden association controls.
- [ ] Verify tests fail before the script exists.
- [ ] Implement deterministic permutation controls and background-only comparisons.
- [ ] Verify tests pass.

### Task 4: Deployment Decision Framework

**Files:**
- Create: `scripts/deployment_decision_framework.py`
- Create: `tests/test_deployment_decision_framework.py`
- Output: `outputs/analysis_outputs/deployment_decision_framework/deployment_decision_rules.csv`
- Output: `outputs/analysis_outputs/deployment_decision_framework/deployment_readiness_by_pair.csv`

- [ ] Write tests for raw/matched/calibration decision categories.
- [ ] Verify tests fail before the script exists.
- [ ] Implement static decision rules and pair-level deployment actions.
- [ ] Verify tests pass.

### Task 5: Final Framework Integration

**Files:**
- Modify: `scripts/make_final_framework_tables_figures.py`
- Output: `outputs/final_framework_outputs/table_12_deployment_decision_rules.csv`
- Output: `outputs/final_framework_outputs/table_13_calibration_summary.csv`
- Output: `outputs/final_framework_outputs/table_14_temporal_reliability_audit.csv`
- Output: `outputs/final_framework_outputs/table_15_falsification_controls.csv`
- Output: `outputs/final_framework_outputs/table_16_deployment_readiness_by_pair.csv`

- [ ] Add generated analysis CSVs to `INPUTS`.
- [ ] Export final numbered tables and update `final_framework_summary.md`.
- [ ] Run the four new scripts and regenerate final framework outputs.
- [ ] Run unit tests and syntax checks before committing.
