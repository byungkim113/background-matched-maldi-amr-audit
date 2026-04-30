# Nature Communications Manuscript Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an Overleaf-ready Nature Communications-style manuscript package for the Background-Matched MALDI-AMR Audit paper.

**Architecture:** The manuscript is assembled from verified CSV outputs already committed in `outputs/`, while new figure PDFs are generated reproducibly by a single script. The LaTeX source keeps the main claims guarded: background-matched evaluation is the core method, cross-resistance structure and public WGS/proteomic analyses are supporting biological evidence.

**Tech Stack:** LaTeX, BibTeX, Python 3, pandas, reportlab.

---

### Task 1: Source and Journal Framing

**Files:**
- Create: `manuscript/references.bib`
- Modify: `manuscript/main.tex`

- [ ] Verify Nature Communications manuscript structure against official author guidance.
- [ ] Verify key MALDI-AMR and population-structure references using primary sources.
- [ ] Add citations for DRIAMS, Weis/Borgwardt MALDI-AMR, public UPEC WGS-linked data, and ST131 MALDI biomarkers.

### Task 2: Publication-Style Figure Builder

**Files:**
- Create: `scripts/make_ncomms_figures.py`
- Create outputs: `manuscript/figures/*.pdf`
- Create outputs: `manuscript/tables/*.tex`

- [ ] Read final audit tables from `outputs/final_framework_outputs/`.
- [ ] Create clean vector figures with consistent sizing, typography, line weights, and colorblind-safe colors.
- [ ] Export LaTeX-ready tables using `booktabs`.

### Task 3: Overleaf Manuscript

**Files:**
- Create: `manuscript/main.tex`

- [ ] Write a concise title, abstract, introduction, results, discussion, methods, data availability, code availability, and limitations.
- [ ] Include figures and tables in a Nature Communications-compatible order.
- [ ] State the claim as an evaluation framework, not definitive proof of DRIAMS clone identity or protein identity.

### Task 4: Verification

**Files:**
- Verify: `scripts/make_ncomms_figures.py`
- Verify: `manuscript/main.tex`

- [ ] Run the figure/table builder.
- [ ] Syntax-check the Python script.
- [ ] Compile LaTeX if `pdflatex` is available; otherwise record that Overleaf compilation is expected.

