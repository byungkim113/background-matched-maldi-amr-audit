# Background-Matched MALDI-AMR Audit

**Can MALDI-TOF mass spectrometry genuinely predict antibiotic resistance — or does it just recognise resistant bacteria as a population?**

This repository provides a model-agnostic audit framework that separates *focal-drug resistance signal* from *resistant-population background* in MALDI-TOF antimicrobial resistance (AMR) prediction models. We apply it to a deep learning CNN and a multi-task LightGBM across four hospital sites, and cross-reference findings against publicly available whole-genome sequencing (WGS) and published ST131 proteomic biomarkers.

> **Key finding.** Ciprofloxacin resistance prediction in *E. coli* survives rigorous background matching and is grounded in detectable ST131 lineage signal (AUC = 0.91). Apparent performance for most beta-lactams collapses after matching, indicating background confounding rather than drug-specific biology.

---

## The Problem

Standard susceptibility testing (AST) takes 2–4 days. MALDI-TOF already identifies bacterial species in minutes from the same sample. A growing body of work claims MALDI models can *also* predict resistance, potentially collapsing the wait to under an hour.

However, resistant bacteria are not a random sample. Bacteria resistant to ciprofloxacin are often also resistant to ceftriaxone, ceftazidime, and cefepime — they cluster together in co-resistance blocks. A model that learns to recognise this resistant-population background will appear to predict each individual drug well, even if it has learned nothing drug-specific. Raw AUC alone cannot distinguish these two cases.

---

## The Framework

![Framework overview](outputs/final_framework_outputs/figure_5_framework_flow.png)

The audit proceeds in four steps:

1. **External transfer audit** — train at the source site, evaluate at held-out temporal and external hospital sites to measure raw portability and quantify random-CV inflation.
2. **Background-matched contrastive audit** — group isolates by their co-resistance profile (excluding the focal drug), then measure prediction AUC *within* each background stratum. This stratum-centered AUC is the key metric.
3. **Cross-resistance ecology network** — characterise the drug-label co-resistance structure that models can exploit, and link each drug-site result to its ecological block.
4. **Public WGS / proteomic cross-reference** — use publicly available WGS-linked MALDI data and published ST131 biomarkers to test whether the surviving residual signal has a plausible biological mechanism.

---

## Co-Resistance Landscape

The phi correlation heatmap below shows why background matching is necessary. Ciprofloxacin and norfloxacin are nearly always co-resistant (φ = 0.98). Ceftriaxone, ceftazidime, and cefepime form a tight ESBL/AmpC block (φ = 0.80–0.88). Any model trained on these labels can exploit this structure.

![Phi correlation heatmap](outputs/final_framework_outputs/figure_3_cross_resistance_phi_heatmap.png)

The 3D phi landscape makes the block structure explicit — two peaks (fluoroquinolone block and cephalosporin/ESBL block) separated by a mid-phi valley.

![Phi surface mesh and contour map](outputs/final_framework_outputs/figure_phi_surface_mesh.png)

---

## Primary Results

### Raw vs Background-Centered AUC

Each row is one site/drug combination. The circle is the raw external AUC; the square is the stratum-centered AUC after background matching. A large gap between the two indicates the raw performance was inflated by co-resistance background. Results are shown for both the CNN and LightGBM to confirm model-family independence.

![Raw vs background-centered AUC](outputs/final_framework_outputs/figure_1_raw_to_background_centered_auc.png)

**Reading the figure:**
- Blue (fluoroquinolone): ciprofloxacin retains substantial signal after matching at all well-powered sites.
- Orange (mixed beta-lactam): amoxicillin-clavulanic acid collapses to near-chance after matching across both model families.
- Red (cephalosporin/ESBL): ceftriaxone, ceftazidime, and cefepime are mostly low-powered or background-driven after matching.

### Signal Drop vs Matched Retention

The scatter below shows how much AUC was lost after matching (y-axis) against how many isolates survived matching (x-axis). Points in the upper-left corner had high raw AUC driven largely by background — few matched pairs survive and most apparent performance disappears.

![Signal drop vs matched retention](outputs/final_framework_outputs/figure_2_drop_vs_matched_retention.png)

### Summary Table

| Drug | Raw AUC (range) | Centered AUC (range) | Verdict |
|---|---|---|---|
| Ciprofloxacin | 0.67 – 0.82 | 0.60 – 0.70 | Survives background matching |
| Norfloxacin | 0.70 – 0.78 | 0.34 – 0.86 | Site-dependent; low matched support |
| Ceftriaxone | 0.62 – 0.87 | 0.59 – 0.73 | Partially survives; strongly attenuated |
| Ceftazidime | 0.61 – 0.84 | 0.37 – 0.72 | Mostly background-driven or underpowered |
| Cefepime | 0.65 – 0.87 | 0.39 – 0.66 | Background-driven or underpowered |
| Amox-Clavulanic acid | 0.54 – 0.69 | 0.49 – 0.58 | Collapses to chance after matching |

---

## Biological Mechanism: ST131 Lineage

To explain *why* ciprofloxacin prediction survives matching, we linked publicly available Basel UPEC MALDI spectra to whole-genome sequencing metadata. The ST131 lineage — a globally dominant high-risk *E. coli* clone — is detectable from MALDI peak features with AUC = 0.91. ST131 carries chromosomal fluoroquinolone resistance mutations at 20× higher odds than non-ST131 isolates.

Discriminative MALDI peak bins for each resistance target are enriched 2–3× for published ST131 proteomic biomarkers, well beyond a mass-matched permutation null.

![Public WGS and proteomic support](outputs/final_framework_outputs/figure_4_public_wgs_proteomic_support.png)

| Target | MALDI AUC | Biomarker enrichment | Empirical p |
|---|---|---|---|
| ST131 lineage | 0.906 | 3.11× | < 0.0001 |
| Ciprofloxacin resistance | 0.739 | 2.24× | 0.0068 |
| Ceftriaxone resistance | 0.650 | 2.64× | 0.0006 |
| All targets combined | — | 2.66× | < 0.0001 |

The interpretation: MALDI detects ST131 through conserved surface protein markers. Because ST131 almost always carries chromosomal fluoroquinolone resistance (QRDR mutations), lineage detection indirectly predicts ciprofloxacin resistance. For cephalosporins, resistance is mediated by mobile ESBL/AmpC plasmids that vary across hospitals — hence the site-dependent, often background-driven signal.

---

## Data Stack

```
DRIAMS (Basel, Switzerland)
  Four hospital sites: A-2018, B, C, D
  ~4,500 E. coli isolates with paired MALDI-TOF spectra and AST labels
  Six drugs: Ciprofloxacin, Norfloxacin, Amox-Clavulanic acid,
             Ceftriaxone, Ceftazidime, Cefepime
  ↓
  Raw spectra: not distributed (controlled access via DRIAMS portal)
  Isolate-level prediction CSVs: generated by export scripts below

Public Basel UPEC dataset (Cuenod et al.)
  407 urinary E. coli isolates
  Paired Bruker MALDI median-peak features + WGS metadata
  ST131 lineage calls, ciprofloxacin/ceftriaxone resistance labels
  ↓
  Bruker median-peak features: data_manifests/Bruker_csv_medianpeaks_df.csv
  WGS metadata bridge: data_manifests/upec_bruker_wgs_bridge.tsv

MARISMa (external stress-test)
  Independent Bruker MALDI E. coli snapshot
  Used to validate that DRIAMS-trained models do not blindly transfer
  ↓
  Snapshot summary: outputs/analysis_outputs/marisma_external_validation/
```

---

## Repository Layout

```
background-matched-maldi-amr-audit/
│
├── Mega_Model.py                        # CNN model engine (training, eval, export)
├── run_background_audit_framework.py    # Model-agnostic audit — the core method
│
├── scripts/
│   ├── run_training_ecoli6.py           # Train CNN on E. coli 6-drug panel
│   ├── run_training_clinical4.py        # Train CNN on clinical 4-pair profile
│   ├── run_lgbm_baselines.py            # Train multi-task LightGBM baselines
│   ├── export_mega_predictions_for_audit.py  # Export isolate-level CNN predictions
│   ├── run_background_audit.py          # Thin wrapper: run audit on default CSV
│   ├── run_public_upec_analysis.py      # WGS-linked MALDI + biomarker enrichment
│   ├── build_cross_resistance_network.py     # Phi/lift co-resistance network
│   ├── make_paper_figures.py            # Reproduce all manuscript figures
│   ├── make_final_framework_tables_figures.py
│   ├── marisma_end_to_end_kaggle.py     # MARISMa external stress-test
│   └── background_matched_contrastive_kaggle.py
│
├── data_manifests/                      # Public UPEC manifest and bridge tables
├── model_checkpoints/                   # Locked CNN checkpoint archive (5 seeds)
├── outputs/
│   ├── final_framework_outputs/         # Paper-ready tables and figures
│   └── analysis_outputs/               # Intermediate results and sub-analyses
├── manuscript/                          # Nature Communications LaTeX draft + figures
├── docs/                                # Reproduction guide, data availability notes
└── tests/                               # Regression tests for audit and model helpers
```

---

## Reproducing the Analysis

Full step-by-step instructions are in [`docs/reproduce.md`](docs/reproduce.md). Data availability and redistribution notes are in [`docs/data_availability.md`](docs/data_availability.md).

**1. Train the CNN model**

```bash
python scripts/run_training_ecoli6.py --data-root /path/to/driams
```

**2. Export isolate-level predictions**

```bash
python scripts/export_mega_predictions_for_audit.py \
  --run-dir runs/exp_ecoli_mechanism6_drugid_mae30
```

**3. Run the background-matched audit**

```bash
python scripts/run_background_audit.py \
  --predictions-csv runs/exp_ecoli_mechanism6_drugid_mae30/metrics/mega_predictions_long.csv \
  --output-dir outputs/background_audit
```

**4. Run the public UPEC WGS / proteomic support analysis**

```bash
python scripts/run_public_upec_analysis.py \
  --median-peaks data_manifests/Bruker_csv_medianpeaks_df.csv
```

**5. Reproduce all paper figures and tables**

```bash
python scripts/make_paper_figures.py
python scripts/make_final_framework_tables_figures.py
```

### Prediction CSV Format

The audit engine accepts any model's output as a long CSV with one row per isolate/drug:

```
isolate_id, site, year, organism, drug, label, prob
```

The framework builds co-resistance background signatures internally from the other drug labels available for the same isolate. An optional `model_name` column enables multi-model comparisons in a single run.

---

## Key Claims and Their Status

| Claim | Status | Supporting output |
|---|---|---|
| Ciprofloxacin prediction survives background matching | Confirmed (3/4 sites, both model families) | Table 1, Figure 1 |
| Amox-clavulanic acid is background-driven | Confirmed (3/4 sites, both model families) | Table 1, Figure 1 |
| Cephalosporin predictions are site-dependent | Partially confirmed | Table 2, Table 3 |
| Fluoroquinolone signal transfers across hospitals | Partially confirmed | Table 4 |
| MALDI encodes ST131 lineage | Confirmed (AUC = 0.91) | Table 7, Figure 4 |
| Resistance-associated peaks enriched for ST131 biomarkers | Confirmed (2–3×, p < 0.01) | Table 9, Figure 4 |
| Source thresholds do not transport to target hospitals | Confirmed | Table 6 |

---

## Included Outputs

Pre-computed outputs are provided so the audit results can be inspected without re-running training:

- `outputs/final_framework_outputs/` — all numbered paper tables (Table 1–10) and figures
- `outputs/analysis_outputs/cross_resistance_network/` — per-site phi heatmaps and network SVGs
- `outputs/analysis_outputs/upec_wgs_validation_outputs/` — ST131 / ciprofloxacin / ceftriaxone classification results
- `outputs/analysis_outputs/updated_proteomic_overlap_outputs/` — biomarker enrichment permutation results
- `model_checkpoints/mega_cnn_archive_2026-04-22/` — five-seed locked CNN checkpoint archive

These are analysis artifacts derived from the DRIAMS and public UPEC datasets. Raw spectra and AST tables are not included.

---

## Citation

If you use this framework or any of its outputs, please cite:

```bibtex
@software{background_matched_maldi_amr_audit,
  author    = {Richard-Y-W},
  title     = {Background-Matched MALDI-AMR Audit},
  year      = {2026},
  url       = {https://github.com/Richard-Y-W/background-matched-maldi-amr-audit}
}
```

See also [`CITATION.cff`](CITATION.cff) for the formal citation metadata.

---

## License

See [`LICENSE`](LICENSE) for terms. Raw DRIAMS data is subject to separate access conditions via the DRIAMS data portal. The public UPEC/Cuenod dataset is distributed under its original terms.
