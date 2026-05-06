# Final Background-Matched MALDI-AMR Framework Outputs

## Core Claim

MALDI-TOF AMR models should be interpreted as background-sensitive predictors: raw external AUC can mix focal-drug signal with resistant-population, lineage, and co-resistance background. Background-matched auditing tests how much apparent resistance prediction survives after that background is controlled.

## Key Results From The Current Artifact Set

- In interpretable Amox-Clav rows, CNN mean raw-minus-centered drop is 0.084; LGBM multi mean drop is 0.116.
- In interpretable Cipro rows, CNN mean background-centered AUC is 0.648; LGBM multi mean background-centered AUC is 0.618.
- Low-retention cephalosporin rows are explicitly flagged rather than overclaimed.
- Co-resistance stratification reports 907 strata, including 166 interpretable strata across burden bins and exact background signatures.
- Deployment readiness assigns 5 action categories across 23 audited pair/site rows.
- Temporal reliability monitoring reports 23/23 rows as single-period only; current DRIAMS export cannot estimate reliability duration without future periods.
- Falsification controls report 0/23 pair/site rows exceeding both burden-only and within-background shuffle controls; 3/23 exceed the shuffle null, while burden-only controls remain competitive in 23/23 rows.
- Public WGS-linked Bruker MALDI data show ST131 AUC=0.906, higher than Cipro-R and Ceftriaxone-R peak-only AUCs.
- Published ST131 biomarker enrichment is strongest for ST131 itself (3.11x) and remains significant for Cipro-R and Ceftriaxone-R discriminative peaks.

## What These Tables/Figures Are For

- Table 1: paper-ready primary background-matched audit with bootstrap CIs and adequacy labels.
- Table 2: CNN vs LGBM multi comparison, showing the effect is not just a neural-network artifact.
- Table 3: ecology-aware interpretation, linking background sensitivity to co-resistance blocks.
- Table 4: locked ecoli_mechanism6 transfer prediction assessment.
- Table 5: top cross-resistance network edges.
- Table 6: block-level source and external transfer behavior.
- Table 7-10 files: public WGS-linked lineage/resistance support and proteomic biomarker enrichment.
- Table 11: co-resistance stratification showing where focal-drug signal lives inside burden/signature strata.
- Table 12: deployment decision rules for pass/fail/underpowered/calibration cases.
- Table 13: calibration, Brier score, ECE, and threshold metrics.
- Table 14: temporal reliability and recalibration monitor.
- Table 15: falsification controls against burden-only and within-background shuffled labels.
- Table 16: pair-level deployment readiness actions.

## Figure Captions

Figure 1. Raw external AUC versus background-centered AUC for CNN and LGBM multi. Lines show how much performance remains after matching/centering by co-resistance background.
Figure 2. Signal drop versus matched retention. This distinguishes interpretable collapses from sparse matched strata.
Figure 3. Cross-resistance phi heatmap. Strong drug-drug blocks show the label ecology that AMR models can exploit.
Figure 4. Public WGS-linked support. ST131 is strongly predictable from MALDI peaks, and resistance-associated peaks are enriched for published ST131 biomarkers.
Figure 5. Framework schematic.

## Cautious Claims

- We can claim background sensitivity and the need for background-matched evaluation.
- We can claim the current evidence supports lineage/co-resistance background as part of the MALDI-AMR signal.
- We should not claim direct protein identity for DRIAMS saliency peaks or prove ST131 detection inside DRIAMS without WGS labels.

## Output Files

- `table_10_top_published_st131_peak_overlaps.csv`
- `table_11_co_resistance_stratification.csv`
- `table_12_deployment_decision_rules.csv`
- `table_13_calibration_summary.csv`
- `table_14_temporal_reliability_audit.csv`
- `table_15_falsification_controls.csv`
- `table_16_deployment_readiness_by_pair.csv`
- `table_1_primary_background_matched_audit.csv`
- `table_2_cnn_vs_lgbm_multi_background_audit.csv`
- `table_3_ecology_interpretation.csv`
- `table_4_transfer_prediction_assessment.csv`
- `table_5_top_cross_resistance_edges.csv`
- `table_6_ecoli_block_site_summary.csv`
- `table_7_public_wgs_maldi_auc.csv`
- `table_8_public_wgs_st131_resistance_associations.csv`
- `table_9_published_st131_biomarker_enrichment.csv`
- `figure_1_raw_to_background_centered_auc.png`
- `figure_2_drop_vs_matched_retention.png`
- `figure_3_cross_resistance_phi_heatmap.png`
- `figure_4_public_wgs_proteomic_support.png`
- `figure_5_framework_flow.png`
