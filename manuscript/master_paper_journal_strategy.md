# Master Paper Journal And Display Strategy

This note records the editorial choices behind `master_research_paper.md`. It is not part of the manuscript.

## Target Interpretation

The user asked about Journal of Clinical Microbiology and "NPJ Biology." There is no single journal called NPJ Biology in the Nature Portfolio npj series. I interpreted the closest relevant target as `npj Systems Biology and Applications`, because the paper is a computational and systems-level analysis of clinical bacterial phenotypes, model behavior, co-resistance networks, and lineage-associated spectral structure.

## Journal Fit

### Journal of Clinical Microbiology

Best fit if the paper is framed as a clinical-laboratory validation guardrail for MALDI-TOF AMR prediction. The strongest JCM pitch is:

> Before MALDI-TOF AMR models are used for clinical decisions, external AUC should be decomposed by co-resistance background to determine whether prediction survives within clinically observable resistant-population contexts.

This route should keep the language practical: model output, AST panel, external site, audit decision, deployment risk. The paper should not over-emphasize method novelty at the expense of clinical microbiology value.

### npj Systems Biology and Applications

Best fit if the paper is framed as a systems-biology/methods contribution:

> MALDI-TOF spectra, AST labels, co-resistance ecology, and bacterial lineage form a coupled system. The background-matched audit quantifies which part of AMR prediction survives after controlling one observable layer of that system.

This route can foreground cross-resistance networks, public WGS-linked lineage evidence, and model-agnostic evaluation. It should still be concise and avoid becoming a software paper.

## Display Rule Applied

Both journal styles favor a clean main story and supplemental support. The main manuscript should stand on its own, while large tables, secondary model checks, stress tests, and implementation details belong in supplementary material.

## Revision Response To External Critique

The Overleaf draft now treats the manuscript as a background-sensitivity audit, not a proof of causal confounding. The title and major claims avoid saying that the model definitively learned the resistance mechanism or that background caused the raw AUC. The central claim is narrower: raw MALDI-TOF AMR performance can be sensitive to resistant-population background, and the audit quantifies how much focal-drug ranking survives matched co-resistance backgrounds.

The closest existing method is Weis et al.'s hierarchical stratification. The revised Introduction and Discussion now compare it directly with this audit: hierarchical stratification is spectrum-level and split-level, whereas this audit is prediction-table-level and post hoc. That makes it complementary and easier for an independent auditor to run when only isolate-level predictions and AST labels are available.

The paper now explicitly frames the audit as a domain-specific extension to TRIPOD+AI/PROBAST-style prediction-model reporting. It also foregrounds Yu et al. as the WGS precedent: they showed population-structure confounding for WGS AMR prediction; this manuscript tests and audits the analogous risk for MALDI-TOF AMR prediction.

The ciprofloxacin/norfloxacin circularity concern was checked empirically in `scripts/revision_audit_checks.py`. The critique's expected pattern was not observed: valid ciprofloxacin strata at A-2018, DRIAMS-C, and DRIAMS-D were not mostly norfloxacin-susceptible; norfloxacin was unknown in those valid strata. This weakens a fluoroquinolone-pair interpretation, but it strengthens the primary cipro result in another way: the interpretable cipro rows are effectively cross-class tests of ciprofloxacin discrimination within beta-lactam and cephalosporin co-resistance backgrounds. The manuscript now presents this as a positive reinterpretation in the Results and Discussion, while still noting that norfloxacin missingness prevents a fully controlled fluoroquinolone-pair analysis.

The n >= 3 threshold is now supported by descriptive sensitivity analyses at n >= 5 and n >= 10. The ciprofloxacin versus amoxicillin-clavulanic acid contrast persists at stricter thresholds, but the manuscript describes DRIAMS-D ciprofloxacin as weak retained signal rather than strong retention. DRIAMS-B ciprofloxacin is treated as statistically uninformative because it has one valid stratum and only 25 matched isolates.

The critique's recommendation to add `Staphylococcus aureus`/oxacillin is scientifically right, but the required Mega/CNN isolate-level export is not available in the current artifact set. The manuscript now names this as the most important next extension rather than pretending the result exists. The small Weis LR official-panel oxacillin outputs should stay supplementary until there is a full primary-model export.

## Recommended Main Display Items

1. **Figure 1: Audit framework and metric logic.**
   - Use `manuscript/figures/figure_1_framework.pdf`.
   - Purpose: define raw AUC, matched AUC, background-centered AUC, pairwise within-background accuracy, and matched retention.

2. **Figure 2: Primary DRIAMS background audit.**
   - Use `manuscript/figures/figure_2_primary_background_audit.pdf`.
   - Purpose: show ciprofloxacin retention and amoxicillin-clavulanic acid collapse without burying the reader in all drugs.

3. **Figure 3: Biological support from co-resistance and public WGS-linked MALDI.**
   - Prefer a combined figure using current `figure_4_cross_resistance_network.pdf` plus `figure_5_public_wgs_proteomic_support.pdf`, or use the current WGS/proteomic figure as the main display and move the full cross-resistance heatmap to supplement.
   - Purpose: demonstrate that exploitable co-resistance blocks exist and that MALDI spectra strongly encode ST131 lineage.

4. **Table 1: Primary interpretable audit rows.**
   - Use a shortened version of `outputs/final_framework_outputs/table_1_primary_background_matched_audit.csv`.
   - Include only interpretable/cautionary primary contrast rows plus the main cephalosporin support rows. Move the full 23-row table to supplement.

## Move To Supplement

- Model-family comparison: important, but secondary to the primary claim. Move current Figure 3 and Table 2 to Supplementary Figure/Table.
- Falsification controls: important guardrail, but too technical for main flow. Move current Figure 6 and Table 15 to supplement; mention the headline in Results and Discussion.
- Weis/Borgwardt official LR parity and six-drug stress test: supplement only. It strengthens reproducibility and audit portability, but the main paper should not look like it is mainly a Weis replication paper.
- MARISMa external stress test: supplement only. It is useful as a boundary condition and failed transfer case, not part of the main positive argument.
- Deployment decision flow: supplement or final discussion schematic only if the target journal wants translational guidance. Otherwise, state the decision rules in text.
- Calibration, temporal reliability, and readiness tables: supplement.

## Claim Guardrails

- Do claim that raw MALDI-TOF AMR performance can contain resistant-population background signal.
- Do claim that background-matched evaluation distinguishes retained within-background signal from background-sensitive transfer.
- Do claim that public WGS-linked MALDI data make lineage-associated spectral shortcuts biologically plausible.
- Do not claim direct ST131 detection inside DRIAMS without DRIAMS-linked WGS/MLST.
- Do not claim protein identity for mass-bin overlaps without MS/MS or LC-MS/MS.
- Do not claim exact Weis replication beyond the official 8-row LR parity panel.
- Do not present MARISMa as successful external validation.

## Online Sources Consulted

Primary journal and editorial sources checked on 2026-05-12:

- Nature Portfolio, `npj Systems Biology and Applications`, For authors: https://www.nature.com/npjsba/for-authors-and-referees/submisions
  - Initial submissions do not need special formatting if the study is suitable for editorial assessment and peer review.
  - Figures may be inserted near the relevant text for reviewer readability.
- Nature Portfolio, `npj Systems Biology and Applications`, Aims and Scope: https://www.nature.com/npjsba/aims
  - The journal considers computational and mathematical approaches to complex biological systems, including network biology and disease applications.
- Nature formatting guide: https://www.nature.com/nature/for-authors/formatting-guide
  - Nature-style Articles use a concise summary paragraph, short subheadings, and a limited main display set.
  - Typical biological/clinical Articles with 5-6 modest display items are around 4,300 words, so the master draft was kept near that scale.
- Nature Portfolio reporting standards for `npj Systems Biology and Applications`: https://www.nature.com/npjsba/editorial-policies/reporting-standards
  - Data, materials, code, and protocols must be made available or restrictions must be disclosed.
  - Data availability statements should identify primary and referenced datasets with access information and accession identifiers where relevant.
- ICMJE manuscript preparation recommendations: https://www.icmje.org/recommendations/browse/manuscript-preparation/preparing-for-submission.html
  - Original biomedical research usually follows IMRAD structure.
  - References should prioritize original sources, and published articles should cite persistent identifiers for datasets.
- ASM announcement that ASM journals accept initial submissions in any format: https://asm.org/press-releases/2018/asm-journals-now-accepting-submissions-in-any-form
  - This supports focusing the first submission on clarity and scientific structure rather than premature formatting.
- ASM scientific writing guidance: https://asm.org/articles/2019/november/get-your-scientific-paper-off-the-ground
  - The paper was rebuilt from results outline to figures, title/abstract, Results, Introduction, then Discussion.
- ASM figure guidance: https://asm.org/articles/2019/april/use-clear-figures-to-tell-a-story-with-your-data
  - The main display plan prioritizes a small number of figures that each carry one conclusion, with detailed audit tables in supplement.
- Journal of Clinical Microbiology archived instructions: https://pmc.ncbi.nlm.nih.gov/articles/PMC2620845/
  - JCM scope is clinical microbiology, diagnosis, and epidemiology of infection.
  - Supplemental material is appropriate for large or complex datasets, but main conclusions should be supported by the manuscript without requiring the supplement.
  - The current journals.asm.org pages were not fully accessible through the crawler, so current ASM-wide author resources were paired with this archived JCM-specific instruction page.
