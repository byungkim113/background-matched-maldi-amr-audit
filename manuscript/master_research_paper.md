# Background-Matched Evaluation Reveals Resistant-Population Confounding In MALDI-TOF AMR Prediction

Running title: Background-matched MALDI-AMR evaluation

Authors: Byung Kim, Yuchen Wang, Nicolas Samuel Khoury-Levy

Affiliations: To be completed before submission.

Correspondence: To be completed before submission.

## Abstract

Machine learning models trained on matrix-assisted laser desorption/ionization time-of-flight (MALDI-TOF) spectra can predict antimicrobial resistance (AMR), but raw external AUC does not reveal whether prediction reflects focal resistance biology or the spectral background of resistant bacterial populations. We developed a model-agnostic background-matched audit that evaluates isolate-level predictions within matched co-resistance strata. Applied to *Escherichia coli* DRIAMS models trained at DRIAMS-A and tested temporally and across external hospitals, ciprofloxacin predictions retained within-background discrimination after co-resistance adjustment, with background-centered AUCs of 0.703, 0.646, and 0.596 at A-2018, DRIAMS-C, and DRIAMS-D. In contrast, amoxicillin-clavulanic acid weakened or collapsed toward chance at the same sites, with centered AUCs of 0.541, 0.497, and 0.486. The same attenuation pattern appeared in convolutional neural network and LightGBM outputs, indicating that background sensitivity is not specific to one architecture. Cross-resistance networks showed strong exploitable label structure, and independent public WGS-linked MALDI data showed that ST131 lineage was more predictable from MALDI peak features than ciprofloxacin or ceftriaxone resistance. Background-matched evaluation provides a practical guardrail for MALDI-TOF AMR studies by distinguishing aggregate transfer from focal-drug signal that survives within resistant-population background.

## Importance

MALDI-TOF instruments are already present in clinical microbiology laboratories, making AMR prediction from routine spectra attractive. The central risk is that a model may appear to predict resistance while partly detecting the resistant clone, co-resistance block, or hospital ecology in which resistance is concentrated. That distinction matters for deployment: a shortcut that works in one hospital can fail when local resistant populations change. This study converts that concern into an operational audit requiring only ordinary model predictions and AST labels. The audit does not prove mechanism, but it shows whether a claimed drug prediction survives within matched co-resistance backgrounds. This makes MALDI-TOF AMR evaluation more interpretable, more reproducible, and safer to compare across sites.

## Introduction

Rapid antimicrobial susceptibility information is central to effective treatment and stewardship. MALDI-TOF mass spectrometry is already used at scale for bacterial identification, and several studies have shown that the same spectra can support machine-learning prediction of AMR. The DRIAMS resource was especially important because it paired clinical MALDI-TOF spectra with susceptibility labels across multiple hospitals and demonstrated that AMR prediction can transfer beyond a single laboratory under some organism-drug settings \cite{weis2022driams,driamsdryad}. That result created a realistic translational question: not whether MALDI-TOF spectra contain predictive information, but what biological and population signals that information represents.

A resistance label is not an isolated molecular state. In clinical bacterial populations, resistance is often concentrated in disseminated lineages, plasmid backgrounds, and co-resistance blocks. These structures vary by hospital and sampling period. A model trained to separate resistant from susceptible isolates can therefore succeed by learning any spectral correlate of resistance, including lineage-associated protein peaks that are not the focal resistance mechanism. This risk is not hypothetical. Population structure has been shown to bias machine-learning AMR prediction from whole-genome sequence data \cite{yu2025population}, and MALDI-TOF spectra can encode phylogenetic and strain background because routine spectra contain abundant conserved proteins whose masses vary across lineages.

For clinical translation, the usual question, "What is the external AUC?", is incomplete. A high external AUC can arise when the training and test sites share a resistant subpopulation, even if the model is using background structure rather than focal-drug biology. Conversely, an external AUC can decline when the resistant population changes, even if the model learned a real signal that is not represented in the new site. AUC alone does not distinguish these cases. Existing efforts, including hierarchical stratification by spectral similarity, have begun to address this issue \cite{weis2022stratification}, but a simple audit that can be applied to any already-generated prediction table has been missing.

We introduce a background-matched MALDI-AMR audit. For each focal drug, the method groups isolates by the AST labels of the other drugs in the panel, then asks whether model predictions still separate resistant from susceptible isolates within those matched co-resistance backgrounds. The audit requires no model retraining, no access to model internals, and no WGS data. It produces raw AUC, matched AUC, background-centered AUC, pairwise within-background accuracy, matched retention, and adequacy flags for each drug-site row.

We apply this framework to an expanded *E. coli* DRIAMS panel and focus the main evidence on a contrast between ciprofloxacin and amoxicillin-clavulanic acid. Ciprofloxacin resistance in *E. coli* is strongly linked to globally disseminated lineages, especially ST131, whereas amoxicillin-clavulanic acid resistance is more ecologically heterogeneous \cite{nicolas2014st131}. We then test whether the same audit behavior appears across model families, examine the co-resistance network that models can exploit, and use independent public WGS-linked MALDI data to ask whether lineage-associated spectral shortcuts are biologically plausible.

## Results

### Raw external performance did not identify the source of signal

We first evaluated raw external performance for a DRIAMS-A-trained *E. coli* MALDI-AMR model across a temporal DRIAMS-A 2018 holdout and external DRIAMS-B, DRIAMS-C, and DRIAMS-D test sites. The six-drug panel included ciprofloxacin, norfloxacin, amoxicillin-clavulanic acid, ceftriaxone, ceftazidime, and cefepime. Raw AUCs were heterogeneous across drugs and sites. At A-2018, AUC ranged from 0.650 for amoxicillin-clavulanic acid to 0.867 for cefepime; ciprofloxacin reached 0.823. At DRIAMS-C and DRIAMS-D, ciprofloxacin raw AUCs were 0.750 and 0.671, while amoxicillin-clavulanic acid raw AUCs were 0.535 and 0.557.

These values describe aggregate transfer, but they do not identify what the model used to transfer. In a panel where resistance phenotypes are correlated, a model can assign higher scores to isolates from a resistant background even when it does not rank resistant isolates above susceptible isolates within the same background. This is the central distinction tested by the background-matched audit.

### The audit measured prediction inside matched co-resistance backgrounds

For each focal drug, we constructed a co-resistance signature from all non-focal drug AST labels in the same organism panel. Isolates with the same signature formed a background stratum. A valid stratum required at least three resistant and three susceptible isolates for the focal drug. Within each site-drug row, we computed four performance summaries. Raw AUC used all available isolates. Matched AUC restricted evaluation to valid strata. Background-centered AUC subtracted the mean model score within each valid stratum before pooling residual scores and recomputing AUC. Pairwise within-background accuracy counted only resistant-susceptible pairs drawn from the same valid stratum, with ties counted as 0.5.

The distinction between the last two statistics is important. Background-centered AUC removes stratum-level score shifts but remains a pooled residual AUC. Pairwise within-background accuracy is the direct within-stratum statistic. We therefore report both, together with matched retention and adequacy flags. Low retention is not interpreted as biological failure; it is reported as an audit limitation because the dataset does not contain enough matched resistant-susceptible comparisons to support a strong conclusion.

Figure 1 illustrates this workflow: isolate-level predictions enter the audit; co-resistance signatures define strata; valid strata support matched and centered evaluation; rows with low matched support are flagged rather than overinterpreted.

### Ciprofloxacin retained signal after matching, whereas amoxicillin-clavulanic acid collapsed externally

The primary comparison was *E. coli*/ciprofloxacin versus *E. coli*/amoxicillin-clavulanic acid. The two drugs were evaluated with the same training site, model pipeline, organism, and external sites, but they differed in resistant-population ecology.

Ciprofloxacin retained residual discrimination after background matching at the main interpretable sites. Background-centered AUC was 0.703 at A-2018, 0.646 at DRIAMS-C, and 0.596 at DRIAMS-D. The direct within-background statistic agreed with this pattern: pairwise within-background accuracy was 0.754, 0.675, and 0.616 at the same sites. DRIAMS-B had a high centered estimate, but it was based on only 25 matched isolates and one valid stratum, so it is treated as cautionary support rather than primary evidence.

Amoxicillin-clavulanic acid behaved differently. At A-2018, the model showed only weak residual signal after centering, with centered AUC 0.541. At DRIAMS-C and DRIAMS-D, centered AUC fell to 0.497 and 0.486 despite high matched retention, and pairwise within-background accuracy fell to chance or below chance. DRIAMS-B again provided only uncertain support. Thus, the contrast was not simply that one drug had higher raw AUC than another. The key result is that ciprofloxacin retained ranking information inside matched backgrounds, whereas amoxicillin-clavulanic acid did not retain meaningful external within-background signal.

Figure 2 shows raw versus background-centered AUC for the primary drug-site rows. Table 1 reports raw AUC, centered AUC, matched retention, valid strata, pairwise within-background accuracy, and interpretation labels for the primary contrast.

### Co-resistance structure provided an exploitable background

The audit results were consistent with the structure of the AST panel. Across the *E. coli* panel, ciprofloxacin and norfloxacin formed an exceptionally strong fluoroquinolone block, with phi = 0.976 and resistant Jaccard similarity = 0.963. Extended-spectrum cephalosporins also formed strong co-resistance blocks: ceftriaxone-ceftazidime phi = 0.884, ceftazidime-cefepime phi = 0.828, and ceftriaxone-cefepime phi = 0.804. These correlations mean that a model trained on focal resistance labels is also exposed to a broader label ecology in which resistance to one drug often identifies resistance to another.

This structure is exactly the kind of background that can inflate or destabilize raw AUC. If a model assigns higher scores to isolates from a co-resistant population, it can perform well in aggregate even when focal-drug discrimination within matched backgrounds is weak. The cephalosporin rows illustrate the same caution. Several raw AUCs were high, but many matched analyses had low retention or strong attenuation after centering, so those rows are reported as secondary or cautionary rather than as clean evidence of focal signal.

### Background sensitivity was not specific to one model architecture

To test whether the audit result was a neural-network artifact, we applied the same background-matched procedure to convolutional neural network predictions and multi-task LightGBM predictions. The qualitative pattern was shared. In interpretable amoxicillin-clavulanic acid rows, the mean raw-minus-centered drop was 0.084 for the CNN and 0.116 for LightGBM. In interpretable ciprofloxacin rows, the mean background-centered AUC was 0.648 for the CNN and 0.618 for LightGBM. At DRIAMS-C and DRIAMS-D, amoxicillin-clavulanic acid collapsed after matching in both model families, while ciprofloxacin retained partial signal.

This result does not prove that the two architectures learned identical features. It shows something more important for evaluation: background sensitivity persisted when the model family changed. The confound therefore appears to originate in the relationship among spectra, AST labels, and resistant-population ecology, not in a single representation learned by one architecture.

The full model-family comparison belongs in the supplement. It is important for robustness, but it is secondary to the primary clinical claim.

### Public WGS-linked MALDI data supported lineage as a plausible spectral shortcut

The DRIAMS isolates used in the primary audit do not have linked WGS labels, so the audit cannot directly prove that a specific lineage drove any DRIAMS prediction. We therefore used independent public WGS-linked Basel uropathogenic *E. coli* data to test a narrower biological premise: can MALDI peak features detect lineage strongly enough to act as a shortcut for resistance prediction?

In 407 isolates with Bruker MALDI peak features and WGS-derived metadata, ST131 lineage identity was predicted from MALDI peak features with AUC 0.932. In the same cohort, ciprofloxacin resistance was predicted with AUC 0.755 and ceftriaxone resistance with AUC 0.689. ST131 was also strongly associated with these resistance phenotypes. Exact sequence-type centering attenuated ciprofloxacin AUC to 0.564 and ceftriaxone AUC to 0.503. Thus, MALDI spectra contained a strong lineage signal, and resistance prediction weakened when fine-grained lineage background was controlled.

We also compared discriminative MALDI bins with published ST131 biomarker masses \cite{nakamura2019st131}. ST131-discriminative bins were enriched for published ST131 biomarker neighborhoods by 3.11-fold. Ciprofloxacin- and ceftriaxone-discriminative bins were also enriched, at 2.24-fold and 2.64-fold, respectively. These are mass-neighborhood overlaps, not protein identifications. They support the plausibility that resistance-associated MALDI features can overlap lineage-associated features, but they do not establish protein identity or prove the same mechanism inside DRIAMS.

Figure 3 combines the strongest biological support: the co-resistance structure available to the models and the public WGS-linked MALDI evidence that ST131 is strongly encoded in spectra. If the figure becomes crowded in journal layout, the full cross-resistance heatmap can move to supplement and the main figure can focus on the WGS-linked MALDI panel.

### Stress tests defined the boundary of the claim

We ran two additional analyses to define what the framework can and cannot claim.

First, falsification controls compared observed model AUC with a background-burden-only score and a within-background shuffled-label null. None of 23 audited pair-site rows exceeded both controls. Three rows exceeded the shuffle null while remaining competitive with the burden-only control. These controls make the interpretation conservative: high raw AUC should not be treated as focal-drug evidence unless it survives background-aware comparisons.

Second, we aggregated MARISMa spot-level predictions to isolate-drug rows and audited the resulting external predictions. The locked DRIAMS-trained model showed weak raw signal across all six targets in this stress test, with raw AUCs near chance. This is a transfer failure, not a successful external validation. It is useful because it shows how the audit should be used in practice: when raw performance is weak, the correct conclusion is failed transfer before any more nuanced background interpretation.

These analyses should remain supplementary. They strengthen the paper by preventing overclaiming, but the main story should stay focused on the DRIAMS primary audit and the public WGS-linked biological support.

## Discussion

This study shows that MALDI-TOF AMR prediction should not be interpreted from raw external AUC alone. In the DRIAMS *E. coli* panel, ciprofloxacin predictions retained measurable within-background signal after matching by co-resistance context, whereas amoxicillin-clavulanic acid weakened or collapsed at external sites. The same qualitative pattern appeared across model families. Cross-resistance networks showed strong label ecology, and independent WGS-linked MALDI data showed that lineage can be predicted from routine MALDI peak features more strongly than resistance in the same cohort. Together, these results support a practical conclusion: MALDI-TOF AMR models are background-sensitive predictors, and their evaluation should report how much performance survives resistant-population background control.

The main contribution is not the observation that bacterial population structure exists. That is well established in microbial genomics and AMR prediction \cite{yu2025population}. The contribution is an operational audit that can be applied to any fixed MALDI-TOF AMR prediction table. The audit converts a vague concern, "maybe the model learned the resistant clone," into row-level quantities: raw AUC, matched AUC, background-centered AUC, pairwise within-background accuracy, matched retention, and adequacy labels. This makes the confound visible without requiring WGS labels or model internals.

The ciprofloxacin result should be interpreted carefully. Retained within-background signal does not prove that the model learned a biochemical fluoroquinolone resistance mechanism. Ciprofloxacin resistance is linked to ST131 and other fluoroquinolone-resistant *E. coli* lineages, and MALDI can encode lineage-associated protein variation. The audit only shows that after matching on the observed co-resistance panel, some ranking signal remains. That remaining signal may reflect focal biology, residual lineage not captured by the AST signature, or both. This is why the paper reports pairwise within-background accuracy and matched retention instead of presenting centered AUC as proof of mechanism.

The amoxicillin-clavulanic acid result is more directly actionable. At external sites with high matched retention, performance approached chance after background adjustment. A model with that behavior should not be interpreted as a robust focal-drug predictor simply because raw AUC is above chance in some sites. For clinical microbiology, this is the point of the audit: it identifies when apparent transfer is likely to depend on the distribution of resistant populations in the test site.

The public WGS-linked UPEC analysis provides biological grounding but not direct DRIAMS lineage proof. It shows that MALDI peak features can identify ST131 very strongly and that resistance predictions are attenuated when sequence-type background is controlled. Because those isolates are independent of the DRIAMS benchmark, this analysis should be read as mechanistic plausibility and external support, not as evidence that the DRIAMS model detected ST131 in the primary data.

Several limitations are important. First, co-resistance background is only a proxy for population background. Isolates with identical AST signatures can differ by sequence type, plasmid content, and resistance mechanism. DRIAMS-linked WGS or MLST would allow direct lineage-centered decomposition. Second, background-centered AUC is a pooled residual statistic, not a purely within-stratum statistic; that is why pairwise within-background accuracy is reported as the strict companion metric. Third, sparse matched strata limit interpretation for some drug-site rows, especially cephalosporins and smaller site panels. Fourth, mass-bin enrichment against published ST131 biomarker masses does not identify proteins. Fifth, the Weis/Borgwardt analysis supports official LR parity only for the limited upstream panel already checked, while the broader six-drug export is a background-audit stress test rather than an exact replication claim. Sixth, the present primary analysis focuses on *E. coli*. Applying the same framework to *Klebsiella pneumoniae*, *Staphylococcus aureus*, and other organisms is a necessary next step.

For the field, the implication is straightforward. MALDI-TOF AMR studies should report background-centered AUC, pairwise within-background accuracy, matched retention, and a cross-resistance summary alongside raw AUC. A model should be considered deployment-ready for a drug-site context only when its predictions survive the observable resistant-population background available in that context. This standard would not eliminate the need for prospective clinical validation, but it would make preclinical and retrospective validation more honest about what the model may have learned.

## Methods

### Primary DRIAMS data and evaluation design

We used the DRIAMS MALDI-TOF and AST resource as the primary clinical benchmark \cite{weis2022driams,driamsdryad}. Models were trained on DRIAMS-A isolates collected before 2018 and evaluated on a temporal DRIAMS-A 2018 holdout and the external DRIAMS-B, DRIAMS-C, and DRIAMS-D sites. The primary organism panel was *E. coli*. The six-drug panel included ciprofloxacin, norfloxacin, amoxicillin-clavulanic acid, ceftriaxone, ceftazidime, and cefepime. Intermediate AST labels and isolates missing the focal-drug label were excluded from the corresponding focal-drug evaluation.

### Prediction models

The primary model was a multi-task convolutional neural network trained on MALDI-TOF spectra represented over the 2,000 to 20,000 m/z range. A multi-task LightGBM model trained on peak-bin features from the same general spectral range was used as an architecture check. The LightGBM analysis was not intended to optimize performance; it tested whether background sensitivity appeared outside the neural-network model. Weis/Borgwardt-style logistic regression outputs were exported separately with an ID-preserving pipeline and audited as supplementary published-workflow compatibility evidence.

### Background signatures

For each isolate and focal drug, the co-resistance background signature was the concatenation of the binary AST labels for all non-focal drugs in the same organism panel. For ciprofloxacin, for example, the signature encoded norfloxacin, amoxicillin-clavulanic acid, ceftriaxone, ceftazidime, and cefepime labels. The focal drug was never included in its own background signature. A valid stratum required at least three resistant and three susceptible isolates for the focal drug. Strata failing this threshold were excluded from matched and centered performance estimates and counted against matched retention.

### Audit metrics

Raw AUC was computed from all available focal-drug prediction rows in a site. Matched AUC was computed after restricting to valid background strata. Background-centered AUC was computed by subtracting the mean model score within each valid stratum from each isolate score, then recomputing AUC from the pooled residual scores. Pairwise within-background accuracy was computed from resistant-susceptible pairs within the same valid stratum and measured the fraction of pairs in which the resistant isolate received the higher model score; ties counted as 0.5. Matched retention was the fraction of focal-drug rows retained in valid strata.

### Uncertainty and adequacy labels

Bootstrap confidence intervals were computed with 500 resamples. Permutation p-values for the raw-to-centered drop were computed by shuffling focal-drug labels within valid background strata and recomputing centered AUC over 500 permutations. Rows with no valid strata were labeled not interpretable. Rows with very small matched sample counts or low matched retention were labeled cautionary and were not used as primary evidence even when centered AUC was high.

### Cross-resistance network

For each pair of drugs in the *E. coli* panel, we computed phi correlation, resistant-resistant lift, and resistant Jaccard similarity from isolate-level binary AST labels among isolates with both labels known. These measures described the label ecology available to AMR prediction models and identified co-resistance blocks that could serve as background signal.

### Public WGS-linked UPEC support analysis

We analyzed public Basel uropathogenic *E. coli* resources linking WGS metadata and Bruker MALDI-derived peak features \cite{cuenod2023papgii}. The analysis included 407 isolates with complete usable peak features and metadata after manifest joining; 55 were ST131 by WGS-derived sequence type. We created 50 Da-binned MALDI peak features and evaluated centroid-direction classifiers for ST131 identity, ciprofloxacin resistance, and ceftriaxone resistance using five cross-validation folds. To estimate lineage-controlled resistance signal, we computed centered AUCs after subtracting WGS-derived background group means, including exact sequence type. Sequence-type-centered analyses retained only groups containing both resistant and susceptible isolates for the focal phenotype.

### Published ST131 biomarker enrichment

Published ST131 MALDI biomarker masses were taken from Nakamura et al. \cite{nakamura2019st131}. For each target, we ranked MALDI bins by absolute mean log-intensity difference between positive and negative isolates and selected the top 75 bins. We counted overlaps within 40 Da of published ST131 biomarker masses and compared observed overlap with a mass-stratified permutation null. Individual mass-bin matches were interpreted as spectral-neighborhood overlaps only. They were not treated as protein identifications.

### Falsification controls

For each audited drug-site row, we compared observed model AUC with a background-burden-only score, defined as the number of resistant non-focal drug labels in the isolate's panel. We also generated a within-background shuffled-label null by permuting focal-drug labels within valid background strata. These controls were used to assess whether observed raw performance exceeded simple background-derived alternatives.

### MARISMa stress test

MARISMa predictions were generated at spot level and aggregated before auditing to avoid treating replicate spectra as independent isolates. Rows were grouped by site, year, isolate identifier, organism, and drug. The isolate-drug probability was the mean of replicate spot probabilities. Groups with conflicting binary labels were excluded and written to a duplicate-handling report. The committed stress-test snapshot collapsed 9,597 spot-level rows to 6,022 isolate-drug rows and excluded 20 conflicting isolate-drug groups before auditing.

### Software and reproducibility

All analyses operate on a flat prediction table with isolate identifier, site, year, organism, drug, binary AST label, and model probability. The main audit engine is `run_background_audit_framework.py`. Publication outputs are generated by `scripts/make_final_framework_tables_figures.py` and related figure scripts. Weis/Borgwardt published-code prediction exports are produced by `scripts/export_weis_predictions_for_audit.py` and parity checks by `scripts/compare_weis_raw_metrics.py`.

## Data Availability

Raw DRIAMS spectra and AST labels are available from the DRIAMS Dryad release \cite{driamsdryad}. Public WGS-linked UPEC resources are available through the Cuenod et al. study resources \cite{cuenod2023papgii}. Derived summary tables and figure source data are included in this repository under `outputs/final_framework_outputs/`, `outputs/analysis_outputs/`, and `manuscript/source_data/`.

## Code Availability

Analysis code is available at `https://github.com/byungkim113/background-matched-maldi-amr-audit`, including the model-agnostic audit engine, prediction exporters, public WGS-linked MALDI support analyses, source-data tables, and figure-generation scripts.

## Acknowledgments

We thank the DRIAMS investigators and the public UPEC WGS/MALDI data generators for releasing resources that make independent evaluation possible.

## Author Contributions

B.K. conceived the audit framework, implemented the analyses, generated outputs, and drafted the manuscript. Additional author contributions should be completed and approved before submission.

## Competing Interests

The authors declare no competing interests. This statement should be confirmed by all authors before submission.

## Main Display Items

### Figure 1. Background-matched audit workflow

The audit accepts isolate-level prediction rows from any MALDI-TOF AMR model. For each focal drug, non-focal AST labels define a co-resistance background signature. Valid strata contain at least three resistant and three susceptible isolates for the focal drug. The audit reports raw AUC, matched AUC, background-centered AUC, pairwise within-background accuracy, matched retention, and adequacy labels.

Recommended file: `manuscript/figures/figure_1_framework.pdf`.

### Figure 2. Primary DRIAMS background-matched audit

Raw and background-centered AUCs for *E. coli*/ciprofloxacin and *E. coli*/amoxicillin-clavulanic acid across temporal and external DRIAMS sites. Ciprofloxacin retains residual within-background discrimination at the main interpretable sites, whereas amoxicillin-clavulanic acid weakens or collapses toward chance externally. Rows with low matched support are flagged rather than treated as primary evidence.

Recommended file: `manuscript/figures/figure_2_primary_background_audit.pdf`.

### Figure 3. Resistant-population background is spectrally and phenotypically structured

Combined display showing the strongest co-resistance blocks in the *E. coli* AST panel and public WGS-linked MALDI evidence. The AST network shows strong fluoroquinolone and cephalosporin blocks. The public UPEC analysis shows ST131 is highly predictable from MALDI peak features and that resistance-associated discriminative bins are enriched for published ST131 biomarker mass neighborhoods.

Recommended source files: `manuscript/figures/figure_4_cross_resistance_network.pdf` and `manuscript/figures/figure_5_public_wgs_proteomic_support.pdf`. If this is too crowded, use the WGS-linked support figure as the main display and move the full cross-resistance heatmap to supplement.

### Table 1. Primary background-matched audit rows

Shortened table reporting raw AUC, background-centered AUC, pairwise within-background accuracy, matched retention, valid strata, adequacy label, and interpretation for the primary ciprofloxacin and amoxicillin-clavulanic acid rows. The full drug-site table should be Supplementary Table 1.

Recommended source file: `outputs/final_framework_outputs/table_1_primary_background_matched_audit.csv`.

## Supplementary Display Items

Supplementary Figure 1: CNN versus LightGBM model-family audit.

Supplementary Figure 2: Full cross-resistance network or heatmap if not used in main Figure 3.

Supplementary Figure 3: Falsification controls against background-burden-only and within-background shuffled-label baselines.

Supplementary Figure 4: Weis/Borgwardt official LR parity and six-drug E. coli background-audit stress test.

Supplementary Figure 5: MARISMa isolate-level external stress test.

Supplementary Figure 6: Deployment decision flow, if not included as a main translational schematic.

Supplementary Tables: Full audit rows, model-family comparison, co-resistance stratification, public WGS-linked AUCs, biomarker enrichment, calibration, temporal reliability, falsification controls, Weis audit, MARISMa duplicate-handling and stress-test rows, and deployment readiness rules.
