# Background-Matched MALDI-AMR Audit: Rough Reviewer Draft

## Central thesis

MALDI-TOF AMR models should not be interpreted as generic resistance detectors from raw AUC alone. Their predictions can reflect focal-drug signal, but they can also reflect the spectral background of resistant subpopulations: lineage, hospital ecology, and co-resistance blocks. The framework asks how much apparent AMR prediction survives after matching isolates by co-resistance background.

## Evidence map

- External transfer audit: shows raw portability varies by organism-drug pair.
- Background-matched contrastive audit: tests whether prediction survives after controlling co-resistance background.
- Cross-resistance network: shows drug labels are organized into exploitable ecology blocks.
- Public WGS/proteomic cross-reference: shows MALDI encodes lineage strongly and resistance peaks overlap published ST131 biomarkers.

## Table 1: table_1_primary_background_matched_audit

Source file: `table_1_primary_background_matched_audit.csv`

| evidence_tier      | pair                  | site     | raw_auc_95ci        | matched_auc_95ci    | stratum_centered_auc_95ci | raw_to_centered_delta | matched_retention_pct |
| ------------------ | --------------------- | -------- | ------------------- | ------------------- | ------------------------- | --------------------- | --------------------- |
| Main contrast      | E. coli / Cipro       | A-2018   | 0.823 (0.799-0.849) | 0.835 (0.801-0.865) | 0.703 (0.664-0.747)       | 0.120                 | 61.900                |
| Cautionary support | E. coli / Cipro       | DRIAMS-B | 0.774 (0.692-0.841) | 0.860 (0.667-1.000) | 0.860 (0.636-1.000)       | -0.086                | 11.800                |
| Main contrast      | E. coli / Cipro       | DRIAMS-C | 0.750 (0.702-0.795) | 0.714 (0.639-0.779) | 0.646 (0.576-0.719)       | 0.104                 | 47.600                |
| Main contrast      | E. coli / Cipro       | DRIAMS-D | 0.671 (0.640-0.700) | 0.667 (0.637-0.696) | 0.596 (0.565-0.624)       | 0.075                 | 98.400                |
| Not main evidence  | E. coli / Norfloxacin | A-2018   | 0.743 (0.678-0.803) |                     |                           |                       | 0.000                 |
| Not main evidence  | E. coli / Norfloxacin | DRIAMS-B | 0.703 (0.631-0.783) | 0.343 (0.167-0.582) | 0.343 (0.176-0.578)       | 0.360                 | 54.300                |
| Not main evidence  | E. coli / Norfloxacin | DRIAMS-C | 0.779 (0.727-0.835) |                     |                           |                       | 0.000                 |
| Main contrast      | E. coli / Amox-Clav   | A-2018   | 0.650 (0.619-0.685) | 0.637 (0.602-0.674) | 0.541 (0.504-0.576)       | 0.109                 | 92.200                |

## Table 2: table_2_cnn_vs_lgbm_multi_background_audit

Source file: `table_2_cnn_vs_lgbm_multi_background_audit.csv`

| site     | drug      | cnn_raw_auc | cnn_centered_auc | cnn_raw_minus_centered | cnn_retention | cnn_adequacy                            | lgbm_raw_auc |
| -------- | --------- | ----------- | ---------------- | ---------------------- | ------------- | --------------------------------------- | ------------ |
| A-2018   | Amox-Clav | 0.650       | 0.541            | 0.109                  | 0.922         | interpretable                           | 0.635        |
| A-2018   | FEP       | 0.867       | 0.659            | 0.208                  | 0.060         | caution_low_n_matched_and_low_retention | 0.886        |
| A-2018   | CAZ       | 0.838       | 0.367            | 0.471                  | 0.037         | caution_low_n_matched_and_low_retention | 0.806        |
| A-2018   | CRO       | 0.863       | 0.645            | 0.218                  | 0.109         | interpretable                           | 0.854        |
| A-2018   | Cipro     | 0.823       | 0.703            | 0.120                  | 0.619         | interpretable                           | 0.828        |
| A-2018   | Norflox   | 0.743       |                  |                        | 0.000         | not_interpretable_no_valid_strata       | 0.744        |
| DRIAMS-B | Amox-Clav | 0.694       | 0.576            | 0.118                  | 0.621         | interpretable                           | 0.749        |
| DRIAMS-B | FEP       | 0.753       |                  |                        | 0.000         | not_interpretable_no_valid_strata       | 0.844        |

## Table 3: table_3_ecology_interpretation

Source file: `table_3_ecology_interpretation.csv`

| site     | drug      | resistance_ecology_block     | strongest_network_partner | partner_phi | partner_lift | cnn_drop | lgbm_drop |
| -------- | --------- | ---------------------------- | ------------------------- | ----------- | ------------ | -------- | --------- |
| A-2018   | Amox-Clav | mixed beta-lactam/background | CAZ                       | 0.379       | 3.014        | 0.109    | 0.153     |
| A-2018   | FEP       | cephalosporin/ESBL block     | CAZ                       | 0.828       | 10.212       | 0.208    | 0.347     |
| A-2018   | CAZ       | cephalosporin/ESBL block     | CRO                       | 0.884       | 7.841        | 0.471    | 0.465     |
| A-2018   | CRO       | cephalosporin/ESBL block     | CAZ                       | 0.884       | 7.841        | 0.218    | 0.129     |
| A-2018   | Cipro     | fluoroquinolone block        | Norflox                   | 0.976       | 3.993        | 0.120    | 0.188     |
| A-2018   | Norflox   | fluoroquinolone block        | Cipro                     | 0.976       | 3.993        |          |           |
| DRIAMS-B | Amox-Clav | mixed beta-lactam/background | CAZ                       | 0.379       | 3.014        | 0.118    | 0.107     |
| DRIAMS-B | FEP       | cephalosporin/ESBL block     | CAZ                       | 0.828       | 10.212       |          |           |

## Table 4: table_4_transfer_prediction_assessment

Source file: `table_4_transfer_prediction_assessment.csv`

| drug      | block                               | locked_prediction | observed_external_category | prediction_assessment    | auc_A-2018 | auc_DRIAMS-B | auc_DRIAMS-C |
| --------- | ----------------------------------- | ----------------- | -------------------------- | ------------------------ | ---------- | ------------ | ------------ |
| CAZ       | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.843      | 0.771        | 0.628        |
| CRO       | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.873      | 0.776        | 0.644        |
| FEP       | ESBL/AmpC cephalosporin block       | weaker_transfer   | transfer                   | contradicted_by_mean_auc | 0.874      | 0.815        | 0.670        |
| Cipro     | Fluoroquinolone block               | transfer          | partial_transfer           | partially_confirmed      | 0.833      | 0.794        | 0.760        |
| Norflox   | Fluoroquinolone block               | transfer          | transfer                   | confirmed                | 0.767      | 0.733        | 0.782        |
| Amox-Clav | Heterogeneous inhibitor beta-lactam | partial_weak      | weak_transfer              | confirmed                | 0.650      | 0.724        | 0.542        |

## Table 5: table_5_top_cross_resistance_edges

Source file: `table_5_top_cross_resistance_edges.csv`

| site | drug_a  | drug_b  | n_both_known | n_rr | rr_lift | phi   | resistant_jaccard |
| ---- | ------- | ------- | ------------ | ---- | ------- | ----- | ----------------- |
| ALL  | Cipro   | Norflox | 1090         | 263  | 3.993   | 0.976 | 0.963             |
| ALL  | CRO     | CAZ     | 4349         | 446  | 7.841   | 0.884 | 0.811             |
| ALL  | CAZ     | FEP     | 3907         | 271  | 10.212  | 0.828 | 0.719             |
| ALL  | CRO     | FEP     | 3956         | 283  | 9.394   | 0.804 | 0.680             |
| ALL  | Norflox | FEP     | 771          | 89   | 2.908   | 0.479 | 0.392             |
| ALL  | Cipro   | CRO     | 4371         | 437  | 3.223   | 0.479 | 0.382             |
| ALL  | Norflox | CRO     | 1079         | 123  | 2.559   | 0.419 | 0.362             |
| ALL  | Cipro   | FEP     | 3873         | 222  | 3.988   | 0.414 | 0.273             |

## Table 6: table_6_ecoli_block_site_summary

Source file: `table_6_ecoli_block_site_summary.csv`

| site     | block                               | mean_auc | mean_aupr | mean_sens | mean_spec | n_pairs | n    |
| -------- | ----------------------------------- | -------- | --------- | --------- | --------- | ------- | ---- |
| A-2018   | ESBL/AmpC cephalosporin block       | 0.863    | 0.562     | 0.671     | 0.823     | 3       | 3877 |
| A-2018   | Fluoroquinolone block               | 0.800    | 0.620     | 0.790     | 0.600     | 2       | 1789 |
| A-2018   | Heterogeneous inhibitor beta-lactam | 0.650    | 0.442     | 0.514     | 0.698     | 1       | 1367 |
| DRIAMS-B | ESBL/AmpC cephalosporin block       | 0.787    | 0.533     | 0.485     | 0.790     | 3       | 639  |
| DRIAMS-B | Fluoroquinolone block               | 0.763    | 0.640     | 0.790     | 0.506     | 2       | 422  |
| DRIAMS-B | Heterogeneous inhibitor beta-lactam | 0.724    | 0.534     | 0.672     | 0.650     | 1       | 198  |
| DRIAMS-C | ESBL/AmpC cephalosporin block       | 0.647    | 0.311     | 0.371     | 0.762     | 3       | 2430 |
| DRIAMS-C | Fluoroquinolone block               | 0.771    | 0.634     | 0.818     | 0.478     | 2       | 1328 |

## Table 7: table_7_public_wgs_maldi_auc

Source file: `table_7_public_wgs_maldi_auc.csv`

| target          | n   | class_0 | class_1 | auc   | folds | model              |
| --------------- | --- | ------- | ------- | ----- | ----- | ------------------ |
| ST131           | 407 | 352     | 55      | 0.932 | 5     | centroid_direction |
| Ciprofloxacin_R | 350 | 297     | 53      | 0.755 | 5     | centroid_direction |
| Ceftriaxone_R   | 360 | 330     | 30      | 0.689 | 5     | centroid_direction |

## Table 8: table_8_public_wgs_st131_resistance_associations

Source file: `table_8_public_wgs_st131_resistance_associations.csv`

| phenotype     | st131_resistant | st131_susceptible | nonst131_resistant | nonst131_susceptible | haldane_odds_ratio |
| ------------- | --------------- | ----------------- | ------------------ | -------------------- | ------------------ |
| Ciprofloxacin | 30              | 17                | 23                 | 280                  | 20.803             |
| Ceftriaxone   | 21              | 27                | 9                  | 303                  | 24.977             |

## Table 9: table_9_published_st131_biomarker_enrichment

Source file: `table_9_published_st131_biomarker_enrichment.csv`

| target          | observed_overlap_count | top_peak_count | null_mean_overlap_count | fold_enrichment | z_score | empirical_p_ge_observed | permutations |
| --------------- | ---------------------- | -------------- | ----------------------- | --------------- | ------- | ----------------------- | ------------ |
| ST131           | 14                     | 75             | 4.506                   | 3.107           | 5.128   | 0.000                   | 10000        |
| Ciprofloxacin_R | 10                     | 75             | 4.474                   | 2.235           | 2.971   | 0.007                   | 10000        |
| Ceftriaxone_R   | 12                     | 75             | 4.539                   | 2.644           | 3.995   | 0.001                   | 10000        |
| ALL_TARGETS     | 36                     | 225            | 13.550                  | 2.657           | 8.739   | 0.000                   | 10000        |

## Table 10: table_10_top_published_st131_peak_overlaps

Source file: `table_10_top_published_st131_peak_overlaps.csv`

| target          | mz_center | published_mz | delta_da | protein              | annotation                                     |
| --------------- | --------- | ------------ | -------- | -------------------- | ---------------------------------------------- |
| Ceftriaxone_R   | 11787.500 | 11783        | 4.500    | cytochrome b562      | soluble cytochrome b562                        |
| Ceftriaxone_R   | 4862.500  | 4857         | 5.500    | HdeA multivalent ion | multivalent ion of HdeA m/z 9710               |
| Ceftriaxone_R   | 8437.500  | 8448         | 10.500   | YnfD                 | uncharacterized protein YnfD                   |
| Ceftriaxone_R   | 6837.500  | 6827         | 10.500   | unidentified         | ST131-specific peak not identified by LC-MS/MS |
| Ceftriaxone_R   | 8337.500  | 8351         | 13.500   | YjbJ                 | UPF0337 protein YjbJ                           |
| Ciprofloxacin_R | 11787.500 | 11783        | 4.500    | cytochrome b562      | soluble cytochrome b562                        |
| Ciprofloxacin_R | 4862.500  | 4857         | 5.500    | HdeA multivalent ion | multivalent ion of HdeA m/z 9710               |
| Ciprofloxacin_R | 7662.500  | 7655         | 7.500    | YahO                 | uncharacterized protein YahO                   |

## Figure 1: figure_1_raw_to_background_centered_auc.png

See `figure_1_raw_to_background_centered_auc.png`.

## Figure 2: figure_2_drop_vs_matched_retention.png

See `figure_2_drop_vs_matched_retention.png`.

## Figure 3: figure_3_cross_resistance_phi_heatmap.png

See `figure_3_cross_resistance_phi_heatmap.png`.

## Figure 4: figure_4_public_wgs_proteomic_support.png

See `figure_4_public_wgs_proteomic_support.png`.

## Figure 5: figure_5_framework_flow.png

See `figure_5_framework_flow.png`.
