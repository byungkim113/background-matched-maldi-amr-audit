# E. coli Mechanism6 Expanded-Panel Assessment

## Pair-Level Prediction Assessment

| organism         | drug                        | block                               | locked_prediction | observed_external_category | prediction_assessment    | auc_A-2018 | auc_DRIAMS-B | auc_DRIAMS-C | auc_DRIAMS-D | auc_external_mean | auc_source_to_external_delta | aupr_A-2018 | aupr_DRIAMS-B | aupr_DRIAMS-C | aupr_DRIAMS-D | aupr_external_mean | aupr_source_to_external_delta | threshold_transport_flag    | prediction_reason                                                            |
| ---------------- | --------------------------- | ----------------------------------- | ----------------- | -------------------------- | ------------------------ | ---------- | ------------ | ------------ | ------------ | ----------------- | ---------------------------- | ----------- | ------------- | ------------- | ------------- | ------------------ | ----------------------------- | --------------------------- | ---------------------------------------------------------------------------- |
| Escherichia coli | Cefepime                    | ESBL/AmpC cephalosporin block       | weaker_transfer   | transfer                   | contradicted_by_mean_auc | 0.8738     | 0.8146       | 0.6696       | 0.7979       | 0.7607            | -0.1131                      | 0.5653      | 0.5319        | 0.3517        | 0.3136        | 0.3991             | -0.1662                       | DRIAMS-B;DRIAMS-C;DRIAMS-D; | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| Escherichia coli | Ceftazidime                 | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.8426     | 0.7709       | 0.6280       | 0.7477       | 0.7155            | -0.1271                      | 0.4304      | 0.5245        | 0.2612        | 0.3679        | 0.3845             | -0.0459                       | DRIAMS-B;DRIAMS-C;DRIAMS-D; | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| Escherichia coli | Ceftriaxone                 | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.8726     | 0.7763       | 0.6437       | 0.7443       | 0.7214            | -0.1512                      | 0.6894      | 0.5423        | 0.3192        | 0.4020        | 0.4212             | -0.2682                       | DRIAMS-B;                   | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| Escherichia coli | Ciprofloxacin               | Fluoroquinolone block               | transfer          | partial_transfer           | partially_confirmed      | 0.8326     | 0.7938       | 0.7601       | 0.6805       | 0.7448            | -0.0878                      | 0.7008      | 0.6540        | 0.5802        | 0.3903        | 0.5415             | -0.1593                       | DRIAMS-B;DRIAMS-C;DRIAMS-D; | fluoroquinolone resistance expected to have more conserved population signal |
| Escherichia coli | Norfloxacin                 | Fluoroquinolone block               | transfer          | transfer                   | confirmed                | 0.7670     | 0.7328       | 0.7820       |              | 0.7574            | -0.0096                      | 0.5384      | 0.6261        | 0.6879        |               | 0.6570             | 0.1186                        | none                        | same fluoroquinolone/co-resistance block as ciprofloxacin                    |
| Escherichia coli | Amoxicillin-Clavulanic acid | Heterogeneous inhibitor beta-lactam | partial_weak      | weak_transfer              | confirmed                | 0.6504     | 0.7236       | 0.5423       | 0.5559       | 0.6073            | -0.0431                      | 0.4425      | 0.5343        | 0.3017        | 0.2549        | 0.3636             | -0.0789                       | none                        | heterogeneous inhibitor beta-lactam resistance expected to be noisy          |

## Block/Site Summary

| site     | block                               | mean_auc | mean_aupr | mean_sens | mean_spec | n_pairs | n    | n_r | source_auc | source_aupr | auc_source_delta | aupr_source_delta |
| -------- | ----------------------------------- | -------- | --------- | --------- | --------- | ------- | ---- | --- | ---------- | ----------- | ---------------- | ----------------- |
| A-2018   | ESBL/AmpC cephalosporin block       | 0.8630   | 0.5617    | 0.6715    | 0.8225    | 3       | 3877 | 443 | 0.8630     | 0.5617      | 0.0000           | 0.0000            |
| A-2018   | Fluoroquinolone block               | 0.7998   | 0.6196    | 0.7899    | 0.5998    | 2       | 1789 | 442 | 0.7998     | 0.6196      | 0.0000           | 0.0000            |
| A-2018   | Heterogeneous inhibitor beta-lactam | 0.6504   | 0.4425    | 0.5144    | 0.6982    | 1       | 1367 | 383 | 0.6504     | 0.4425      | 0.0000           | 0.0000            |
| DRIAMS-B | ESBL/AmpC cephalosporin block       | 0.7873   | 0.5329    | 0.4850    | 0.7898    | 3       | 639  | 133 | 0.8630     | 0.5617      | -0.0757          | -0.0288           |
| DRIAMS-B | Fluoroquinolone block               | 0.7633   | 0.6401    | 0.7902    | 0.5062    | 2       | 422  | 127 | 0.7998     | 0.6196      | -0.0365          | 0.0205            |
| DRIAMS-B | Heterogeneous inhibitor beta-lactam | 0.7236   | 0.5343    | 0.6721    | 0.6496    | 1       | 198  | 61  | 0.6504     | 0.4425      | 0.0732           | 0.0918            |
| DRIAMS-C | ESBL/AmpC cephalosporin block       | 0.6471   | 0.3107    | 0.3711    | 0.7624    | 3       | 2430 | 400 | 0.8630     | 0.5617      | -0.2159          | -0.2510           |
| DRIAMS-C | Fluoroquinolone block               | 0.7711   | 0.6341    | 0.8180    | 0.4781    | 2       | 1328 | 306 | 0.7998     | 0.6196      | -0.0288          | 0.0145            |
| DRIAMS-C | Heterogeneous inhibitor beta-lactam | 0.5423   | 0.3017    | 0.4957    | 0.5774    | 1       | 902  | 230 | 0.6504     | 0.4425      | -0.1081          | -0.1408           |
| DRIAMS-D | ESBL/AmpC cephalosporin block       | 0.7633   | 0.3612    | 0.5267    | 0.7852    | 3       | 5904 | 456 | 0.8630     | 0.5617      | -0.0997          | -0.2005           |
| DRIAMS-D | Fluoroquinolone block               | 0.6805   | 0.3903    | 0.8302    | 0.3450    | 1       | 1939 | 371 | 0.7998     | 0.6196      | -0.1193          | -0.2293           |
| DRIAMS-D | Heterogeneous inhibitor beta-lactam | 0.5559   | 0.2549    | 0.4940    | 0.5836    | 1       | 1994 | 332 | 0.6504     | 0.4425      | -0.0945          | -0.1876           |

## Macro AUC

| site     | macro_auc | n_pairs |
| -------- | --------- | ------- |
| A-2018   | 0.8065    | 6       |
| DRIAMS-B | 0.7687    | 6       |
| DRIAMS-C | 0.6709    | 6       |
| DRIAMS-D | 0.7053    | 6       |

## Random-CV Evaluation Critique

| row_type | model       | organism   | drug         | random_site             | temporal_site | random_auc | temporal_auc | auc_inflation | random_aupr | temporal_aupr | aupr_inflation | n_random | n_r_random |
| -------- | ----------- | ---------- | ------------ | ----------------------- | ------------- | ---------- | ------------ | ------------- | ----------- | ------------- | -------------- | -------- | ---------- |
| macro    | lgbm_multi  | Macro mean | active pairs | DRIAMS-A-random-holdout | A-2018        | 0.9047     | 0.7924       | 0.1123        | 0.7885      | 0.5268        | 0.2617         | 5014     | 974        |
| macro    | lgbm_single | Macro mean | active pairs | DRIAMS-A-random-holdout | A-2018        | 0.8220     | 0.7800       | 0.0420        | 0.6242      | 0.5195        | 0.1047         | 5014     | 974        |

## Interpretation

- The simple prediction that fluoroquinolones would always outperform the ESBL/AmpC cephalosporin block externally was only partially supported.
- DRIAMS-C strongly supports the expected pattern: fluoroquinolone block AUC 0.771 versus ESBL/AmpC block AUC 0.647.
- DRIAMS-B and DRIAMS-D show that the ESBL/AmpC block can transfer by AUC at some target hospitals, so transferability is site-ecology dependent rather than a fixed drug-mechanism property.
- Amoxicillin-Clavulanic acid remains the consistently weak/noisy task across shifted sites.
- Several beta-lactam tasks show strong threshold imbalance despite reasonable AUC, supporting a deployment claim that source-derived thresholds do not transport reliably.

## What This Is Not

This is not the true background-matched contrastive analysis; that requires per-isolate prediction scores plus the full AST matrix, which were not exported in this run folder.
