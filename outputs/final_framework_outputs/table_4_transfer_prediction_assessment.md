# Locked Transfer Prediction Assessment

Prospective transfer predictions compared with observed external behavior.

| drug      | block                               | locked_prediction | observed_external_category | prediction_assessment    | auc_A-2018 | auc_DRIAMS-B | auc_DRIAMS-C | auc_DRIAMS-D | auc_external_mean | prediction_reason                                                            |
| --------- | ----------------------------------- | ----------------- | -------------------------- | ------------------------ | ---------- | ------------ | ------------ | ------------ | ----------------- | ---------------------------------------------------------------------------- |
| CAZ       | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.843      | 0.771        | 0.628        | 0.748        | 0.716             | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| CRO       | ESBL/AmpC cephalosporin block       | weaker_transfer   | partial_transfer           | partially_confirmed      | 0.873      | 0.776        | 0.644        | 0.744        | 0.721             | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| FEP       | ESBL/AmpC cephalosporin block       | weaker_transfer   | transfer                   | contradicted_by_mean_auc | 0.874      | 0.815        | 0.670        | 0.798        | 0.761             | ESBL/AmpC block expected to be more heterogeneous/site-specific              |
| Cipro     | Fluoroquinolone block               | transfer          | partial_transfer           | partially_confirmed      | 0.833      | 0.794        | 0.760        | 0.680        | 0.745             | fluoroquinolone resistance expected to have more conserved population signal |
| Norflox   | Fluoroquinolone block               | transfer          | transfer                   | confirmed                | 0.767      | 0.733        | 0.782        |              | 0.757             | same fluoroquinolone/co-resistance block as ciprofloxacin                    |
| Amox-Clav | Heterogeneous inhibitor beta-lactam | partial_weak      | weak_transfer              | confirmed                | 0.650      | 0.724        | 0.542        | 0.556        | 0.607             | heterogeneous inhibitor beta-lactam resistance expected to be noisy          |
