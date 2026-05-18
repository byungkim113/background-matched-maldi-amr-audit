# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Erythromycin | Clindamycin | 0.790504 | 5.29668 | 198 | 1550 |
| Oxacillin | Ciprofloxacin | 0.324788 | 4.00021 | 86 | 2255 |
| Oxacillin | Gentamicin | 0.24996 | 5.79254 | 31 | 2249 |
| Oxacillin | Penicillin | 0.218309 | 1.40025 | 236 | 2218 |
| Ciprofloxacin | Erythromycin | 0.209814 | 2.57178 | 111 | 3172 |
| Oxacillin | Erythromycin | 0.184177 | 2.15226 | 61 | 1493 |
| Ciprofloxacin | Clindamycin | 0.180545 | 2.54297 | 61 | 2198 |
| Ciprofloxacin | Gentamicin | 0.141978 | 4.87864 | 26 | 4360 |
| Oxacillin | Clindamycin | 0.139624 | 2.07026 | 59 | 2141 |
| Clindamycin | Gentamicin | 0.120138 | 3.04292 | 20 | 2244 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| A-2018 | Staphylococcus aureus | Ciprofloxacin | 0.760738 | 0.702911 | 0.857143 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Clindamycin | 0.625962 | 0.526846 | 0.749126 | interpretable | background_driven_collapse |
| A-2018 | Staphylococcus aureus | Erythromycin | 0.630668 | 0.539817 | 0.834495 | interpretable | background_driven_collapse |
| A-2018 | Staphylococcus aureus | Fusidic acid | 0.646576 | 0.628557 | 0.755738 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Gentamicin | 0.765831 | 0.659305 | 0.517597 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Oxacillin | 0.838238 | 0.762878 | 0.600816 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Penicillin | 0.626836 | 0.568614 | 0.777778 | interpretable | partially_retained_or_uncertain |
| DRIAMS-B | Staphylococcus aureus | Ciprofloxacin | 0.710736 | 0.702652 | 0.532738 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Clindamycin | 0.585904 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Erythromycin | 0.605093 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Fusidic acid | 0.660429 | 0.569945 | 0.557803 | interpretable | partially_retained_or_uncertain |
| DRIAMS-B | Staphylococcus aureus | Gentamicin | 0.580147 | 0.315018 | 0.270115 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-B | Staphylococcus aureus | Oxacillin | 0.775238 | 0.863636 | 0.526012 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Penicillin | 0.531911 | 0.52678 | 0.83908 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Ciprofloxacin | 0.576531 | 0.590147 | 0.825203 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Clindamycin | 0.529565 | 0.51574 | 0.930876 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Fusidic acid | 0.459733 | 0.504992 | 0.814132 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Gentamicin | 0.488704 | 0.471779 | 0.607046 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Oxacillin | 0.725199 | 0.699471 | 0.596206 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Penicillin | 0.552208 | 0.537408 | 0.915989 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Ciprofloxacin | 0.656042 | 0.63711 | 0.844054 | interpretable | focal_signal_retained |
| DRIAMS-D | Staphylococcus aureus | Clindamycin | 0.483088 | 0.184615 | 0.166667 | caution_low_n_matched_and_low_pairwise | caution_low_matched_support |
| DRIAMS-D | Staphylococcus aureus | Erythromycin | 0.576706 | 0.5686 | 0.959078 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Gentamicin | 0.515826 | 0.441553 | 0.547179 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Penicillin | 0.572733 | 0.571057 | 0.940523 | interpretable | weak_raw_signal |
