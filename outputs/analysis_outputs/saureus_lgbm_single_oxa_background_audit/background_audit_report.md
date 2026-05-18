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
| A-2018 | Staphylococcus aureus | Ciprofloxacin | 0.765865 | 0.726317 | 0.857143 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Clindamycin | 0.588432 | 0.502709 | 0.749126 | interpretable | weak_raw_signal |
| A-2018 | Staphylococcus aureus | Erythromycin | 0.607152 | 0.46552 | 0.834495 | interpretable | background_driven_collapse |
| A-2018 | Staphylococcus aureus | Fusidic acid | 0.621753 | 0.613102 | 0.755738 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Gentamicin | 0.643705 | 0.501807 | 0.517597 | interpretable | background_driven_collapse |
| A-2018 | Staphylococcus aureus | Oxacillin | 0.839466 | 0.730665 | 0.600816 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Penicillin | 0.729213 | 0.72194 | 0.777778 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Ciprofloxacin | 0.627773 | 0.710227 | 0.532738 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Clindamycin | 0.612584 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Erythromycin | 0.602793 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Fusidic acid | 0.558436 | 0.544262 | 0.557803 | interpretable | weak_raw_signal |
| DRIAMS-B | Staphylococcus aureus | Gentamicin | 0.536765 | 0.347985 | 0.270115 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-B | Staphylococcus aureus | Oxacillin | 0.824762 | 0.871212 | 0.526012 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Penicillin | 0.658536 | 0.665789 | 0.83908 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Ciprofloxacin | 0.676682 | 0.705891 | 0.825203 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Clindamycin | 0.51611 | 0.484 | 0.930876 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Fusidic acid | 0.596307 | 0.629357 | 0.814132 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Gentamicin | 0.48935 | 0.518856 | 0.607046 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Oxacillin | 0.71764 | 0.733367 | 0.596206 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Penicillin | 0.705811 | 0.689898 | 0.915989 | interpretable | focal_signal_retained |
| DRIAMS-D | Staphylococcus aureus | Ciprofloxacin | 0.695527 | 0.662356 | 0.844054 | interpretable | focal_signal_retained |
| DRIAMS-D | Staphylococcus aureus | Clindamycin | 0.527941 | 0.338462 | 0.166667 | caution_low_n_matched_and_low_pairwise | caution_low_matched_support |
| DRIAMS-D | Staphylococcus aureus | Erythromycin | 0.598503 | 0.578011 | 0.959078 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Gentamicin | 0.578081 | 0.485674 | 0.547179 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Penicillin | 0.72553 | 0.716922 | 0.940523 | interpretable | focal_signal_retained |
