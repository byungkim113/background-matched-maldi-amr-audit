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
| A-2018 | Staphylococcus aureus | Ciprofloxacin | 0.721608 | 0.670416 | 0.857143 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Clindamycin | 0.568608 | 0.50024 | 0.749126 | interpretable | weak_raw_signal |
| A-2018 | Staphylococcus aureus | Erythromycin | 0.591863 | 0.490138 | 0.834495 | interpretable | weak_raw_signal |
| A-2018 | Staphylococcus aureus | Fusidic acid | 0.570416 | 0.521387 | 0.755738 | interpretable | weak_raw_signal |
| A-2018 | Staphylococcus aureus | Gentamicin | 0.689534 | 0.504216 | 0.517597 | interpretable | background_driven_collapse |
| A-2018 | Staphylococcus aureus | Oxacillin | 0.794721 | 0.690489 | 0.600816 | interpretable | focal_signal_retained |
| A-2018 | Staphylococcus aureus | Penicillin | 0.640891 | 0.596414 | 0.777778 | interpretable | partially_retained_or_uncertain |
| DRIAMS-B | Staphylococcus aureus | Ciprofloxacin | 0.683673 | 0.638258 | 0.532738 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Clindamycin | 0.564178 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Erythromycin | 0.577585 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Staphylococcus aureus | Fusidic acid | 0.66319 | 0.59235 | 0.557803 | interpretable | partially_retained_or_uncertain |
| DRIAMS-B | Staphylococcus aureus | Gentamicin | 0.61875 | 0.501832 | 0.270115 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-B | Staphylococcus aureus | Oxacillin | 0.774799 | 0.811553 | 0.526012 | interpretable | focal_signal_retained |
| DRIAMS-B | Staphylococcus aureus | Penicillin | 0.554957 | 0.562693 | 0.83908 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Ciprofloxacin | 0.630045 | 0.639709 | 0.825203 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Clindamycin | 0.542777 | 0.500144 | 0.930876 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Fusidic acid | 0.474021 | 0.418112 | 0.814132 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Gentamicin | 0.556662 | 0.409517 | 0.607046 | interpretable | weak_raw_signal |
| DRIAMS-C | Staphylococcus aureus | Oxacillin | 0.655772 | 0.622116 | 0.596206 | interpretable | focal_signal_retained |
| DRIAMS-C | Staphylococcus aureus | Penicillin | 0.584588 | 0.566544 | 0.915989 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Ciprofloxacin | 0.758416 | 0.721268 | 0.844054 | interpretable | focal_signal_retained |
| DRIAMS-D | Staphylococcus aureus | Clindamycin | 0.695956 | 0.615385 | 0.166667 | caution_low_n_matched_and_low_pairwise | caution_low_matched_support |
| DRIAMS-D | Staphylococcus aureus | Erythromycin | 0.572875 | 0.529463 | 0.959078 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Gentamicin | 0.53284 | 0.487614 | 0.547179 | interpretable | weak_raw_signal |
| DRIAMS-D | Staphylococcus aureus | Penicillin | 0.602131 | 0.596248 | 0.940523 | interpretable | partially_retained_or_uncertain |
