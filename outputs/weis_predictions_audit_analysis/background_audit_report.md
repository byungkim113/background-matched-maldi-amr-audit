# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Ceftriaxone | Ceftazidime | 0.855315 | 5.84706 | 28 | 213 |
| Ceftriaxone | Cefepime | 0.770991 | 6.89474 | 12 | 131 |
| Ceftazidime | Cefepime | 0.752066 | 6.2125 | 14 | 142 |
| Ciprofloxacin | Ceftazidime | 0.572652 | 4.98039 | 10 | 127 |
| Amoxicillin-Clavulanic acid | Ceftazidime | 0.547357 | 4.5422 | 12 | 148 |
| Ciprofloxacin | Ceftriaxone | 0.515491 | 4.125 | 9 | 110 |
| Amoxicillin-Clavulanic acid | Ceftriaxone | 0.454284 | 3.37857 | 11 | 129 |
| Amoxicillin-Clavulanic acid | Cefepime | 0.432075 | 3.16667 | 8 | 95 |
| Ciprofloxacin | Cefepime | 0.396931 | 3.27315 | 7 | 101 |
| Norfloxacin | Ceftazidime | 0.299183 | 1.61538 | 6 | 35 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| DRIAMS-B | Escherichia coli | Amoxicillin-Clavulanic acid | 0.461905 | 0.485507 | 0.674419 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-B | Escherichia coli | Cefepime | 0.676471 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ceftazidime | 0.434641 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ceftriaxone | 0.683007 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ciprofloxacin | 0.75 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Norfloxacin | 0.666667 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-C | Escherichia coli | Amoxicillin-Clavulanic acid | 0.626248 | 0.578182 | 0.745856 | interpretable | partially_retained_or_uncertain |
| DRIAMS-C | Escherichia coli | Cefepime | 0.4748 | 0.2 | 0.072 | caution_low_n_matched_and_low_retention_and_low_pairwise | caution_low_matched_support |
| DRIAMS-C | Escherichia coli | Ceftazidime | 0.533139 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-C | Escherichia coli | Ceftriaxone | 0.644589 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-C | Escherichia coli | Ciprofloxacin | 0.625 | 0.816176 | 0.438202 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-C | Escherichia coli | Norfloxacin | 0.694154 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-D | Escherichia coli | Amoxicillin-Clavulanic acid | 0.652425 | 0.524015 | 0.924812 | interpretable | background_driven_collapse |
| DRIAMS-D | Escherichia coli | Cefepime | 0.622082 | 0.5625 | 0.02 | caution_low_n_matched_and_low_retention_and_low_pairwise | caution_low_matched_support |
| DRIAMS-D | Escherichia coli | Ceftazidime | 0.653278 | 0.388889 | 0.132832 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-D | Escherichia coli | Ceftriaxone | 0.627507 | 0.612179 | 0.117794 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-D | Escherichia coli | Ciprofloxacin | 0.667111 | 0.619607 | 0.914948 | interpretable | focal_signal_retained |
