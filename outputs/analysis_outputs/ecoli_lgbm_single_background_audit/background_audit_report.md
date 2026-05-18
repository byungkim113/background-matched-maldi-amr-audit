# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Ciprofloxacin | Norfloxacin | 0.975562 | 3.99267 | 263 | 1090 |
| Ceftriaxone | Ceftazidime | 0.883538 | 7.84082 | 446 | 4349 |
| Ceftazidime | Cefepime | 0.827912 | 10.2122 | 271 | 3907 |
| Ceftriaxone | Cefepime | 0.803566 | 9.39376 | 283 | 3956 |
| Norfloxacin | Cefepime | 0.479034 | 2.9082 | 89 | 771 |
| Ciprofloxacin | Ceftriaxone | 0.478556 | 3.22292 | 437 | 4371 |
| Norfloxacin | Ceftriaxone | 0.41936 | 2.55865 | 123 | 1079 |
| Ciprofloxacin | Cefepime | 0.414441 | 3.98808 | 222 | 3873 |
| Ciprofloxacin | Ceftazidime | 0.402416 | 3.26275 | 310 | 4246 |
| Amoxicillin-Clavulanic acid | Ceftazidime | 0.378758 | 3.01416 | 318 | 4304 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| A-2018 | Escherichia coli | Amoxicillin-Clavulanic acid | 0.613805 | 0.465149 | 0.922458 | interpretable | background_driven_collapse |
| A-2018 | Escherichia coli | Cefepime | 0.887818 | 0.575735 | 0.0601626 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| A-2018 | Escherichia coli | Ceftazidime | 0.804485 | 0.410417 | 0.0373906 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| A-2018 | Escherichia coli | Ceftriaxone | 0.861078 | 0.825749 | 0.109353 | interpretable | focal_signal_retained |
| A-2018 | Escherichia coli | Ciprofloxacin | 0.827574 | 0.652469 | 0.619259 | interpretable | focal_signal_retained |
| A-2018 | Escherichia coli | Norfloxacin | 0.684953 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Amoxicillin-Clavulanic acid | 0.603686 | 0.451507 | 0.621212 | interpretable | background_driven_collapse |
| DRIAMS-B | Escherichia coli | Cefepime | 0.798906 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ceftazidime | 0.817063 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ceftriaxone | 0.835317 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-B | Escherichia coli | Ciprofloxacin | 0.815159 | 0.8 | 0.117925 | caution_low_n_matched | caution_low_matched_support |
| DRIAMS-B | Escherichia coli | Norfloxacin | 0.753006 | 0.436697 | 0.542857 | interpretable | background_driven_collapse |
| DRIAMS-C | Escherichia coli | Amoxicillin-Clavulanic acid | 0.559045 | 0.532529 | 0.892461 | interpretable | weak_raw_signal |
| DRIAMS-C | Escherichia coli | Cefepime | 0.564802 | 0.272727 | 0.0771757 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| DRIAMS-C | Escherichia coli | Ceftazidime | 0.605888 | 0.292308 | 0.0748899 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| DRIAMS-C | Escherichia coli | Ceftriaxone | 0.649929 | 0.768519 | 0.0361446 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| DRIAMS-C | Escherichia coli | Ciprofloxacin | 0.804403 | 0.682482 | 0.476298 | interpretable | focal_signal_retained |
| DRIAMS-C | Escherichia coli | Norfloxacin | 0.785991 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| DRIAMS-D | Escherichia coli | Amoxicillin-Clavulanic acid | 0.559239 | 0.486938 | 0.974423 | interpretable | weak_raw_signal |
| DRIAMS-D | Escherichia coli | Cefepime | 0.749695 | 0.476339 | 0.049402 | caution_low_n_matched_and_low_retention | caution_low_matched_support |
| DRIAMS-D | Escherichia coli | Ceftazidime | 0.703359 | 0.472323 | 0.83392 | interpretable | background_driven_collapse |
| DRIAMS-D | Escherichia coli | Ceftriaxone | 0.730968 | 0.544648 | 0.821966 | interpretable | background_driven_collapse |
| DRIAMS-D | Escherichia coli | Ciprofloxacin | 0.712116 | 0.591114 | 0.984012 | interpretable | partially_retained_or_uncertain |
