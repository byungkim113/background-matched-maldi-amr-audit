# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Ceftriaxone | Ceftazidime | 1 | 5.15487 | 226 | 1165 |
| Ceftriaxone | Cefepime | 0.961092 | 5.18304 | 210 | 1161 |
| Ceftazidime | Cefepime | 0.961092 | 5.18304 | 210 | 1161 |
| Ciprofloxacin | Norfloxacin | 0.93268 | 2.31304 | 106 | 266 |
| Ciprofloxacin | Ceftriaxone | 0.535618 | 2.46143 | 186 | 1096 |
| Ciprofloxacin | Ceftazidime | 0.535618 | 2.46143 | 186 | 1096 |
| Ciprofloxacin | Cefepime | 0.519869 | 2.47766 | 175 | 1093 |
| Norfloxacin | Ceftriaxone | 0.454334 | 2.13496 | 44 | 279 |
| Norfloxacin | Ceftazidime | 0.454334 | 2.13496 | 44 | 279 |
| Norfloxacin | Cefepime | 0.447555 | 2.1328 | 43 | 279 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| MARISMa | Escherichia coli | Amoxicillin-Clavulanic acid | 0.50355 | 0.526842 | 0.961739 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Cefepime | 0.465519 | 0.430702 | 0.0920826 | caution_low_retention | caution_low_matched_support |
| MARISMa | Escherichia coli | Ceftazidime | 0.468164 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| MARISMa | Escherichia coli | Ceftriaxone | 0.46709 |  | 0 | not_interpretable_no_valid_strata | insufficient_matched_overlap |
| MARISMa | Escherichia coli | Ciprofloxacin | 0.54469 | 0.537698 | 0.789282 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Norfloxacin | 0.54161 | 0.458753 | 0.55914 | interpretable | weak_raw_signal |
