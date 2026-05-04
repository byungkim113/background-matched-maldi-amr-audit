# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Ceftriaxone | Ceftazidime | 1 | 5.06494 | 231 | 1170 |
| Ceftriaxone | Cefepime | 0.961792 | 5.0917 | 215 | 1166 |
| Ceftazidime | Cefepime | 0.961792 | 5.0917 | 215 | 1166 |
| Ciprofloxacin | Norfloxacin | 0.93268 | 2.31304 | 106 | 266 |
| Ciprofloxacin | Ceftriaxone | 0.536979 | 2.44153 | 190 | 1101 |
| Ciprofloxacin | Ceftazidime | 0.536979 | 2.44153 | 190 | 1101 |
| Ciprofloxacin | Cefepime | 0.52142 | 2.45702 | 179 | 1098 |
| Norfloxacin | Ceftriaxone | 0.454334 | 2.13496 | 44 | 279 |
| Norfloxacin | Ceftazidime | 0.454334 | 2.13496 | 44 | 279 |
| Norfloxacin | Cefepime | 0.447555 | 2.1328 | 43 | 279 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| MARISMa | Escherichia coli | Amoxicillin-Clavulanic acid | 0.518595 | 0.537513 | 0.963745 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Cefepime | 0.472523 | 0.415942 | 0.204787 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Ceftazidime | 0.479139 | 0.562264 | 0.143541 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Ceftriaxone | 0.477515 | 0.535849 | 0.143541 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Ciprofloxacin | 0.518688 | 0.506902 | 0.862422 | interpretable | weak_raw_signal |
| MARISMa | Escherichia coli | Norfloxacin | 0.541468 | 0.474696 | 0.502874 | interpretable | weak_raw_signal |
