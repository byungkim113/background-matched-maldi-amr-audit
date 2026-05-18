# Background-Matched Transfer Audit Report

## What This Tests

The audit asks whether focal-drug prediction survives after matching isolates by co-resistance background.
Raw AUC can be high because a model learned resistant-population background; stratum-centered AUC is the stricter within-background test.

## Strongest Co-Resistance Edges

| Drug A | Drug B | Phi | Lift | n RR | n |
| --- | --- | ---: | ---: | ---: | ---: |
| Ceftriaxone | Amoxicillin-Clavulanic acid | 0.748201 | 4.04707 | 47 | 300 |
| Ciprofloxacin | Norfloxacin | 0.744667 | 2.89066 | 69 | 300 |
| Ciprofloxacin | Ceftriaxone | 0.161308 | 1.51192 | 26 | 300 |
| Norfloxacin | Ceftriaxone | 0.159728 | 1.44439 | 30 | 300 |
| Norfloxacin | Amoxicillin-Clavulanic acid | 0.0929152 | 1.30273 | 21 | 300 |
| Ciprofloxacin | Amoxicillin-Clavulanic acid | 0.0131714 | 1.04895 | 14 | 300 |

## Audit Summary

| Site | Organism | Drug | Raw AUC | Centered AUC | Retention | Adequacy | Interpretation |
| --- | --- | --- | ---: | ---: | ---: | --- | --- |
| Hospital-A | Escherichia coli | Amoxicillin-Clavulanic acid | 0.519511 | 0.405375 | 0.753333 | interpretable | weak_raw_signal |
| Hospital-A | Escherichia coli | Ceftriaxone | 0.623414 | 0.502457 | 0.886667 | interpretable | background_driven_collapse |
| Hospital-A | Escherichia coli | Ciprofloxacin | 1 | 1 | 0.253333 | caution_low_n_matched | caution_low_matched_support |
| Hospital-A | Escherichia coli | Norfloxacin | 0.998928 | 0.99867 | 0.68 | interpretable | focal_signal_retained |
| Hospital-B | Escherichia coli | Amoxicillin-Clavulanic acid | 0.565867 | 0.416149 | 0.2 | caution_low_n_matched_and_low_pairwise | caution_low_matched_support |
| Hospital-B | Escherichia coli | Ceftriaxone | 0.521171 | 0.425685 | 0.706667 | interpretable | weak_raw_signal |
| Hospital-B | Escherichia coli | Ciprofloxacin | 1 | 0.97547 | 0.833333 | interpretable | focal_signal_retained |
| Hospital-B | Escherichia coli | Norfloxacin | 0.998384 | 0.895544 | 0.88 | interpretable | focal_signal_retained |
