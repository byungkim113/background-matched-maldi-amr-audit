# Cross-Resistance Network Summary

Built from the same isolate/drug label matrix used by the background-matched audit.

## Overall Resistance Prevalence

| Drug | Known n | Resistant n | Prevalence |
| --- | ---: | ---: | ---: |
| Ciprofloxacin | 4387 | 972 | 0.222 |
| Norfloxacin | 1091 | 274 | 0.251 |
| Amoxicillin-Clavulanic acid | 4461 | 1006 | 0.226 |
| Ceftriaxone | 4510 | 649 | 0.144 |
| Ceftazidime | 4365 | 482 | 0.110 |
| Cefepime | 3975 | 301 | 0.076 |

## Strongest Positive Phi Edges

| Drug A | Drug B | n | RR observed | RR expected | Lift | Phi | Resistant Jaccard |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Ciprofloxacin | Norfloxacin | 1090 | 0.241 | 0.060 | 3.99 | 0.976 | 0.963 |
| Ceftriaxone | Ceftazidime | 4349 | 0.103 | 0.013 | 7.84 | 0.884 | 0.811 |
| Ceftazidime | Cefepime | 3907 | 0.069 | 0.007 | 10.21 | 0.828 | 0.719 |
| Ceftriaxone | Cefepime | 3956 | 0.072 | 0.008 | 9.39 | 0.804 | 0.680 |
| Norfloxacin | Cefepime | 771 | 0.115 | 0.040 | 2.91 | 0.479 | 0.392 |
| Ciprofloxacin | Ceftriaxone | 4371 | 0.100 | 0.031 | 3.22 | 0.479 | 0.382 |
| Norfloxacin | Ceftriaxone | 1079 | 0.114 | 0.045 | 2.56 | 0.419 | 0.362 |
| Ciprofloxacin | Cefepime | 3873 | 0.057 | 0.014 | 3.99 | 0.414 | 0.273 |

## Strongest Co-Resistance Lift Edges

| Drug A | Drug B | n | Lift | Phi | n RR |
| --- | --- | ---: | ---: | ---: | ---: |
| Ceftazidime | Cefepime | 3907 | 10.21 | 0.828 | 271 |
| Ceftriaxone | Cefepime | 3956 | 9.39 | 0.804 | 283 |
| Ceftriaxone | Ceftazidime | 4349 | 7.84 | 0.884 | 446 |
| Ciprofloxacin | Norfloxacin | 1090 | 3.99 | 0.976 | 263 |
| Ciprofloxacin | Cefepime | 3873 | 3.99 | 0.414 | 222 |
| Ciprofloxacin | Ceftazidime | 4246 | 3.26 | 0.402 | 310 |
| Ciprofloxacin | Ceftriaxone | 4371 | 3.22 | 0.479 | 437 |
| Amoxicillin-Clavulanic acid | Cefepime | 3917 | 3.18 | 0.332 | 209 |

## Paper Use

Use this as the ecological/background layer underneath the audit: drugs connected by strong co-resistance edges share resistant subpopulation structure, so high raw AUC may reflect that shared background rather than focal-drug biology.
