# Background Audit With Resistance Ecology

This table links the background-matched audit to the co-resistance network inferred from the same isolate label matrix.

| site | drug | resistance_ecology_block | strongest_network_partner | partner_phi | cnn_raw_auc | cnn_centered_auc | lgbm_raw_auc | lgbm_centered_auc | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A-2018 | Amoxicillin-Clavulanic acid | mixed beta-lactam/background | Ceftazidime | 0.379 | 0.650 | 0.541 | 0.635 | 0.482 | background-driven in both model families |
| A-2018 | Cefepime | cephalosporin/ESBL block | Ceftazidime | 0.828 | 0.867 | 0.659 | 0.886 | 0.539 | partial residual focal signal |
| A-2018 | Ceftazidime | cephalosporin/ESBL block | Ceftriaxone | 0.884 | 0.838 | 0.367 | 0.806 | 0.342 | background-driven in both model families |
| A-2018 | Ceftriaxone | cephalosporin/ESBL block | Ceftazidime | 0.884 | 0.863 | 0.645 | 0.854 | 0.726 | focal signal survives background matching |
| A-2018 | Ciprofloxacin | fluoroquinolone block | Norfloxacin | 0.976 | 0.823 | 0.703 | 0.828 | 0.640 | focal signal survives background matching |
| A-2018 | Norfloxacin | fluoroquinolone block | Ciprofloxacin | 0.976 | 0.743 |  | 0.744 |  | insufficient matched support |
| DRIAMS-B | Amoxicillin-Clavulanic acid | mixed beta-lactam/background | Ceftazidime | 0.379 | 0.694 | 0.576 | 0.749 | 0.642 | partial residual focal signal |
| DRIAMS-B | Cefepime | cephalosporin/ESBL block | Ceftazidime | 0.828 | 0.753 |  | 0.844 |  | insufficient matched support |
| DRIAMS-B | Ceftazidime | cephalosporin/ESBL block | Ceftriaxone | 0.884 | 0.735 |  | 0.802 |  | insufficient matched support |
| DRIAMS-B | Ceftriaxone | cephalosporin/ESBL block | Ceftazidime | 0.884 | 0.743 |  | 0.816 |  | insufficient matched support |
| DRIAMS-B | Ciprofloxacin | fluoroquinolone block | Norfloxacin | 0.976 | 0.774 | 0.860 | 0.816 | 0.510 | partial residual focal signal |
| DRIAMS-B | Norfloxacin | fluoroquinolone block | Ciprofloxacin | 0.976 | 0.703 | 0.343 | 0.763 | 0.628 | partial residual focal signal |
| DRIAMS-C | Amoxicillin-Clavulanic acid | mixed beta-lactam/background | Ceftazidime | 0.379 | 0.535 | 0.497 | 0.593 | 0.518 | background-driven in both model families |
| DRIAMS-C | Cefepime | cephalosporin/ESBL block | Ceftazidime | 0.828 | 0.648 | 0.394 | 0.684 | 0.336 | background-driven in both model families |
| DRIAMS-C | Ceftazidime | cephalosporin/ESBL block | Ceftriaxone | 0.884 | 0.609 | 0.503 | 0.637 | 0.246 | background-driven in both model families |
| DRIAMS-C | Ceftriaxone | cephalosporin/ESBL block | Ceftazidime | 0.884 | 0.623 | 0.727 | 0.659 | 0.718 | focal signal survives background matching |
| DRIAMS-C | Ciprofloxacin | fluoroquinolone block | Norfloxacin | 0.976 | 0.750 | 0.646 | 0.814 | 0.637 | focal signal survives background matching |
| DRIAMS-C | Norfloxacin | fluoroquinolone block | Ciprofloxacin | 0.976 | 0.779 |  | 0.837 |  | insufficient matched support |
| DRIAMS-D | Amoxicillin-Clavulanic acid | mixed beta-lactam/background | Ceftazidime | 0.379 | 0.557 | 0.486 | 0.580 | 0.449 | background-driven in both model families |
| DRIAMS-D | Cefepime | cephalosporin/ESBL block | Ceftazidime | 0.828 | 0.769 | 0.605 | 0.812 | 0.557 | partial residual focal signal |
| DRIAMS-D | Ceftazidime | cephalosporin/ESBL block | Ceftriaxone | 0.884 | 0.727 | 0.581 | 0.757 | 0.438 | partial residual focal signal |
| DRIAMS-D | Ceftriaxone | cephalosporin/ESBL block | Ceftazidime | 0.884 | 0.726 | 0.592 | 0.750 | 0.520 | partial residual focal signal |
| DRIAMS-D | Ciprofloxacin | fluoroquinolone block | Norfloxacin | 0.976 | 0.671 | 0.596 | 0.723 | 0.576 | partial residual focal signal |
