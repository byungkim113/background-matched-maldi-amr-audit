# Background-Matched Transfer Audit: CNN vs Multi-task LGBM

| site | drug | cnn_raw_auc | cnn_centered_auc | cnn_raw_minus_centered | lgbm_raw_auc | lgbm_centered_auc | lgbm_raw_minus_centered | model_family_consensus |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A-2018 | Amoxicillin-Clavulanic acid | 0.650 | 0.541 | 0.109 | 0.635 | 0.482 | 0.153 | collapses after matching |
| A-2018 | Cefepime | 0.867 | 0.659 | 0.208 | 0.886 | 0.539 | 0.347 | model-dependent / mixed |
| A-2018 | Ceftazidime | 0.838 | 0.367 | 0.471 | 0.806 | 0.342 | 0.465 | low support |
| A-2018 | Ceftriaxone | 0.863 | 0.645 | 0.218 | 0.854 | 0.726 | 0.129 | survives background matching |
| A-2018 | Ciprofloxacin | 0.823 | 0.703 | 0.120 | 0.828 | 0.640 | 0.188 | survives background matching |
| A-2018 | Norfloxacin | 0.743 |  |  | 0.744 |  |  | not matched |
| DRIAMS-B | Amoxicillin-Clavulanic acid | 0.694 | 0.576 | 0.118 | 0.749 | 0.642 | 0.107 | residual signal in both models |
| DRIAMS-B | Cefepime | 0.753 |  |  | 0.844 |  |  | not matched |
| DRIAMS-B | Ceftazidime | 0.735 |  |  | 0.802 |  |  | not matched |
| DRIAMS-B | Ceftriaxone | 0.743 |  |  | 0.816 |  |  | not matched |
| DRIAMS-B | Ciprofloxacin | 0.774 | 0.860 | -0.086 | 0.816 | 0.510 | 0.306 | model-dependent / mixed |
| DRIAMS-B | Norfloxacin | 0.703 | 0.343 | 0.360 | 0.763 | 0.628 | 0.136 | model-dependent / mixed |
| DRIAMS-C | Amoxicillin-Clavulanic acid | 0.535 | 0.497 | 0.038 | 0.593 | 0.518 | 0.074 | collapses after matching |
| DRIAMS-C | Cefepime | 0.648 | 0.394 | 0.254 | 0.684 | 0.336 | 0.348 | low support |
| DRIAMS-C | Ceftazidime | 0.609 | 0.503 | 0.107 | 0.637 | 0.246 | 0.391 | low support |
| DRIAMS-C | Ceftriaxone | 0.623 | 0.727 | -0.104 | 0.659 | 0.718 | -0.058 | suggestive, low support |
| DRIAMS-C | Ciprofloxacin | 0.750 | 0.646 | 0.104 | 0.814 | 0.637 | 0.177 | survives background matching |
| DRIAMS-C | Norfloxacin | 0.779 |  |  | 0.837 |  |  | not matched |
| DRIAMS-D | Amoxicillin-Clavulanic acid | 0.557 | 0.486 | 0.070 | 0.580 | 0.449 | 0.131 | collapses after matching |
| DRIAMS-D | Cefepime | 0.769 | 0.605 | 0.164 | 0.812 | 0.557 | 0.254 | model-dependent / mixed |
| DRIAMS-D | Ceftazidime | 0.727 | 0.581 | 0.146 | 0.757 | 0.438 | 0.320 | model-dependent / mixed |
| DRIAMS-D | Ceftriaxone | 0.726 | 0.592 | 0.133 | 0.750 | 0.520 | 0.230 | model-dependent / mixed |
| DRIAMS-D | Ciprofloxacin | 0.671 | 0.596 | 0.075 | 0.723 | 0.576 | 0.147 | partial/weak residual signal |
