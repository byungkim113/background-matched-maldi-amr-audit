# Calibration Summary

Brier score, expected calibration error, and threshold metrics for each pair/site.

| site     | drug      | n    | n_r | prevalence | auc   | brier | expected_calibration_error | threshold | sensitivity | specificity | ppv   | npv   | balanced_accuracy | calibration_label      |
| -------- | --------- | ---- | --- | ---------- | ----- | ----- | -------------------------- | --------- | ----------- | ----------- | ----- | ----- | ----------------- | ---------------------- |
| A-2018   | Amox-Clav | 1367 | 383 | 0.280      | 0.650 | 0.230 | 0.194                      | 0.500     | 0.525       | 0.685       | 0.393 | 0.787 | 0.605             | borderline_calibration |
| A-2018   | FEP       | 1230 | 65  | 0.053      | 0.867 | 0.070 | 0.121                      | 0.500     | 0.692       | 0.920       | 0.326 | 0.982 | 0.806             | well_calibrated        |
| A-2018   | CAZ       | 1257 | 120 | 0.095      | 0.838 | 0.101 | 0.127                      | 0.500     | 0.658       | 0.894       | 0.397 | 0.961 | 0.776             | well_calibrated        |
| A-2018   | CRO       | 1390 | 258 | 0.186      | 0.863 | 0.130 | 0.144                      | 0.500     | 0.748       | 0.838       | 0.513 | 0.936 | 0.793             | well_calibrated        |
| A-2018   | Cipro     | 1350 | 355 | 0.263      | 0.823 | 0.171 | 0.155                      | 0.500     | 0.704       | 0.770       | 0.522 | 0.879 | 0.737             | borderline_calibration |
| A-2018   | Norflox   | 439  | 87  | 0.198      | 0.743 | 0.194 | 0.177                      | 0.500     | 0.563       | 0.747       | 0.355 | 0.874 | 0.655             | borderline_calibration |
| DRIAMS-B | Amox-Clav | 198  | 61  | 0.308      | 0.694 | 0.224 | 0.178                      | 0.500     | 0.623       | 0.672       | 0.458 | 0.800 | 0.647             | borderline_calibration |
| DRIAMS-B | FEP       | 213  | 43  | 0.202      | 0.753 | 0.143 | 0.073                      | 0.500     | 0.395       | 0.900       | 0.500 | 0.855 | 0.648             | well_calibrated        |
| DRIAMS-B | CAZ       | 213  | 45  | 0.211      | 0.735 | 0.149 | 0.076                      | 0.500     | 0.444       | 0.887       | 0.513 | 0.856 | 0.666             | well_calibrated        |
| DRIAMS-B | CRO       | 213  | 45  | 0.211      | 0.743 | 0.166 | 0.117                      | 0.500     | 0.578       | 0.827       | 0.473 | 0.880 | 0.703             | well_calibrated        |
| DRIAMS-B | Cipro     | 212  | 58  | 0.274      | 0.774 | 0.192 | 0.162                      | 0.500     | 0.672       | 0.727       | 0.481 | 0.855 | 0.700             | borderline_calibration |
| DRIAMS-B | Norflox   | 210  | 69  | 0.329      | 0.703 | 0.228 | 0.178                      | 0.500     | 0.609       | 0.716       | 0.512 | 0.789 | 0.663             | borderline_calibration |
| DRIAMS-C | Amox-Clav | 902  | 230 | 0.255      | 0.535 | 0.265 | 0.246                      | 0.500     | 0.487       | 0.549       | 0.270 | 0.758 | 0.518             | borderline_calibration |
| DRIAMS-C | FEP       | 609  | 113 | 0.186      | 0.648 | 0.175 | 0.132                      | 0.500     | 0.345       | 0.857       | 0.355 | 0.852 | 0.601             | well_calibrated        |
| DRIAMS-C | CAZ       | 908  | 139 | 0.153      | 0.609 | 0.180 | 0.161                      | 0.500     | 0.338       | 0.826       | 0.260 | 0.873 | 0.582             | borderline_calibration |
| DRIAMS-C | CRO       | 913  | 148 | 0.162      | 0.623 | 0.215 | 0.229                      | 0.500     | 0.459       | 0.735       | 0.251 | 0.875 | 0.597             | borderline_calibration |
| DRIAMS-C | Cipro     | 886  | 188 | 0.212      | 0.750 | 0.216 | 0.233                      | 0.500     | 0.702       | 0.659       | 0.357 | 0.891 | 0.681             | borderline_calibration |
| DRIAMS-C | Norflox   | 442  | 118 | 0.267      | 0.779 | 0.215 | 0.213                      | 0.500     | 0.737       | 0.667       | 0.446 | 0.874 | 0.702             | borderline_calibration |
| DRIAMS-D | Amox-Clav | 1994 | 332 | 0.166      | 0.557 | 0.259 | 0.333                      | 0.500     | 0.512       | 0.552       | 0.186 | 0.850 | 0.532             | poorly_calibrated      |
| DRIAMS-D | FEP       | 1923 | 80  | 0.042      | 0.769 | 0.107 | 0.183                      | 0.500     | 0.550       | 0.859       | 0.145 | 0.978 | 0.704             | borderline_calibration |
| DRIAMS-D | CAZ       | 1987 | 178 | 0.090      | 0.727 | 0.135 | 0.171                      | 0.500     | 0.506       | 0.833       | 0.229 | 0.945 | 0.669             | borderline_calibration |
| DRIAMS-D | CRO       | 1994 | 198 | 0.099      | 0.726 | 0.174 | 0.234                      | 0.500     | 0.561       | 0.762       | 0.206 | 0.940 | 0.661             | borderline_calibration |
| DRIAMS-D | Cipro     | 1939 | 371 | 0.191      | 0.671 | 0.240 | 0.260                      | 0.500     | 0.606       | 0.636       | 0.283 | 0.872 | 0.621             | poorly_calibrated      |
