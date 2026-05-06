# Temporal Reliability Audit

Period-wise reliability monitor for deciding when recalibration or retraining should be reviewed.

| site     | drug      | period | n_periods_observed | n    | n_r | prevalence | auc   | brier | expected_calibration_error | mean_background_resistant_count | support_label | reliability_status   | recommended_action     |
| -------- | --------- | ------ | ------------------ | ---- | --- | ---------- | ----- | ----- | -------------------------- | ------------------------------- | ------------- | -------------------- | ---------------------- |
| A-2018   | Amox-Clav | 2018   | 1                  | 1367 | 383 | 0.280      | 0.650 | 0.230 | 0.194                      | 0.633                           | adequate      | insufficient_periods | collect_future_periods |
| A-2018   | FEP       | 2018   | 1                  | 1230 | 65  | 0.053      | 0.867 | 0.070 | 0.121                      | 0.633                           | adequate      | insufficient_periods | collect_future_periods |
| A-2018   | CAZ       | 2018   | 1                  | 1257 | 120 | 0.095      | 0.838 | 0.101 | 0.127                      | 0.696                           | adequate      | insufficient_periods | collect_future_periods |
| A-2018   | CRO       | 2018   | 1                  | 1390 | 258 | 0.186      | 0.863 | 0.130 | 0.144                      | 0.722                           | adequate      | insufficient_periods | collect_future_periods |
| A-2018   | Cipro     | 2018   | 1                  | 1350 | 355 | 0.263      | 0.823 | 0.171 | 0.155                      | 0.643                           | adequate      | insufficient_periods | collect_future_periods |
| A-2018   | Norflox   | 2018   | 1                  | 439  | 87  | 0.198      | 0.743 | 0.194 | 0.177                      | 0.588                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | Amox-Clav | 2018   | 1                  | 198  | 61  | 0.308      | 0.694 | 0.224 | 0.178                      | 1.288                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | FEP       | 2018   | 1                  | 213  | 43  | 0.202      | 0.753 | 0.143 | 0.073                      | 1.305                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | CAZ       | 2018   | 1                  | 213  | 45  | 0.211      | 0.735 | 0.149 | 0.076                      | 1.296                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | CRO       | 2018   | 1                  | 213  | 45  | 0.211      | 0.743 | 0.166 | 0.117                      | 1.296                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | Cipro     | 2018   | 1                  | 212  | 58  | 0.274      | 0.774 | 0.192 | 0.162                      | 1.217                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-B | Norflox   | 2018   | 1                  | 210  | 69  | 0.329      | 0.703 | 0.228 | 0.178                      | 1.200                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | Amox-Clav | 2018   | 1                  | 902  | 230 | 0.255      | 0.535 | 0.265 | 0.246                      | 0.783                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | FEP       | 2018   | 1                  | 609  | 113 | 0.186      | 0.648 | 0.175 | 0.132                      | 1.034                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | CAZ       | 2018   | 1                  | 908  | 139 | 0.153      | 0.609 | 0.180 | 0.161                      | 0.859                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | CRO       | 2018   | 1                  | 913  | 148 | 0.162      | 0.623 | 0.215 | 0.229                      | 0.852                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | Cipro     | 2018   | 1                  | 886  | 188 | 0.212      | 0.750 | 0.216 | 0.233                      | 0.837                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-C | Norflox   | 2018   | 1                  | 442  | 118 | 0.267      | 0.779 | 0.215 | 0.213                      | 1.129                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-D | Amox-Clav | 2018   | 1                  | 1994 | 332 | 0.166      | 0.557 | 0.259 | 0.333                      | 0.415                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-D | FEP       | 2018   | 1                  | 1923 | 80  | 0.042      | 0.769 | 0.107 | 0.183                      | 0.454                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-D | CAZ       | 2018   | 1                  | 1987 | 178 | 0.090      | 0.727 | 0.135 | 0.171                      | 0.485                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-D | CRO       | 2018   | 1                  | 1994 | 198 | 0.099      | 0.726 | 0.174 | 0.234                      | 0.480                           | adequate      | insufficient_periods | collect_future_periods |
| DRIAMS-D | Cipro     | 2018   | 1                  | 1939 | 371 | 0.191      | 0.671 | 0.240 | 0.260                      | 0.381                           | adequate      | insufficient_periods | collect_future_periods |
