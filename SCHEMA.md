# Background-Matched Audit — Input/Output Schema

Any model's predictions can be audited by providing a single long CSV file.
No DRIAMS data, no PyTorch, and no MALDI-specific code are required.

---

## Input: Predictions CSV

One row per isolate-drug prediction. All drug predictions for the same isolate
must share the same `isolate_id`, `site`, `year`, and `organism`.

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `isolate_id` | string | Yes | Unique identifier for the isolate/sample |
| `site` | string | Yes | Hospital or lab site (e.g. `DRIAMS-A`, `Hospital-B`) |
| `year` | string | Yes | Collection year (e.g. `2018`) |
| `organism` | string | Yes | Species name (e.g. `Escherichia coli`) |
| `drug` | string | Yes | Antibiotic name (e.g. `Ciprofloxacin`) |
| `label` | int or string | Yes | Resistance outcome: `R` or `1` = resistant; `S` or `0` = susceptible |
| `prob` | float [0, 1] | Yes | Model probability of resistance |
| `model_name` | string | No | Model identifier (used when combining multiple models) |
| `background_signature` | string | No | Precomputed co-resistance signature (see below); if absent, derived automatically |

### Label encoding

The `label` column accepts: `R`, `S`, `1`, `0`, `resistant`, `susceptible`,
`True`, `False`, `Yes`, `No`. Intermediate (`I`), unknown (`U`, `-`, `NA`),
or non-parseable values are dropped.

### Precomputed background signatures (optional)

If supplied, `background_signature` must be a pipe-delimited string of
`Drug=STATUS` pairs for all **other** drugs tested on the same isolate,
excluding the focal drug. Status values: `R`, `S`, or `U` (unknown).

Example: `Norfloxacin=R|Ceftriaxone=S|Amoxicillin-Clavulanic acid=U`

If this column is absent, the framework derives signatures automatically from
other rows that share the same `isolate_id`, `site`, `year`, and `organism`.

---

## Output files

### `background_matched_audit_summary.csv`

One row per (site × organism × drug) combination.

| Column | Type | Description |
|--------|------|-------------|
| `site` | string | Evaluation site |
| `organism` | string | Species |
| `drug` | string | Focal antibiotic |
| `raw_auc` | float | AUROC without background matching |
| `raw_aupr` | float | AUPR without background matching |
| `matched_auc` | float | AUROC restricted to matched strata |
| `stratum_centered_auc` | float | AUROC after subtracting per-stratum mean prediction; the primary audit metric |
| `pairwise_accuracy` | float | Fraction of R/S pairs ranked correctly within the same background stratum |
| `pairwise_comparisons` | int | Number of within-stratum R/S pairs evaluated |
| `matched_retention` | float | Fraction of isolates in strata that contain both R and S isolates |
| `n_total` | int | All isolates for this site/drug |
| `n_r` | int | Resistant isolates |
| `n_matched` | int | Isolates in matched strata |
| `n_matched_r` | int | Resistant isolates in matched strata |
| `n_valid_strata` | int | Number of strata meeting the min-pos and min-neg thresholds |
| `min_pos_per_stratum` | int | Minimum resistant isolates required per stratum |
| `min_neg_per_stratum` | int | Minimum susceptible isolates required per stratum |
| `adequacy` | string | `ok` or `caution` (low matched support) |

### `background_matched_retained_rows.csv`

Individual isolate rows within matched strata. Includes all input columns plus:

| Column | Description |
|--------|-------------|
| `background_signature` | Co-resistance profile (all background drugs, focal excluded) |
| `background_known_count` | Number of background drugs with a known R/S label |
| `background_resistant_count` | Number of resistant background drugs |
| `centered_prob` | `prob` minus per-stratum mean probability |
| `matched_valid_stratum` | Always `True` in this file |

### `background_matched_sensitivity.csv` *(from `sensitivity_sweep.py`)*

One row per minimum-stratum threshold.

| Column | Description |
|--------|-------------|
| `min_stratum` | Minimum R and S isolates required per stratum |
| `n_drug_site_pairs` | Total site × drug pairs evaluated |
| `n_adequate` | Pairs with at least one valid stratum |
| `mean_matched_retention` | Mean fraction of isolates retained across pairs |
| `macro_raw_auc` | Mean raw AUC across all pairs |
| `macro_centered_auc` | Mean stratum-centered AUC across adequate pairs |
| `mean_delta` | Mean (raw − centered) attenuation across adequate pairs |

### `cross_resistance_edges.csv`

One row per drug-pair co-resistance edge (all sites combined and per site).

| Column | Description |
|--------|-------------|
| `site` | Site or `ALL` |
| `drug_a`, `drug_b` | Drug pair (alphabetical order) |
| `n_both_known` | Isolates with known R/S labels for both drugs |
| `n_rr` | Double-resistant isolates |
| `rr_lift` | Observed double-resistance rate / (prev_a × prev_b); > 1 means positive co-resistance |
| `phi` | Phi (Matthews) correlation between resistance labels |
| `resistant_jaccard` | Jaccard similarity of resistant-isolate sets |

---

## Minimal example

```
isolate_id,site,year,organism,drug,label,prob
ISO_001,Hospital-A,2021,Escherichia coli,Ciprofloxacin,R,0.83
ISO_001,Hospital-A,2021,Escherichia coli,Ceftriaxone,S,0.41
ISO_002,Hospital-A,2021,Escherichia coli,Ciprofloxacin,S,0.22
ISO_002,Hospital-A,2021,Escherichia coli,Ceftriaxone,S,0.29
```

Run:
```bash
python run_background_audit_framework.py \
    --predictions-csv example_predictions.csv \
    --output-dir outputs/my_audit
```

---

## Interpretation guide

| stratum_centered_auc | Interpretation |
|----------------------|---------------|
| ≥ 0.65 | Focal-drug signal likely survives background matching |
| 0.55 – 0.65 | Partial retention; signal exists but attenuates |
| ≤ 0.55 | Near-chance after matching; model likely exploiting co-resistance background |

The delta (raw − centered) quantifies attenuation. A large positive delta (> 0.10)
with low centered AUC (< 0.60) is the clearest evidence of background-driven inflation.

`matched_retention` below 10 % means the strata constraint is too strict for the available
data and results should be interpreted with caution (the `adequacy` column will say `caution`).
