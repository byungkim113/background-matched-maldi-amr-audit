# Weis LR Official-Panel Parity Results

This directory records the official Weis/Borgwardt LR parity check for the
site-transfer rows available in the upstream stored result JSONs.

## Scope

- Model: LR
- Seed: 35
- Train site: DRIAMS-A
- Test sites: DRIAMS-B, DRIAMS-C, DRIAMS-D where an upstream reference row exists
- Targets:
  - Escherichia coli / Ceftriaxone
  - Klebsiella pneumoniae / Ceftriaxone
  - Staphylococcus aureus / Oxacillin

## Parity Result

- Rows compared: 8
- Reference rows found: 8
- Rows within tolerance: 8
- Maximum absolute metric difference: 5.551115123125783e-17

This supports exact metric parity for the LR official subset represented here.
The associated per-target audit summaries are included as compatibility outputs,
but these single-drug parity runs do not provide a substantive co-resistance
background stress test because the exported rows have no non-focal background
drugs.

## Files

- `weis_metric_parity.csv`: row-level raw metric comparison against upstream
  stored Weis/Borgwardt JSON results.
- `weis_metric_parity.md`: compact parity summary.
- `combined_weis_raw_results.csv` and `combined_weis_raw_results.json`: raw
  metrics from the local rerun before parity comparison.
- `audit_summaries/`: per-target background-audit summaries for the official
  parity rows.
