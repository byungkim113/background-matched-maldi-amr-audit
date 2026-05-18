# Public UPEC WGS-Linked MALDI Support

This project includes two complementary public Basel/Cuenod UPEC analyses. They are not duplicate analyses:

1. `scripts/upec_wgs_validation_analysis.py`
   - Links public UPEC WGS metadata to Bruker MALDI median-peak features.
   - Tests whether MALDI peak features encode lineage and resistance labels.
   - Main outputs:
     - `outputs/analysis_outputs/upec_wgs_validation_outputs/validation_summary.md`
     - `outputs/analysis_outputs/upec_wgs_validation_outputs/centroid_binary_cv_results.csv`
     - `outputs/analysis_outputs/upec_wgs_validation_outputs/st131_resistance_associations.csv`

2. `scripts/updated_proteomic_overlap_analysis.py`
   - Uses the same public UPEC MALDI/WGS bridge.
   - Takes discriminative peak bins for ST131, ciprofloxacin resistance, and ceftriaxone resistance.
   - Tests whether those bins are enriched for published ST131 MALDI biomarker masses using a mass-stratified permutation null.
   - Main outputs:
     - `outputs/analysis_outputs/updated_proteomic_overlap_outputs/updated_proteomic_overlap_summary.md`
     - `outputs/analysis_outputs/updated_proteomic_overlap_outputs/updated_proteomic_overlap_permutation_enrichment.csv`
     - `outputs/analysis_outputs/updated_proteomic_overlap_outputs/updated_published_st131_proteomic_overlap.csv`

## Current Best Results

WGS-linked MALDI validation:

| Target | AUC |
| --- | ---: |
| ST131 | 0.9318 |
| Ciprofloxacin resistance | 0.7548 |
| Ceftriaxone resistance | 0.6888 |

Published ST131 biomarker enrichment:

| Target | Observed overlaps | Null mean | Fold enrichment | Empirical p |
| --- | ---: | ---: | ---: | ---: |
| ST131 | 14/75 | 4.5065 | 3.1066 | 0.0001 |
| Ciprofloxacin resistance | 10/75 | 4.4738 | 2.2352 | 0.006799 |
| Ceftriaxone resistance | 12/75 | 4.5390 | 2.6438 | 0.0006 |
| All targets | 36/225 | 13.5504 | 2.6567 | 0.0001 |

## Interpretation Guardrail

These analyses support biological plausibility: public Bruker MALDI spectra encode ST131 lineage strongly, and resistance-associated discriminative peak bins are enriched for published ST131 biomarker masses.

They do not prove that the DRIAMS model directly detects ST131, nor do they identify the exact proteins behind every DRIAMS saliency peak. Exact protein identification would require MS/MS or related proteomic validation.

## Re-run

```bash
python scripts/run_public_upec_analysis.py \
  --median-peaks data_manifests/Bruker_csv_medianpeaks_df.csv
```

The wrapper runs both the WGS-linked validation and the updated proteomic overlap analysis.
