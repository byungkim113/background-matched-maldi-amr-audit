# MARISMa v2 Metadata Feasibility Audit

Input: `/Users/byungkim/Downloads/AMR.csv`

## Main Takeaways

- Eligible organism-drug targets by metadata: **10 / 10**.
- Target-level manifest rows: **14,239**.
- Unique spectrum paths to extract for all eligible targets: **4,564**.
- Unique isolate identifiers represented in selected targets: **3,841**.
- `E. coli / Ceftriaxone` is mapped to MARISMa `E. coli / Cefotaxime` as an analogous third-generation cephalosporin phenotype.
- This is a metadata-only audit; model AUC requires extracting/preprocessing spectra for the manifest rows.

## Feasibility Table

| organism | paper_drug | marisma_drug | relationship | n_s | n_i | n_r | n_sr | r_prevalence | n_with_min_background_labels | eligible_metadata_only | prediction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Escherichia coli | Ciprofloxacin | Ciprofloxacin | exact | 952 | 77 | 807 | 1759 | 0.459 | 1757 | True | transfer_retained_or_partial |
| Escherichia coli | Norfloxacin | Norfloxacin | exact | 168 | 0 | 180 | 348 | 0.517 | 348 | True | transfer_retained_or_partial |
| Escherichia coli | Levofloxacin | Levofloxacin | external_extension | 1029 | 1 | 810 | 1839 | 0.440 | 1837 | True | transfer_retained_or_partial |
| Escherichia coli | Amoxicillin-Clavulanic acid | Amoxicillin/Clavulanic acid | spelling_alias | 1117 | 0 | 731 | 1848 | 0.396 | 1844 | True | background_sensitive_or_weak |
| Escherichia coli | Ceftriaxone | Cefotaxime | third_generation_cephalosporin_analog | 1262 | 1 | 619 | 1881 | 0.329 | 1881 | True | background_sensitive_or_mixed |
| Escherichia coli | Ceftazidime | Ceftazidime | exact | 1262 | 1 | 619 | 1881 | 0.329 | 1881 | True | background_sensitive_or_mixed |
| Escherichia coli | Cefepime | Cefepime | exact | 1296 | 4 | 584 | 1880 | 0.311 | 1878 | True | background_sensitive_or_mixed |
| Escherichia coli | Cotrimoxazole | Trimethoprim/Sulfamethoxazole | name_alias | 1134 | 10 | 711 | 1845 | 0.385 | 1843 | True | exploratory |
| Staphylococcus aureus | Oxacillin | Oxacillin | exact | 1928 | 0 | 726 | 2654 | 0.274 | 1995 | True | transfer_retained_or_partial |
| Staphylococcus epidermidis | Erythromycin | Erythromycin | exact | 184 | 0 | 458 | 642 | 0.713 | 615 | True | background_sensitive_or_weak |

## Top Co-Resistance Edges

| organism | drug_a | drug_b | n_both_known | n_rr | rr_lift | phi |
| --- | --- | --- | --- | --- | --- | --- |
| Escherichia coli | Cefotaxime | Ceftazidime | 1881 | 619 | 3.039 | 1.000 |
| Escherichia coli | Ciprofloxacin | Levofloxacin | 1758 | 806 | 2.178 | 0.999 |
| Klebsiella pneumoniae | Cefotaxime | Ceftazidime | 16659 | 4649 | 3.566 | 0.997 |
| Klebsiella pneumoniae | Cefotaxime | Cefepime | 16769 | 4523 | 3.604 | 0.981 |
| Klebsiella pneumoniae | Ceftazidime | Cefepime | 16652 | 4523 | 3.572 | 0.979 |
| Klebsiella pneumoniae | Ciprofloxacin | Levofloxacin | 14042 | 4308 | 3.153 | 0.976 |
| Escherichia coli | Cefotaxime | Cefepime | 1877 | 584 | 3.042 | 0.960 |
| Escherichia coli | Ceftazidime | Cefepime | 1877 | 584 | 3.042 | 0.960 |
| Escherichia coli | Ciprofloxacin | Norfloxacin | 334 | 158 | 1.976 | 0.936 |
| Escherichia coli | Norfloxacin | Levofloxacin | 342 | 158 | 1.943 | 0.900 |
| Klebsiella pneumoniae | Ciprofloxacin | Norfloxacin | 11619 | 2828 | 3.163 | 0.834 |
| Klebsiella pneumoniae | Ciprofloxacin | Cefotaxime | 16293 | 3833 | 3.101 | 0.799 |
| Klebsiella pneumoniae | Ciprofloxacin | Ceftazidime | 16161 | 3826 | 3.069 | 0.796 |
| Klebsiella pneumoniae | Ciprofloxacin | Cefepime | 16277 | 3779 | 3.121 | 0.796 |
| Klebsiella pneumoniae | Norfloxacin | Levofloxacin | 10281 | 2654 | 2.811 | 0.794 |
| Klebsiella pneumoniae | Levofloxacin | Cefotaxime | 14464 | 3571 | 2.788 | 0.758 |
| Klebsiella pneumoniae | Levofloxacin | Ceftazidime | 14351 | 3564 | 2.760 | 0.754 |
| Klebsiella pneumoniae | Levofloxacin | Cefepime | 14442 | 3511 | 2.807 | 0.754 |
| Staphylococcus aureus | Erythromycin | Clindamycin | 1817 | 539 | 2.295 | 0.735 |
| Klebsiella pneumoniae | Cefepime | Trimethoprim/Sulfamethoxazole | 16779 | 3633 | 2.888 | 0.712 |

## Next Step

Extract MARISMa spectra listed in `marisma_spectrum_manifest_for_selected_targets.csv`, preprocess them to the DRIAMS 6000-bin representation, export model predictions, and run the background-matched audit.
