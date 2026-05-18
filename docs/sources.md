# Sources And Upstream Data

This file records the sources that the repository depends on or cites. It is meant to make the paper repository auditable: code and derived outputs live here, while raw data remain with their original providers unless explicitly noted.

## DRIAMS

- Dataset/project: DRIAMS MALDI-TOF antimicrobial resistance benchmark.
- Primary paper: Weis et al., "Direct antimicrobial resistance prediction from clinical MALDI-TOF mass spectra using machine learning," *Nature Medicine*.
- Upstream use in this repository: source spectra/AST labels for the DRIAMS-A training site, temporal A-2018 testing, and DRIAMS-B/C/D external-site evaluation.
- Redistribution note: raw DRIAMS spectra and raw AST exports are not included in this repository. Users must obtain them independently from the official DRIAMS data release.

Useful links:

- Nature Medicine article: <https://www.nature.com/articles/s41591-021-01619-9>
- DRIAMS Dryad dataset DOI: <https://doi.org/10.5061/dryad.bzkh1899q>

## Weis/Borgwardt MALDI-AMR Code And Stratification Work

- Repository: BorgwardtLab `maldi_amr`.
- Use in this repository: comparison target for Weis-style prediction export and model-agnostic audit compatibility.
- Related preprint: "Improved MALDI-TOF MS based antimicrobial resistance prediction through hierarchical stratification," bioRxiv DOI `10.1101/2022.04.13.488198`.
- Redistribution note: this repository does not vendor the Weis/Borgwardt repository; the Kaggle notebook clones it at runtime when needed.

Useful links:

- GitHub: <https://github.com/BorgwardtLab/maldi_amr>
- bioRxiv DOI: <https://doi.org/10.1101/2022.04.13.488198>

## Public Basel/Cuenod UPEC WGS-Linked MALDI Data

- ENA project: `PRJEB55855`
- Secondary study accession: `ERP140793`
- Study title: "Bacterial genome wide association study substantiates papGII of E. coli as a patient independent driver of urosepsis."
- GitHub repository: `acuenod111/UPEC`.
- OSF project used for MALDI files and Bruker-derived tables: `vmqc5`.
- Use in this repository: public WGS-linked support analysis showing that Bruker MALDI peak features encode ST131 lineage and resistance-associated signals.
- Redistribution note: this repository includes selected processed manifest/feature tables used by the current analysis snapshot. It does not include raw FASTQ files, raw Bruker FID archives, or raw OSF tarballs.

Useful links:

- ENA project: <https://www.ebi.ac.uk/ena/browser/view/PRJEB55855>
- ENA study: <https://www.ebi.ac.uk/ena/browser/view/ERP140793>
- GitHub: <https://github.com/acuenod111/UPEC>
- OSF: <https://osf.io/vmqc5/>

## ST131 MALDI Biomarker Reference

- Reference: Scientific Reports 2019 ST131 MALDI biomarker paper, DOI `10.1038/s41598-019-45051-z`.
- Use in this repository: published ST131 marker masses are used as a literature-defined biomarker set for permutation enrichment against discriminative MALDI peak bins.
- Important limitation: overlap with published biomarker masses is supporting evidence, not direct protein identification of DRIAMS peaks. Definitive peak identity would require MS/MS, MALDI-TOF/TOF, LC-MS/MS, or comparable proteomic validation.

Useful links:

- Article DOI: <https://doi.org/10.1038/s41598-019-45051-z>

## MARISMa External Bruker Snapshot

- Dataset: MARISMa Bruker MALDI-TOF spectra and AMR labels.
- Kaggle snapshot used in the current workflow: `bfdf121/marisma`.
- Use in this repository: independent external stress test for a locked
  DRIAMS-trained Mega/CNN checkpoint. The derived prediction table and audit
  outputs are stored under
  `outputs/analysis_outputs/marisma_external_validation/`.
- Important limitation: the AMR file used for this snapshot contains 2018-2024
  records overall, but the three target organisms evaluated here appear only in
  2024. Treat the current result as an external dataset stress test, not a
  longitudinal MARISMa deployment study.
- Redistribution note: raw MARISMa spectra and raw AMR labels are not included
  in this repository. Users should obtain them from the upstream provider or
  Kaggle snapshot.

Useful links:

- Kaggle snapshot: <https://www.kaggle.com/datasets/bfdf121/marisma>
- Zenodo record discussed for MARISMa v2: <https://zenodo.org/records/17201597>

## Notes Before Public Release

When preparing this repository for public release alongside the manuscript,
confirm upstream redistribution terms for each data source listed above and
add the manuscript DOI to CITATION.cff once the paper is accepted.
