# Literature Positioning For The Manuscript

## What Prior MALDI-AMR Work Established

**DRIAMS / Weis et al. Nature Medicine.**
The DRIAMS work established that MALDI-TOF spectra plus machine learning can
predict antimicrobial resistance for many organism-drug pairs under careful
temporal and external-site evaluation. It also made clear that performance is
not uniform across sites and pairs. Our manuscript builds on that benchmark but
changes the question from "what is the AUC?" to "what part of the AUC survives
after controlling for resistant-population background?"

**Weis et al. hierarchical stratification preprint.**
This work explicitly raised the possibility that MALDI spectra can carry
phylogenetic/background structure and used hierarchical clustering to modify
evaluation splits. Our work should cite it directly. The difference is that we
turn the confounder into an operational audit: given predictions from any model,
we quantify raw AUC, matched AUC, background-centered AUC, matched retention and
cross-resistance structure.

**WGS-AMR population-structure literature.**
Yu, Wheeler and Barquist show that bacterial population structure can confound
machine-learning AMR prediction from genomic data. Our manuscript is best framed
as extending the same principle to routine MALDI-TOF AMR evaluation, where
lineage and co-resistance can be spectral shortcuts.

**ST131 MALDI biomarker literature.**
Published ST131 MALDI biomarker studies show that MALDI spectra can encode
lineage-associated peaks. Our public UPEC analysis uses this as biological
support: lineage is strongly predictable from MALDI peak features, and
resistance-associated discriminative peaks are enriched for published ST131
biomarker masses. This does not prove DRIAMS detects ST131.

## Novelty Claim

The novelty is not "population structure exists." The novelty is the
model-agnostic audit that measures whether focal-drug prediction remains after
co-resistance background is controlled.

The strongest claim:

> Background-matched evaluation exposes when MALDI-TOF AMR performance is
> retained within co-resistance strata versus when it is largely explained by
> resistant-population background.

## Nature Communications Framing

For Nature Communications, the paper needs to feel broader than a single
DRIAMS model. The current draft frames the contribution as:

1. A general audit method.
2. A DRIAMS CNN application.
3. A model-family replication with LightGBM.
4. A cross-resistance network showing the exploitable label structure.
5. Public WGS-linked MALDI evidence that lineage is strongly encoded.
6. A reusable software package for other MALDI-AMR models.

## Limitations And Open Questions

- No WGS/MLST labels linked directly to the DRIAMS isolates.
- No direct MS/MS or LC-MS/MS identification of discriminative peaks.
- Some drug-site rows have low matched retention and must be flagged.
- Public UPEC WGS support is biological plausibility evidence, not a full
  second deployment validation of the DRIAMS-trained model.
- Weis-style model auditing is useful supplementary evidence but should not be
  overclaimed unless rerun as a locked full-row export.
