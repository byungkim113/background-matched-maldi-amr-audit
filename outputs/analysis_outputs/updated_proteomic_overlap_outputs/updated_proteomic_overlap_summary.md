# Updated Proteomic Overlap Analysis

Bruker MALDI rows: 409
TGNRs with peak features: 407
Published ST131 peak overlaps: 36 (9 high-confidence, 15 putative, 12 loose)

## Matches by Confidence Tier

Mass-matching tolerance is 40 Da (=bin width). Confidence tiers reflect MALDI-TOF
linear-mode accuracy at 4–12 kDa. Individual protein identities should only be
claimed for high-confidence matches; the permutation enrichment result holds for all tiers.

### High confidence (≤10 Da): within typical MALDI-TOF linear-mode accuracy

- ST131 mz_9700_9725 center=9712.5 matched published m/z 9710 (HdeA, delta=2.5 Da)
- ST131 mz_4850_4875 center=4862.5 matched published m/z 4857 (HdeA multivalent ion, delta=5.5 Da)
- ST131 mz_11775_11800 center=11787.5 matched published m/z 11783 (cytochrome b562, delta=4.5 Da)
- ST131 mz_7650_7675 center=7662.5 matched published m/z 7655 (YahO, delta=7.5 Da)
- Ciprofloxacin_R mz_11775_11800 center=11787.5 matched published m/z 11783 (cytochrome b562, delta=4.5 Da)
- Ciprofloxacin_R mz_4850_4875 center=4862.5 matched published m/z 4857 (HdeA multivalent ion, delta=5.5 Da)
- Ciprofloxacin_R mz_7650_7675 center=7662.5 matched published m/z 7655 (YahO, delta=7.5 Da)
- Ceftriaxone_R mz_11775_11800 center=11787.5 matched published m/z 11783 (cytochrome b562, delta=4.5 Da)
- Ceftriaxone_R mz_4850_4875 center=4862.5 matched published m/z 4857 (HdeA multivalent ion, delta=5.5 Da)

### Putative (10–20 Da): borderline, individual protein identity uncertain

- ST131 mz_8425_8450 center=8437.5 matched published m/z 8448 (YnfD, delta=10.5 Da)
- ST131 mz_4150_4175 center=4162.5 matched published m/z 4176 (YjbJ multivalent ion, delta=13.5 Da)
- ST131 mz_6825_6850 center=6837.5 matched published m/z 6827 (unidentified, delta=10.5 Da)
- ST131 mz_6800_6825 center=6812.5 matched published m/z 6827 (unidentified, delta=14.5 Da)
- ST131 mz_4825_4850 center=4837.5 matched published m/z 4857 (HdeA multivalent ion, delta=19.5 Da)
- ST131 mz_8350_8375 center=8362.5 matched published m/z 8351 (YjbJ, delta=11.5 Da)
- ST131 mz_8325_8350 center=8337.5 matched published m/z 8351 (YjbJ, delta=13.5 Da)
- Ciprofloxacin_R mz_8325_8350 center=8337.5 matched published m/z 8351 (YjbJ, delta=13.5 Da)
- Ciprofloxacin_R mz_8425_8450 center=8437.5 matched published m/z 8448 (YnfD, delta=10.5 Da)
- Ciprofloxacin_R mz_6825_6850 center=6837.5 matched published m/z 6827 (unidentified, delta=10.5 Da)
- Ciprofloxacin_R mz_5350_5375 center=5362.5 matched published m/z 5381 (unidentified, delta=18.5 Da)
- Ceftriaxone_R mz_8425_8450 center=8437.5 matched published m/z 8448 (YnfD, delta=10.5 Da)
- Ceftriaxone_R mz_8325_8350 center=8337.5 matched published m/z 8351 (YjbJ, delta=13.5 Da)
- Ceftriaxone_R mz_6825_6850 center=6837.5 matched published m/z 6827 (unidentified, delta=10.5 Da)
- Ceftriaxone_R mz_4150_4175 center=4162.5 matched published m/z 4176 (YjbJ multivalent ion, delta=13.5 Da)

### Loose (>20 Da): exceeds typical accuracy — enrichment result valid but protein identity not claimed

- ST131 mz_6850_6875 center=6862.5 matched published m/z 6827 (unidentified, delta=35.5 Da)
- ST131 mz_9725_9750 center=9737.5 matched published m/z 9710 (HdeA, delta=27.5 Da)
- ST131 mz_8300_8325 center=8312.5 matched published m/z 8351 (YjbJ, delta=38.5 Da)
- Ciprofloxacin_R mz_6850_6875 center=6862.5 matched published m/z 6827 (unidentified, delta=35.5 Da)
- Ciprofloxacin_R mz_9725_9750 center=9737.5 matched published m/z 9710 (HdeA, delta=27.5 Da)
- Ciprofloxacin_R mz_8300_8325 center=8312.5 matched published m/z 8351 (YjbJ, delta=38.5 Da)
- Ceftriaxone_R mz_6850_6875 center=6862.5 matched published m/z 6827 (unidentified, delta=35.5 Da)
- Ceftriaxone_R mz_9725_9750 center=9737.5 matched published m/z 9710 (HdeA, delta=27.5 Da)
- Ceftriaxone_R mz_8300_8325 center=8312.5 matched published m/z 8351 (YjbJ, delta=38.5 Da)
- Ceftriaxone_R mz_4200_4225 center=4212.5 matched published m/z 4176 (YjbJ multivalent ion, delta=36.5 Da)
- Ceftriaxone_R mz_8375_8400 center=8387.5 matched published m/z 8351 (YjbJ, delta=36.5 Da)
- Ceftriaxone_R mz_7675_7700 center=7687.5 matched published m/z 7655 (YahO, delta=32.5 Da)

## Mass-Matched Permutation Enrichment

The null model samples the same number of peak bins with the same coarse m/z-stratum
distribution, then counts overlaps with the published ST131 MALDI biomarker list.
Enrichment is computed over all 36 matches (all tiers); the result is robust because
the null distribution is matched by the same 40 Da tolerance.

- ST131: observed=14/75, null_mean=4.5065, fold=3.106624, p_ge=0.0001
- Ciprofloxacin_R: observed=10/75, null_mean=4.4738, fold=2.235236, p_ge=0.006799
- Ceftriaxone_R: observed=12/75, null_mean=4.539, fold=2.643754, p_ge=0.0006
- ALL_TARGETS: observed=36/225, null_mean=13.5504, fold=2.656748, p_ge=0.0001
