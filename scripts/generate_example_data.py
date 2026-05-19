#!/usr/bin/env python3
"""Generate synthetic example prediction data for the background-matched audit.

Produces example_predictions.csv at the repository root with two organisms,
demonstrating three archetypal audit outcomes:

  E. coli / Ciprofloxacin   — POSITIVE CONTROL: drug-specific signal survives
                               background matching (chromosomal gyrA/parC mutations).
  E. coli / Ceftriaxone     — NEGATIVE CONTROL: signal collapses after matching
                               (model exploits fluoroquinolone co-resistance block).
  S. aureus / Oxacillin     — SECOND ORGANISM POSITIVE CONTROL: mecA/PBP2a structural
                               phenotype gives focal-drug signal even when matched on
                               Penicillin and Ciprofloxacin co-resistance background.
  S. aureus / Erythromycin  — SECOND ORGANISM NEGATIVE CONTROL: mobile erm/msr
                               resistance co-segregates with the MRSA background, so
                               the model exploits background rather than focal signal.

The data are fully synthetic (fake isolate IDs, simulated probabilities).
No real patient data are included. Run:

    python scripts/generate_example_data.py
    python run_background_audit_framework.py \\
        --predictions-csv example_predictions.csv \\
        --output-dir outputs/example_run
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path


SEED = 20260427
N_ECOLI_PER_SITE = 150
N_SAUREUS_PER_SITE = 160
SITES = [("Hospital-A", "2021"), ("Hospital-B", "2022")]
OUT = Path(__file__).resolve().parents[1] / "example_predictions.csv"
FIELDNAMES = ["isolate_id", "site", "year", "organism", "drug", "label", "prob", "model_name"]

ECOLI_DRUGS = [
    "Ciprofloxacin",
    "Norfloxacin",
    "Ceftriaxone",
    "Amoxicillin-Clavulanic acid",
]

SAUREUS_DRUGS = [
    "Oxacillin",
    "Penicillin",
    "Ciprofloxacin",
    "Erythromycin",
]


def sigmoid(x: float) -> float:
    x = max(-20.0, min(20.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def clamp(x: float, lo: float = 0.02, hi: float = 0.98) -> float:
    return max(lo, min(hi, x))


def generate_ecoli_isolate(rng: random.Random, site: str, year: str, idx: int) -> list[dict]:
    isolate_id = f"EC_{site[:3].upper()}_{idx:04d}"

    # Fluoroquinolone block (~30% cipro-R; norfloxacin almost always co-resistant).
    cipro_r = int(rng.random() < 0.30)
    norflox_r = cipro_r if rng.random() < 0.92 else 1 - cipro_r

    # Beta-lactam / ESBL block (~20% ceftriaxone-R; amox-clav highly co-resistant).
    # Weak cross-block association (OR ~2).
    cef_prev = 0.28 if cipro_r else 0.15
    cef_r = int(rng.random() < cef_prev)
    amox_r = cef_r if rng.random() < 0.80 else int(rng.random() < 0.08)

    labels = {
        "Ciprofloxacin":             cipro_r,
        "Norfloxacin":               norflox_r,
        "Ceftriaxone":               cef_r,
        "Amoxicillin-Clavulanic acid": amox_r,
    }

    # Ciprofloxacin: drug-specific — signal survives within co-resistance strata.
    cipro_prob  = clamp(sigmoid(2.8 * cipro_r  - 1.2 + rng.gauss(0, 0.55)))
    # Norfloxacin: same chromosomal block — also survives.
    norflox_prob = clamp(sigmoid(2.5 * norflox_r - 1.1 + rng.gauss(0, 0.55)))
    # Ceftriaxone: background-driven — tracks cipro/norflox resistance, not label.
    cef_prob    = clamp(sigmoid(2.3 * cipro_r + 1.9 * norflox_r - 1.6 + rng.gauss(0, 0.6)))
    # Amox-Clav: similarly background-driven.
    amox_prob   = clamp(sigmoid(1.9 * cipro_r + 1.6 * norflox_r - 1.3 + rng.gauss(0, 0.65)))

    probs = {
        "Ciprofloxacin":             cipro_prob,
        "Norfloxacin":               norflox_prob,
        "Ceftriaxone":               cef_prob,
        "Amoxicillin-Clavulanic acid": amox_prob,
    }

    return [
        dict(isolate_id=isolate_id, site=site, year=year, organism="Escherichia coli",
             drug=drug, label=labels[drug], prob=round(probs[drug], 4),
             model_name="SyntheticModel-v1")
        for drug in ECOLI_DRUGS
    ]


def generate_saureus_isolate(rng: random.Random, site: str, year: str, idx: int) -> list[dict]:
    """
    S. aureus resistance structure:
    - mecA carrier (~22%): causes oxacillin-R (MRSA) + often penicillin-R + elevated cipro-R.
    - blaZ carrier (~72%): causes penicillin-R but NOT oxacillin-R.
    - erm/msr carrier (~32%): causes erythromycin-R; co-segregates with mecA (OR ~2.5).

    Model design:
    - Oxacillin:   strong focal-drug signal (mecA is a large, stable structural phenotype).
    - Erythromycin: background-driven (model latches onto MRSA background, not erm signal).
    - Penicillin:  near-universal (blaZ), model score mostly uninformative.
    - Ciprofloxacin: moderate focal-drug signal with some mecA cross-signal.
    """
    isolate_id = f"SA_{site[:3].upper()}_{idx:04d}"

    meca = int(rng.random() < 0.22)            # MRSA carrier
    blaz = int(rng.random() < 0.72)            # blaZ penicillinase
    erm  = int(rng.random() < (0.45 if meca else 0.24))  # erm/msr (enriched in MRSA)

    oxa_r  = meca
    pen_r  = meca or blaz                       # both mechanisms cause pen-R
    cip_r  = int(rng.random() < (0.40 if meca else 0.12))
    ery_r  = erm

    labels = {
        "Oxacillin":    oxa_r,
        "Penicillin":   int(pen_r),
        "Ciprofloxacin": cip_r,
        "Erythromycin": ery_r,
    }

    # Oxacillin: strong drug-specific signal from mecA structural phenotype.
    # Within any background stratum it still discriminates.
    oxa_prob = clamp(sigmoid(3.1 * meca - 1.4 + rng.gauss(0, 0.50)))

    # Erythromycin: model exploits the MRSA/mecA background proteomic signature.
    # No direct focal-drug coefficient — collapses after matching because all
    # discriminating information is shared with the oxa/pen co-resistance background.
    ery_prob = clamp(sigmoid(2.2 * meca + 0.4 * pen_r - 1.0 + rng.gauss(0, 0.70)))

    # Penicillin: near-universal, very high baseline probability, weak focal signal.
    pen_prob = clamp(sigmoid(0.5 * meca + 0.7 * blaz + 0.3 + rng.gauss(0, 0.70)))

    # Ciprofloxacin: moderate focal-drug signal with partial mecA cross-signal.
    cip_prob = clamp(sigmoid(1.8 * cip_r + 0.8 * meca - 1.3 + rng.gauss(0, 0.60)))

    probs = {
        "Oxacillin":    oxa_prob,
        "Penicillin":   pen_prob,
        "Ciprofloxacin": cip_prob,
        "Erythromycin": ery_prob,
    }

    return [
        dict(isolate_id=isolate_id, site=site, year=year,
             organism="Staphylococcus aureus",
             drug=drug, label=labels[drug], prob=round(probs[drug], 4),
             model_name="SyntheticModel-v1")
        for drug in SAUREUS_DRUGS
    ]


def main() -> None:
    rng = random.Random(SEED)
    rows: list[dict] = []

    for site, year in SITES:
        for i in range(N_ECOLI_PER_SITE):
            rows.extend(generate_ecoli_isolate(rng, site, year, i + 1))
        for i in range(N_SAUREUS_PER_SITE):
            rows.extend(generate_saureus_isolate(rng, site, year, i + 1))

    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    n_ecoli  = N_ECOLI_PER_SITE  * len(SITES)
    n_sa     = N_SAUREUS_PER_SITE * len(SITES)
    n_drugs  = len(ECOLI_DRUGS) + len(SAUREUS_DRUGS)
    print(f"Wrote {len(rows)} rows ({n_ecoli} E. coli + {n_sa} S. aureus isolates "
          f"x {len(ECOLI_DRUGS)}/{len(SAUREUS_DRUGS)} drugs) -> {OUT}")
    print()
    print("Expected audit outcomes:")
    print("  E. coli / Ciprofloxacin   [POSITIVE CONTROL]  survives matching (chromosomal gyrA/parC)")
    print("  E. coli / Norfloxacin     [POSITIVE CONTROL]  survives (same fluoroquinolone block)")
    print("  E. coli / Ceftriaxone     [NEGATIVE CONTROL]  collapses (background-driven via ESBL block)")
    print("  E. coli / Amox-Clav       [NEGATIVE CONTROL]  collapses (heterogeneous beta-lactamase context)")
    print("  S. aureus / Oxacillin     [POSITIVE CONTROL]  survives (mecA/PBP2a structural phenotype)")
    print("  S. aureus / Ciprofloxacin [INTERMEDIATE]      partial retention with mecA cross-signal")
    print("  S. aureus / Erythromycin  [NEGATIVE CONTROL]  collapses (model exploits MRSA background)")
    print("  S. aureus / Penicillin    [NEAR-UNINFORMATIVE] near-universal, weak focal signal")
    print()
    print("Run the audit with:")
    print(f"  python run_background_audit_framework.py \\")
    print(f"    --predictions-csv {OUT.name} \\")
    print(f"    --output-dir outputs/example_run")
    print()
    print("Run sensitivity sweep:")
    print(f"  python scripts/sensitivity_sweep.py \\")
    print(f"    --predictions-csv {OUT.name} \\")
    print(f"    --output-dir outputs/sensitivity_sweep \\")
    print(f"    --site-robustness")


if __name__ == "__main__":
    main()
