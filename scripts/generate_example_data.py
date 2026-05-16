#!/usr/bin/env python3
"""Generate synthetic example prediction data for the background-matched audit.

Produces example_predictions.csv at the repository root, demonstrating two
archetypal outcomes:

  Ciprofloxacin  — drug-specific signal survives background matching.
  Ceftriaxone    — signal collapses after matching (background-driven).

The data are fully synthetic (fake isolate IDs, simulated probabilities).
No real patient data are included. Run:

    python scripts/generate_example_data.py
    python run_background_audit_framework.py \
        --predictions-csv example_predictions.csv \
        --output-dir outputs/example_run
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path


SEED = 20260427
N_PER_SITE = 150
SITES = [("Hospital-A", "2021"), ("Hospital-B", "2022")]
ORGANISM = "Escherichia coli"
DRUGS = [
    "Ciprofloxacin",
    "Norfloxacin",
    "Ceftriaxone",
    "Amoxicillin-Clavulanic acid",
]
OUT = Path(__file__).resolve().parents[1] / "example_predictions.csv"
FIELDNAMES = ["isolate_id", "site", "year", "organism", "drug", "label", "prob", "model_name"]


def sigmoid(x: float) -> float:
    x = max(-20.0, min(20.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def clamp(x: float, lo: float = 0.02, hi: float = 0.98) -> float:
    return max(lo, min(hi, x))


def generate_isolate(rng: random.Random, site: str, year: str, idx: int) -> list[dict]:
    isolate_id = f"ISO_{site[:3].upper()}_{idx:04d}"

    # Fluoroquinolone block: ~30% cipro-resistant; norfloxacin almost always co-resistant
    cipro_r = int(rng.random() < 0.30)
    norflox_r = cipro_r if rng.random() < 0.92 else 1 - cipro_r

    # Beta-lactam/ESBL block: ~20% ceftriaxone-resistant; amox-clav highly co-resistant.
    # Weak cross-block association (OR ~2).
    cef_prev = 0.28 if cipro_r else 0.15
    cef_r = int(rng.random() < cef_prev)
    amox_r = cef_r if rng.random() < 0.80 else int(rng.random() < 0.08)

    labels = {
        "Ciprofloxacin": cipro_r,
        "Norfloxacin": norflox_r,
        "Ceftriaxone": cef_r,
        "Amoxicillin-Clavulanic acid": amox_r,
    }

    # --- Model probabilities ---
    # Ciprofloxacin: drug-specific — discriminates even within co-resistance strata.
    cipro_prob = clamp(sigmoid(2.8 * cipro_r - 1.2 + rng.gauss(0, 0.55)))

    # Norfloxacin: also drug-specific (same chromosomal mechanism as cipro).
    norflox_prob = clamp(sigmoid(2.5 * norflox_r - 1.1 + rng.gauss(0, 0.55)))

    # Ceftriaxone: background-driven — score tracks cipro/norfloxacin resistance,
    # not the actual ceftriaxone label.  Collapses within matching strata.
    cef_prob = clamp(sigmoid(2.3 * cipro_r + 1.9 * norflox_r - 1.6 + rng.gauss(0, 0.6)))

    # Amoxicillin-Clavulanic acid: similarly background-driven.
    amox_prob = clamp(sigmoid(1.9 * cipro_r + 1.6 * norflox_r - 1.3 + rng.gauss(0, 0.65)))

    probs = {
        "Ciprofloxacin": cipro_prob,
        "Norfloxacin": norflox_prob,
        "Ceftriaxone": cef_prob,
        "Amoxicillin-Clavulanic acid": amox_prob,
    }

    rows = []
    for drug in DRUGS:
        rows.append({
            "isolate_id": isolate_id,
            "site": site,
            "year": year,
            "organism": ORGANISM,
            "drug": drug,
            "label": labels[drug],
            "prob": round(probs[drug], 4),
            "model_name": "SyntheticModel-v1",
        })
    return rows


def main() -> None:
    rng = random.Random(SEED)
    rows: list[dict] = []

    for site, year in SITES:
        for i in range(N_PER_SITE):
            rows.extend(generate_isolate(rng, site, year, i + 1))

    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    n_isolates = N_PER_SITE * len(SITES)
    print(f"Wrote {len(rows)} rows ({n_isolates} isolates x {len(DRUGS)} drugs) -> {OUT}")
    print()
    print("Expected audit outcomes:")
    print("  Ciprofloxacin        -- survives background matching (drug-specific signal)")
    print("  Norfloxacin          -- survives (same fluoroquinolone block)")
    print("  Ceftriaxone          -- collapses after matching (background-driven)")
    print("  Amoxicillin-Clav     -- collapses after matching (background-driven)")
    print()
    print("Run the audit with:")
    print(f"  python run_background_audit_framework.py \\")
    print(f"    --predictions-csv {OUT.name} \\")
    print(f"    --output-dir outputs/example_run")


if __name__ == "__main__":
    main()
