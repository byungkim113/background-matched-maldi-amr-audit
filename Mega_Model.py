"""
Mega_Model.py — consolidated MALDI-TOF AMR multi-drug pipeline.

Notebook default:
  clinical4 diagnostic = E.coli/Cipro, E.coli/Amox-Clav, S.aureus/Oxacillin,
                         S.epidermidis/Erythromycin
  with primary early stopping, all seeds, no prevalence shift, and no BN adaptation

Also available:
  run14          = original 2-pair overlap set
  clinical4      = clinical5 without S.aureus/Penicillin
  clinical5      = clinical4 plus S.aureus/Penicillin
  gram_negative6 = older 6-pair gram-negative-heavy expansion
  ecoli_drug_panel = E.coli-only drug-chemistry panel for fingerprint FiLM

Best-of fixes included:
  1. Per-pair Youden thresholds with finite/clipped ROC threshold handling
  2. Capped prevalence prior shifting with source-only calibration
  3. Macro-AUC early stopping across active pairs
  4. DANN gradient reversal lambda applied exactly once
  5. Pair-and-label balanced sampling plus within-pair mixup
  6. MIN_R_TRAIN guard reserves MIN_R_VAL resistant validation samples
  7. Empty same-site test split guard
  8. MAE pretraining excludes same-site test-year samples
 9. Run 14 overlap summary for direct E.coli/Cipro + S.aureus/Oxacillin comparison
10. Paper-ready core result tables plus occlusion saliency exports
11. Mechanism/detectability framing for chromosomal-structural vs mobile-heterogeneous resistance
 12. Optional drug-structure conditioning for E.coli panel experiments

Usage:
  python Mega_Model.py --pair-profile clinical4 --mae-epochs 30 \
      --experiment exp_clinical4_primary_noprev_nobn_allseeds_mae30
  python Mega_Model.py --pair-profile clinical5 --experiment exp_mega_clinical5
  python Mega_Model.py --pair-profile run14 --experiment exp_run14_repro
  python Mega_Model.py --pair-profile gram_negative6 --experiment exp_gn6
  python Mega_Model.py --pair-profile ecoli_drug_panel --drug-conditioning morgan \
      --with-leave-one-drug-out --experiment exp_ecoli_panel_morgan_loodo
"""

import argparse
import contextlib
import hashlib
import json
import logging
import math
import os
import pickle
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
try:
    import scipy.optimize as scipy_optimize
except ImportError:
    scipy_optimize = None
try:
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
    sklearn_metrics_available = True
except ImportError:
    sklearn_metrics_available = False

    def roc_auc_score(labels, probs):
        labels = np.asarray(labels).astype(int)
        probs = np.asarray(probs, dtype=float)
        n_pos = int(labels.sum())
        n_neg = len(labels) - n_pos
        if n_pos == 0 or n_neg == 0:
            raise ValueError("roc_auc_score is undefined with one class")
        order = np.argsort(probs, kind="mergesort")
        sorted_probs = probs[order]
        ranks = np.empty(len(probs), dtype=float)
        i = 0
        while i < len(probs):
            j = i
            while j + 1 < len(probs) and sorted_probs[j + 1] == sorted_probs[i]:
                j += 1
            avg_rank = 0.5 * (i + j) + 1.0
            ranks[order[i:j + 1]] = avg_rank
            i = j + 1
        pos_rank_sum = float(ranks[labels == 1].sum())
        return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)

    def average_precision_score(labels, probs):
        labels = np.asarray(labels).astype(int)
        probs = np.asarray(probs, dtype=float)
        n_pos = int(labels.sum())
        if n_pos == 0:
            return 0.0
        order = np.argsort(-probs, kind="mergesort")
        y = labels[order]
        tp = np.cumsum(y)
        precision = tp / (np.arange(len(y)) + 1)
        return float(precision[y == 1].sum() / n_pos)

    def roc_curve(labels, probs):
        labels = np.asarray(labels).astype(int)
        probs = np.asarray(probs, dtype=float)
        order = np.argsort(-probs, kind="mergesort")
        y = labels[order]
        s = probs[order]
        distinct = np.where(np.diff(s))[0]
        idx = np.r_[distinct, len(y) - 1]
        tps = np.cumsum(y)[idx]
        fps = 1 + idx - tps
        tps = np.r_[0, tps]
        fps = np.r_[0, fps]
        thresholds = np.r_[np.inf, s[idx]]
        n_pos = tps[-1]
        n_neg = fps[-1]
        tpr = tps / n_pos if n_pos > 0 else np.zeros_like(tps, dtype=float)
        fpr = fps / n_neg if n_neg > 0 else np.zeros_like(fps, dtype=float)
        return fpr, tpr, thresholds

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler, ConcatDataset

# ═══════════════════════════════════════════════════════════════════════════════
# §1  ORGANISM-DRUG PAIRS & PREVALENCE TABLE
# ═══════════════════════════════════════════════════════════════════════════════

PAIR_PROFILES = {
    "run14": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Staphylococcus aureus", "Oxacillin"),
    ],
    "clinical5": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
        ("Staphylococcus aureus", "Oxacillin"),
        ("Staphylococcus aureus", "Penicillin"),
        ("Staphylococcus epidermidis", "Erythromycin"),
    ],
    "saureus_panel": [
        # Main second-organism panel for background-matched Oxacillin audit.
        # The audit interprets Oxacillin against same-organism non-focal AST labels.
        # Pairs with insufficient source resistant isolates are dropped automatically.
        ("Staphylococcus aureus", "Oxacillin"),
        ("Staphylococcus aureus", "Penicillin"),
        ("Staphylococcus aureus", "Ciprofloxacin"),
        ("Staphylococcus aureus", "Erythromycin"),
        ("Staphylococcus aureus", "Clindamycin"),
        ("Staphylococcus aureus", "Gentamicin"),
        ("Staphylococcus aureus", "Fusidic acid"),
    ],
    "clinical4": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
        ("Staphylococcus aureus", "Oxacillin"),
        ("Staphylococcus epidermidis", "Erythromycin"),
    ],
    "gram_negative6": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Escherichia coli", "Ceftriaxone"),
        ("Escherichia coli", "Ceftazidime"),
        ("Klebsiella pneumoniae", "Ciprofloxacin"),
        ("Klebsiella pneumoniae", "Ceftriaxone"),
        ("Staphylococcus aureus", "Oxacillin"),
    ],
    "ecoli_drug_panel": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
        ("Escherichia coli", "Ceftriaxone"),
        ("Escherichia coli", "Ceftazidime"),
        ("Escherichia coli", "Piperacillin-Tazobactam"),
        ("Escherichia coli", "Gentamicin"),
        ("Escherichia coli", "Trimethoprim-Sulfamethoxazole"),
    ],
    "ecoli_mechanism6": [
        ("Escherichia coli", "Ciprofloxacin"),
        ("Escherichia coli", "Norfloxacin"),
        ("Escherichia coli", "Amoxicillin-Clavulanic acid"),
        ("Escherichia coli", "Ceftriaxone"),
        ("Escherichia coli", "Ceftazidime"),
        ("Escherichia coli", "Cefepime"),
    ],
}

RUN14_OVERLAP_PAIRS = list(PAIR_PROFILES["run14"])
RUN14_PAIRS = PAIR_PROFILES["run14"]
MULTIDRUG_PAIRS = PAIR_PROFILES["gram_negative6"]

# Per-(organism, drug) per-site resistance prevalence from published antibiograms.
# Source: Weis 2022, Viollier AG, ECDC. None = unknown → no prevalence shift.
ALL_SITE_PREVALENCE = {
    ("Escherichia coli",      "Ciprofloxacin"):                 {"DRIAMS-A": 0.28, "DRIAMS-B": 0.27, "DRIAMS-C": 0.21, "DRIAMS-D": 0.19},
    ("Escherichia coli",      "Amoxicillin-Clavulanic acid"):   {"DRIAMS-A": 0.24, "DRIAMS-B": 0.31, "DRIAMS-C": 0.25, "DRIAMS-D": 0.17},
    ("Staphylococcus aureus", "Oxacillin"):                      {"DRIAMS-A": 0.19, "DRIAMS-B": 0.06, "DRIAMS-C": 0.06, "DRIAMS-D": None},
    ("Staphylococcus aureus", "Penicillin"):                     {"DRIAMS-A": 0.73, "DRIAMS-B": 0.65, "DRIAMS-C": 0.74, "DRIAMS-D": 0.88},
    ("Staphylococcus epidermidis", "Erythromycin"):              {"DRIAMS-A": 0.61, "DRIAMS-B": 0.45, "DRIAMS-C": None, "DRIAMS-D": 0.68},
    ("Escherichia coli",      "Ceftriaxone"):                    {"DRIAMS-A": 0.14, "DRIAMS-B": 0.14, "DRIAMS-C": 0.13, "DRIAMS-D": None},
    ("Escherichia coli",      "Ceftazidime"):                    {"DRIAMS-A": 0.08, "DRIAMS-B": 0.08, "DRIAMS-C": 0.08, "DRIAMS-D": None},
    ("Escherichia coli",      "Norfloxacin"):                    {"DRIAMS-A": None, "DRIAMS-B": None, "DRIAMS-C": None, "DRIAMS-D": None},
    ("Escherichia coli",      "Cefepime"):                       {"DRIAMS-A": None, "DRIAMS-B": None, "DRIAMS-C": None, "DRIAMS-D": None},
    ("Escherichia coli",      "Piperacillin-Tazobactam"):         {"DRIAMS-A": None, "DRIAMS-B": None, "DRIAMS-C": None, "DRIAMS-D": None},
    ("Escherichia coli",      "Gentamicin"):                      {"DRIAMS-A": None, "DRIAMS-B": None, "DRIAMS-C": None, "DRIAMS-D": None},
    ("Escherichia coli",      "Trimethoprim-Sulfamethoxazole"):   {"DRIAMS-A": None, "DRIAMS-B": None, "DRIAMS-C": None, "DRIAMS-D": None},
    ("Klebsiella pneumoniae", "Ciprofloxacin"):                  {"DRIAMS-A": 0.25, "DRIAMS-B": 0.25, "DRIAMS-C": 0.25, "DRIAMS-D": None},
    ("Klebsiella pneumoniae", "Ceftriaxone"):                    {"DRIAMS-A": 0.25, "DRIAMS-B": 0.25, "DRIAMS-C": 0.25, "DRIAMS-D": None},
}

# Curated drug structure strings for optional drug-aware FiLM. When RDKit is
# available these are converted to Morgan fingerprints; otherwise a deterministic
# hashed-SMILES fallback keeps the pipeline runnable on plain Kaggle images.
DRUG_SMILES = {
    "Ciprofloxacin": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    "Norfloxacin": "CCn1cc(C(=O)O)c(=O)c2cc(F)c(N3CCNCC3)cc21",
    "Amoxicillin-Clavulanic acid": (
        "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O."
        "O=C(O)[C@@H]1C[C@H]2C(=O)N2C1=O"
    ),
    "Ceftriaxone": "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(CSc3nnnn3C)CS[C@H]12)c4csc(N)n4",
    "Ceftazidime": "CC(C)(O/N=C(/C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(C[N+]3(C)CCCC3)CS[C@H]12)c4csc(N)n4)C(=O)O",
    "Cefepime": "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(C[N+]3(C)CCCC3)CS[C@H]12)c4csc(N)n4",
    "Piperacillin-Tazobactam": (
        "CCN1CCN(C(=O)N[C@@H](C(=O)N2[C@H](C(=O)O)C(C)(C)S[C@@H]2C(=O)O)c3ccccc3)CC1."
        "CC1(C)S[C@@H]2[C@H](NC(=O)C3N4C=CC(C4=O)=C3C(=O)O)C(=O)N2[C@H]1C(=O)O"
    ),
    "Gentamicin": "CN[C@@H]1[C@@H](O)[C@H](O)[C@@H](N)C(O)[C@H]1O",
    "Trimethoprim-Sulfamethoxazole": (
        "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC."
        "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    ),
    "Oxacillin": "Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]3C(=O)N4[C@H](C(=O)O)C(C)(C)S[C@H]34",
    "Penicillin": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "Erythromycin": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@@H](C)C[C@@H](C)[C@H](O)[C@@H](C)C(=O)O1",
}

MECHANISM_DETECTABILITY = {
    ("Escherichia coli", "Ciprofloxacin"): {
        "expected_detectability": "more_detectable",
        "mechanism_class": "chromosomal_structural",
        "mechanism_label": "More detectable: chromosomal target mutations",
        "mechanism_rationale": (
            "Fluoroquinolone resistance is often driven by chromosomal gyrA/parC "
            "target mutations and associated stress-response phenotypes, which are "
            "more likely to leave stable MALDI-detectable proteomic structure."
        ),
    },
    ("Escherichia coli", "Norfloxacin"): {
        "expected_detectability": "more_detectable",
        "mechanism_class": "chromosomal_structural",
        "mechanism_label": "More detectable: fluoroquinolone block",
        "mechanism_rationale": (
            "Norfloxacin labels are highly co-resistant with ciprofloxacin in "
            "DRIAMS-A; treat this as a fluoroquinolone-block sanity check rather "
            "than an independent phenotype."
        ),
    },
    ("Escherichia coli", "Ceftriaxone"): {
        "expected_detectability": "less_detectable",
        "mechanism_class": "mobile_heterogeneous",
        "mechanism_label": "Less detectable: ESBL/AmpC beta-lactam block",
        "mechanism_rationale": (
            "Third/fourth generation cephalosporin resistance is often mediated "
            "by ESBL/AmpC and co-resistance blocks; the signal is expected to be "
            "class-structured but less phenotype-specific than fluoroquinolones."
        ),
    },
    ("Escherichia coli", "Ceftazidime"): {
        "expected_detectability": "less_detectable",
        "mechanism_class": "mobile_heterogeneous",
        "mechanism_label": "Less detectable: ESBL/AmpC beta-lactam block",
        "mechanism_rationale": (
            "Ceftazidime is part of the E. coli ESBL/AmpC co-resistance block, "
            "so transfer should be interpreted at the mechanism-block level."
        ),
    },
    ("Escherichia coli", "Cefepime"): {
        "expected_detectability": "less_detectable",
        "mechanism_class": "mobile_heterogeneous",
        "mechanism_label": "Less detectable: ESBL/AmpC beta-lactam block",
        "mechanism_rationale": (
            "Cefepime clusters strongly with ceftriaxone/ceftazidime in source "
            "co-resistance audits, making it useful for ESBL/AmpC block testing."
        ),
    },
    ("Staphylococcus aureus", "Oxacillin"): {
        "expected_detectability": "more_detectable",
        "mechanism_class": "chromosomal_structural",
        "mechanism_label": "More detectable: mecA/PBP2a structural phenotype",
        "mechanism_rationale": (
            "Methicillin resistance is tied to mecA/PBP2a and cell-wall remodeling, "
            "a larger structural phenotype that is plausibly more visible in spectra."
        ),
    },
    ("Escherichia coli", "Amoxicillin-Clavulanic acid"): {
        "expected_detectability": "less_detectable",
        "mechanism_class": "mobile_heterogeneous",
        "mechanism_label": "Less detectable: heterogeneous beta-lactamase context",
        "mechanism_rationale": (
            "Amoxicillin-clavulanate resistance can arise through plasmid beta-lactamases, "
            "porin/regulatory changes, inhibitor-resistant enzymes, and strain context, "
            "so the proteomic signal is expected to be less consistent."
        ),
    },
    ("Staphylococcus epidermidis", "Erythromycin"): {
        "expected_detectability": "less_detectable",
        "mechanism_class": "mobile_heterogeneous",
        "mechanism_label": "Less detectable: mobile ribosomal methylation/efflux",
        "mechanism_rationale": (
            "Macrolide resistance is often mediated by mobile erm/msr determinants "
            "or efflux/ribosomal methylation, which may not produce a stable core "
            "proteome signature across sites."
        ),
    },
}

MECHANISM_GROUP_LABELS = {
    "more_detectable": "More detectable (chromosomal/structural)",
    "less_detectable": "Less detectable (mobile/heterogeneous)",
    "unknown": "Unknown / not pre-specified",
}

SITE_CONTEXT = {
    "DRIAMS-A": {
        "site_context": "source_inpatient_training_site",
        "interpretation": "Source training hospital.",
    },
    "DRIAMS-B": {
        "site_context": "external_hospital_site",
        "interpretation": "External hospital distribution shift.",
    },
    "DRIAMS-C": {
        "site_context": "external_hospital_site",
        "interpretation": "External hospital distribution shift.",
    },
    "DRIAMS-D": {
        "site_context": "community_acquired_enriched",
        "interpretation": (
            "Largest degradation at DRIAMS-D is interpreted as consistent with a "
            "community-acquired distribution shift, not just generic model failure."
        ),
    },
}

# ═══════════════════════════════════════════════════════════════════════════════
# §2  CONFIG CONSTANTS (mode-independent)
# ═══════════════════════════════════════════════════════════════════════════════

DATA_ROOT   = "/kaggle/input/datasets/drscarlat/driams"
TRAIN_SITE  = "DRIAMS-A"
TEST_SITES  = ["DRIAMS-B", "DRIAMS-C", "DRIAMS-D"]
TEST_YEAR   = "2018"
VAL_HOLDOUT = 0.10

EPOCHS         = 60
BATCH_SIZE     = 64
LR             = 1e-4
PATIENCE       = 10
GRL_LAMBDA_MAX = 0.3
N_SITES        = 4
N_SEEDS        = 5
TTA_PASSES     = 15
MIXUP_ALPHA    = 0.2
CHECK_INTERVAL = 5
LABEL_SMOOTH   = 0.05

MAE_EPOCHS     = 30
MAE_MASK_RATIO = 0.75
MAE_PATCH_SIZE = 50
MAE_LR         = 1e-3

N_BINS   = 6000
PROB_EPS = 1e-6
MAX_PREVALENCE_ODDS_SHIFT = 2.0
DEVICE   = torch.device("cuda" if torch.cuda.is_available() else "cpu")

DRUG_FINGERPRINT_DIM = 2048
DRUG_CONDITIONING_CHOICES = ("task_id", "drug_id", "morgan", "drug_id+morgan")
DRUG_CONDITIONING_EMBED_DIM = 32
DRUG_FINGERPRINT_EMBED_DIM = 64

MIN_R_TRAIN  = 20
MIN_R_VAL    = 5
MIN_TOTAL_SOURCE_SAMPLES_WARN = 1000
EPOCHS_ABL   = 40
PATIENCE_ABL = 8
N_SEEDS_ABL  = 3

SEED_POLICY = "topk"      # "all", "threshold", or "topk"
TOP_K_SEEDS = 3
MIN_SEED_MACRO_AUC = 0.65
NEAR_CHANCE_MACRO_AUC = 0.55

SALIENCY_WINDOW = 200
SALIENCY_STRIDE = 100
SALIENCY_MAX_SAMPLES_PER_PAIR = 50
SALIENCY_TOP_K = 10
SALIENCY_STABILITY_TOP_K = 5
SALIENCY_STABILITY_TOP_KS = (5, 10)
SALIENCY_EXCLUDE_LOW_BINS_BELOW = 200
SALIENCY_BROAD_WINDOW = 500
SALIENCY_NULL_N = 250

RANDOM_CV_HOLDOUT = 0.20
RANDOM_CV_SEED = 2026
STAT_BOOTSTRAP_N = 500
STAT_BOOTSTRAP_SEED = 2026

# Mutable globals — set by init_config() after argparse
ORGANISM_DRUG_PAIRS: list = []
N_ORGANISMS:         int  = 0
PRIMARY_PAIR_IDX:    int  = 0

# ═══════════════════════════════════════════════════════════════════════════════
# NOTEBOOK / KAGGLE CONFIG — edit these when running as a notebook cell
# (ignored when running as a CLI script)
# ═══════════════════════════════════════════════════════════════════════════════
NOTEBOOK_PAIR_PROFILE = "saureus_panel"        # run14, clinical4, clinical5, saureus_panel, gram_negative6, ecoli_mechanism6
NOTEBOOK_MODE        = None                    # optional legacy alias: run14 or multidrug
NOTEBOOK_EXPERIMENT  = "exp_saureus_panel_oxa_background_mae30"
NOTEBOOK_CKPT_DIR    = None                    # None = auto (output/experiment/models)
NOTEBOOK_CALIB_JSON  = None                    # path to saved calibration JSON, or None
NOTEBOOK_DATA_ROOT   = None                    # None = use DATA_ROOT above
NOTEBOOK_OUTPUT_DIR  = "/kaggle/working/runs"
NOTEBOOK_WITH_ABLATION = False                 # True = train no_mae/no_dann/no_film (slow)
NOTEBOOK_NO_LGBM     = True                    # True = skip LightGBM baselines; keep Sa/Oxa audit run faster
NOTEBOOK_NO_BN_ADAPT = True                    # True = disable BN adaptation on ext sites
NOTEBOOK_EARLY_STOP  = "macro"                 # "macro" or "primary"
NOTEBOOK_SEED_POLICY = "all"                   # "all", "threshold", or "topk"
NOTEBOOK_TOP_K_SEEDS = TOP_K_SEEDS
NOTEBOOK_MIN_SEED_MACRO_AUC = MIN_SEED_MACRO_AUC
NOTEBOOK_PREVALENCE_SHIFT = "none"             # "capped" or "none"
NOTEBOOK_STRICT_MAE_SOURCE_ONLY = False        # True = MAE sees source train/val only
NOTEBOOK_NO_SALIENCY = True                    # True = skip occlusion saliency export; audit does not need saliency
NOTEBOOK_MAE_EPOCHS = 30                       # MAE-30 performance run default
NOTEBOOK_WITH_RANDOM_CV = False                # True = run LightGBM random-CV inflation diagnostic
NOTEBOOK_DRUG_CONDITIONING = "task_id"         # task_id, drug_id, morgan, drug_id+morgan
NOTEBOOK_WITH_LEAVE_ONE_DRUG_OUT = False       # E.coli panel zero-shot drug transfer diagnostic


def init_config(pair_profile: str):
    global ORGANISM_DRUG_PAIRS, N_ORGANISMS, PRIMARY_PAIR_IDX
    if pair_profile not in PAIR_PROFILES:
        raise ValueError(f"Unknown pair profile: {pair_profile}")
    ORGANISM_DRUG_PAIRS = list(PAIR_PROFILES[pair_profile])
    N_ORGANISMS      = len(ORGANISM_DRUG_PAIRS)
    PRIMARY_PAIR_IDX = 0
    print(f"Pair profile={pair_profile}  pairs={N_ORGANISMS}  device={DEVICE}")
    for i, (org, drug) in enumerate(ORGANISM_DRUG_PAIRS):
        print(f"  [{i}] {org} / {drug}" + (" ← primary" if i == PRIMARY_PAIR_IDX else ""))


def _hashed_smiles_fingerprint(smiles, n_bits=DRUG_FINGERPRINT_DIM):
    """Deterministic fallback fingerprint for environments without RDKit."""
    bits = np.zeros(n_bits, dtype=np.float32)
    text = str(smiles or "")
    if not text:
        return bits
    for width in (2, 3, 4, 5):
        if len(text) < width:
            continue
        for i in range(len(text) - width + 1):
            token = f"{width}:{text[i:i + width]}"
            digest = hashlib.sha256(token.encode("utf-8")).hexdigest()
            bits[int(digest[:8], 16) % n_bits] = 1.0
    return bits


def smiles_to_morgan_fingerprint(smiles, n_bits=DRUG_FINGERPRINT_DIM):
    """
    Convert a SMILES string to a Morgan fingerprint when RDKit is available.
    Falls back to hashed SMILES n-grams so the pipeline remains portable.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit import DataStructs
    except ImportError:
        return _hashed_smiles_fingerprint(smiles, n_bits)

    mol = Chem.MolFromSmiles(str(smiles or ""))
    if mol is None:
        return _hashed_smiles_fingerprint(smiles, n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.astype(np.float32)


def build_drug_feature_matrix(pairs=None, conditioning="morgan"):
    """
    Build per-task drug features indexed by the existing pair/task id.
    task_id mode keeps the legacy architecture and does not need this matrix.
    """
    pairs = list(pairs if pairs is not None else ORGANISM_DRUG_PAIRS)
    if conditioning not in DRUG_CONDITIONING_CHOICES:
        raise ValueError(f"Unknown drug conditioning mode: {conditioning}")
    if conditioning == "task_id":
        return None

    rows = []
    for idx, (_, drug) in enumerate(pairs):
        parts = []
        if "morgan" in conditioning:
            parts.append(smiles_to_morgan_fingerprint(
                DRUG_SMILES.get(drug, drug), DRUG_FINGERPRINT_DIM))
        rows.append(np.concatenate(parts).astype(np.float32) if parts
                    else np.zeros(0, dtype=np.float32))
    return np.stack(rows, axis=0) if rows else np.zeros((0, 0), dtype=np.float32)


# ═══════════════════════════════════════════════════════════════════════════════
# §3  SAFETY GATE  (Validation Protocol §3.4)
# ═══════════════════════════════════════════════════════════════════════════════

_EXTERNAL_EVAL_OPEN = False


def open_external_eval():
    global _EXTERNAL_EVAL_OPEN
    _EXTERNAL_EVAL_OPEN = True


def assert_source_only(fn_name: str):
    if _EXTERNAL_EVAL_OPEN:
        raise RuntimeError(
            f"SAFETY VIOLATION: {fn_name}() called after open_external_eval(). "
            "Calibration must use source validation data only."
        )


# ═══════════════════════════════════════════════════════════════════════════════
# §4  DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_labels(root, site, organism, drug):
    id_root = Path(root) / site / "id"
    frames = []
    for year_dir in sorted(id_root.iterdir()):
        clean = year_dir / f"{year_dir.name}_clean.csv"
        if clean.exists():
            frames.append(pd.read_csv(clean, low_memory=False))
    if not frames:
        return {}
    df = pd.concat(frames, ignore_index=True)
    species_col = next((c for c in df.columns
                        if "species" in c.lower() or "organism" in c.lower()), None)
    if species_col:
        df = df[df[species_col].str.contains(organism, case=False, na=False)]
    if df.empty:
        return {}
    if drug not in df.columns:
        candidates = [c for c in df.columns if drug.lower() in c.lower()]
        if not candidates:
            return {}
        drug_col = candidates[0]
    else:
        drug_col = drug
    id_col = next((c for c in df.columns
                   if c.lower() in ("sample_id", "code", "uuid", "id")), df.columns[0])
    remap = {"S": 0, "R": 1, "I": -1, "-": -1}
    return {str(row[id_col]).strip(): remap.get(str(row[drug_col]).strip(), -1)
            for _, row in df.iterrows()}


def build_spectrum_index(root, site):
    spec_root = Path(root) / site / "binned_6000"
    if not spec_root.exists():
        return None
    all_files = sorted(spec_root.rglob("*.txt"))
    stem_to_paths: dict = {}
    for txt in all_files:
        stem_to_paths.setdefault(txt.stem, []).append(txt)
    return all_files, stem_to_paths


def samples_from_label_map(labels, org_id, spectrum_index):
    _, stem_to_paths = spectrum_index
    out = []
    for uid in sorted(labels):
        if labels[uid] == -1:
            continue
        for txt in stem_to_paths.get(uid, []):
            out.append((str(txt), labels[uid], org_id, txt.parent.name))
    return out


def load_all_organisms(root, site, active_pairs, spectrum_index=None):
    """Returns (path, label, org_id, year) 4-tuples for active_pairs only."""
    print(f"  {site}:")
    if spectrum_index is None:
        spectrum_index = build_spectrum_index(root, site)
    if spectrum_index is None:
        print("    binned_6000 not found — skipping")
        return []
    all_files, _ = spectrum_index
    print(f"    {len(all_files)} spectrum files")
    all_samples = []
    for org_id, organism, drug in active_pairs:
        labels = load_labels(root, site, organism, drug)
        if not labels:
            print(f"    [{org_id}] {organism[:22]:22s}/{drug}: no labels")
            continue
        org_samples = samples_from_label_map(labels, org_id, spectrum_index)
        if not org_samples:
            print(f"    [{org_id}] {organism[:22]:22s}/{drug}: n=0")
            continue
        r = sum(1 for _, l, _, _ in org_samples if l == 1)
        print(f"    [{org_id}] {organism[:22]:22s}/{drug}: n={len(org_samples)}  R={r}")
        all_samples.extend(org_samples)
    print(f"    total={len(all_samples)}")
    return all_samples


def required_pretest_resistant(min_r_train=MIN_R_TRAIN, min_r_val=MIN_R_VAL):
    return MIN_R_TRAIN + MIN_R_VAL if (min_r_train, min_r_val) == (MIN_R_TRAIN, MIN_R_VAL) else min_r_train + min_r_val


def has_enough_pretest_resistance(n_r_pretest, min_r_train=MIN_R_TRAIN,
                                  min_r_val=MIN_R_VAL):
    return n_r_pretest >= required_pretest_resistant(min_r_train, min_r_val)


def screen_active_pairs(root, pair_profile):
    """
    Apply MIN_R_TRAIN guard after reserving MIN_R_VAL resistant validation samples.
    In run14 profile, return all pairs unchanged so the exact baseline can be reproduced.
    """
    global PRIMARY_PAIR_IDX
    base = [(i, org, drug) for i, (org, drug) in enumerate(ORGANISM_DRUG_PAIRS)]
    if pair_profile == "run14":
        return base

    spectrum_index = build_spectrum_index(root, TRAIN_SITE)
    if spectrum_index is None:
        raise RuntimeError(f"{TRAIN_SITE}/binned_6000 not found.")

    active, dropped = [], []
    for org_id, organism, drug in base:
        lmap = load_labels(root, TRAIN_SITE, organism, drug)
        if not lmap:
            dropped.append((org_id, organism, drug, "no labels"))
            continue
        samps = samples_from_label_map(lmap, org_id, spectrum_index)
        n_r_pretest = sum(l for _, l, _, y in samps if y < TEST_YEAR)
        if not has_enough_pretest_resistance(n_r_pretest):
            needed = required_pretest_resistant()
            dropped.append((org_id, organism, drug,
                            f"n_r_pretest={n_r_pretest}<{needed} "
                            "(MIN_R_TRAIN + MIN_R_VAL)"))
        else:
            active.append((org_id, organism, drug))

    print(f"\nActive pairs ({len(active)}):")
    for oid, org, drug in active:
        print(f"  [{oid}] {org} / {drug}")
    if dropped:
        print(f"Dropped ({len(dropped)}):")
        for oid, org, drug, reason in dropped:
            print(f"  [{oid}] {org} / {drug}: {reason}")
    if not active:
        raise RuntimeError("No active pairs passed MIN_R_TRAIN guard.")

    active_ids = {oid for oid, _, _ in active}
    if PRIMARY_PAIR_IDX not in active_ids:
        PRIMARY_PAIR_IDX = active[0][0]
        print(f"  WARNING: primary pair dropped — falling back to [{PRIMARY_PAIR_IDX}]")
    return active


def split_train_val_by_pair_label(samples, val_holdout=VAL_HOLDOUT, min_r_val=MIN_R_VAL,
                                   random_state=42):
    """Stratify by (org_id, label) cell; guarantee min_r_val resistant val samples per pair."""
    rng = np.random.default_rng(random_state)
    grouped = defaultdict(list)
    for sample in samples:
        _, label, org_id = sample
        grouped[(org_id, label)].append(sample)

    train_s, val_s = [], []
    for (org_id, label), group in sorted(grouped.items()):
        if len(group) < 2:
            train_s.extend(group)
            continue
        indices = rng.permutation(len(group))
        n_val = max(1, int(round(len(group) * val_holdout)))
        if label == 1 and len(group) >= min_r_val:
            n_val = max(n_val, min_r_val)
        n_val = min(n_val, len(group) - 1)
        val_idx = set(indices[:n_val].tolist())
        for i, s in enumerate(group):
            (val_s if i in val_idx else train_s).append(s)

    rng.shuffle(train_s)
    rng.shuffle(val_s)
    return train_s, val_s


def make_split(all_samples, active_pairs):
    year_counts: dict = {}
    for _, _, _, y in all_samples:
        year_counts[y] = year_counts.get(y, 0) + 1
    print(f"  Year dist: { {k: year_counts[k] for k in sorted(year_counts)} }")

    pre_test = [(p, l, o) for p, l, o, y in all_samples if y < TEST_YEAR]
    test_s   = [(p, l, o) for p, l, o, y in all_samples if y == TEST_YEAR]
    if not pre_test:
        raise RuntimeError(f"No samples before {TEST_YEAR}.")

    train_s, val_s = split_train_val_by_pair_label(pre_test)

    # Warn pairs with few resistant val samples
    for oid, organism, drug in active_pairs:
        n_r = sum(1 for _, l, o in val_s if o == oid and l == 1)
        if n_r < MIN_R_VAL:
            print(f"  WARNING [{oid}] {organism[:22]:22s}/{drug}: only {n_r} R in val")

    r_tr  = sum(l for _, l, _ in train_s)
    r_val = sum(l for _, l, _ in val_s)
    r_te  = sum(l for _, l, _ in test_s)
    print(f"  Train={len(train_s)}(R={r_tr})  Val={len(val_s)}(R={r_val})  "
          f"Test({TEST_YEAR})={len(test_s)}(R={r_te})")
    return train_s, val_s, test_s


# ═══════════════════════════════════════════════════════════════════════════════
# §5  SPECTRUM LOADING + AUGMENTATION + MIXUP
# ═══════════════════════════════════════════════════════════════════════════════

def load_spectrum(path):
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"{path}: malformed spectrum array with shape {getattr(data, 'shape', None)}; "
            "expected two columns"
        )
    if data.shape[0] < N_BINS:
        raise ValueError(
            f"{path}: only {data.shape[0]} rows, expected at least {N_BINS}"
        )
    x = np.log1p(np.maximum(data[:N_BINS, 1], 0.0)).astype(np.float32)
    return (x - x.mean()) / (x.std() + 1e-8)


def load_mz_axis(path):
    """Load the m/z axis for saliency; map DRIAMS bin indices to approximate Da."""
    try:
        data = np.loadtxt(path, skiprows=1)
        if data.ndim == 2 and data.shape[0] >= N_BINS:
            mz = data[:N_BINS, 0].astype(float)
            if mz[0] < 10 and mz[-1] < N_BINS + 10:
                mz = 2000.0 + mz * 3.0
            return mz
    except Exception:
        pass
    return 2000.0 + np.arange(N_BINS, dtype=float) * 3.0


def augment(x):
    x = x.copy()
    if np.random.rand() < 0.7:
        x += np.random.normal(0, np.random.uniform(0.01, 0.08), x.shape).astype(np.float32)
    if np.random.rand() < 0.5:
        shift = int(np.random.randint(-10, 11))
        if shift != 0:
            shifted = np.zeros_like(x)
            if shift > 0:
                shifted[shift:] = x[:-shift]
            else:
                shifted[:shift] = x[-shift:]
            x = shifted
    if np.random.rand() < 0.5:
        x *= np.random.uniform(0.85, 1.15)
    if np.random.rand() < 0.4:
        t = np.linspace(0, 1, len(x))
        x += np.polyval(np.random.uniform(-0.05, 0.05, 3), t).astype(np.float32)
    return x.astype(np.float32)


def intraclass_mixup(x, y, org_id, alpha=MIXUP_ALPHA):
    """Mix only within (same label AND same organism) — fixes run14_core bug."""
    if alpha <= 0:
        return x, y
    mixed = x.clone()
    for i in range(x.size(0)):
        same = ((y == y[i]) & (org_id == org_id[i])).nonzero(as_tuple=True)[0]
        if len(same) > 1:
            j = same[torch.randint(len(same), (1,)).item()].item()
            if j != i:
                lam = float(np.random.beta(alpha, alpha))
                mixed[i] = lam * x[i] + (1 - lam) * x[j]
    return mixed, y


# ═══════════════════════════════════════════════════════════════════════════════
# §6  DATASETS + LOADERS
# ═══════════════════════════════════════════════════════════════════════════════

class MALDIDataset(Dataset):
    def __init__(self, samples, training=True):
        self.samples = samples; self.training = training

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        path, label, org_id = self.samples[idx]
        x = load_spectrum(path)
        if self.training:
            x = augment(x)
        return (torch.from_numpy(x).unsqueeze(0),
                torch.tensor(label,  dtype=torch.float32),
                torch.tensor(org_id, dtype=torch.long))


class DomainDataset(Dataset):
    def __init__(self, paths, domain_id):
        self.paths = paths; self.domain_id = domain_id

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        x = augment(load_spectrum(self.paths[idx]))
        return torch.from_numpy(x).unsqueeze(0), torch.tensor(self.domain_id, dtype=torch.long)


class MAEDataset(Dataset):
    """Unlabeled dataset used for MAE pretraining and BN adaptation."""
    def __init__(self, paths):
        self.paths = paths

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        try:
            x = augment(load_spectrum(self.paths[idx]))
        except Exception:
            x = np.zeros(N_BINS, dtype=np.float32)
        return torch.from_numpy(x).unsqueeze(0)


def _sample_weights_by_pair_label(samples):
    """Weight inversely by (org_id, label) cell count — correct for multi-pair."""
    cell_counts = Counter((org_id, label) for _, label, org_id in samples)
    return [1.0 / cell_counts[(org_id, label)] for _, label, org_id in samples]


def make_loader(samples, training=True):
    ds = MALDIDataset(samples, training)
    if training:
        weights = _sample_weights_by_pair_label(samples)
        sampler = WeightedRandomSampler(weights, len(weights))
        return DataLoader(ds, batch_size=BATCH_SIZE, sampler=sampler,
                          num_workers=2, pin_memory=True)
    return DataLoader(ds, batch_size=BATCH_SIZE, shuffle=False,
                      num_workers=2, pin_memory=True)


def make_tta_loader(samples):
    return DataLoader(MALDIDataset(samples, training=True), batch_size=BATCH_SIZE,
                      shuffle=False, num_workers=2, pin_memory=True)


def make_domain_loader(target_site_data):
    SITE_DOMAIN = {"DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3}
    ds_list, all_domain_ids = [], []
    for site, domain_id in SITE_DOMAIN.items():
        if site not in target_site_data or not target_site_data[site]:
            continue
        paths = [s[0] for s in target_site_data[site]]
        ds_list.append(DomainDataset(paths, domain_id))
        all_domain_ids.extend([domain_id] * len(paths))
    if not ds_list:
        return None
    counts  = np.bincount(all_domain_ids, minlength=N_SITES)
    weights = [1.0 / max(counts[d], 1) for d in all_domain_ids]
    sampler = WeightedRandomSampler(weights, len(weights))
    return DataLoader(ConcatDataset(ds_list), batch_size=BATCH_SIZE,
                      sampler=sampler, num_workers=2, pin_memory=True, drop_last=True)


# ═══════════════════════════════════════════════════════════════════════════════
# §7  GRADIENT REVERSAL
# ═══════════════════════════════════════════════════════════════════════════════

class _GradRev(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, lam):
        ctx.lam = lam
        return x.clone()

    @staticmethod
    def backward(ctx, grad):
        return -ctx.lam * grad, None


def grad_reverse(x, lam):
    return _GradRev.apply(x, lam)


# ═══════════════════════════════════════════════════════════════════════════════
# §8  ARCHITECTURE
# ═══════════════════════════════════════════════════════════════════════════════

class SEBlock(nn.Module):
    def __init__(self, channels, reduction=8):
        super().__init__()
        self.se = nn.Sequential(
            nn.AdaptiveAvgPool1d(1), nn.Flatten(),
            nn.Linear(channels, channels // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(channels // reduction, channels, bias=False),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return x * self.se(x).unsqueeze(-1)


class ResBlock(nn.Module):
    def __init__(self, ch_in, ch_out, stride=1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(ch_in, ch_out, 7, stride=stride, padding=3, bias=False),
            nn.BatchNorm1d(ch_out), nn.ReLU(inplace=True),
            nn.Conv1d(ch_out, ch_out, 7, padding=3, bias=False),
            nn.BatchNorm1d(ch_out),
        )
        self.se   = SEBlock(ch_out)
        self.skip = (nn.Sequential(
            nn.Conv1d(ch_in, ch_out, 1, stride=stride, bias=False),
            nn.BatchNorm1d(ch_out),
        ) if (stride != 1 or ch_in != ch_out) else nn.Identity())
        self.act = nn.ReLU(inplace=True)

    def forward(self, x):
        return self.act(self.se(self.net(x)) + self.skip(x))


class MAEWrapper(nn.Module):
    def __init__(self):
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv1d(1, 32, 15, stride=2, padding=7, bias=False),
            nn.BatchNorm1d(32), nn.ReLU(inplace=True),
            nn.MaxPool1d(3, stride=2, padding=1),
        )
        self.layers = nn.Sequential(
            ResBlock(32,  64,  stride=2),
            ResBlock(64,  128, stride=2),
            ResBlock(128, 256, stride=2),
        )
        self.gap     = nn.AdaptiveAvgPool1d(1)
        self.decoder = nn.Sequential(
            nn.Linear(256, 512), nn.ReLU(inplace=True), nn.Linear(512, N_BINS),
        )

    def encode(self, x):
        return self.gap(self.layers(self.stem(x))).flatten(1)

    def forward(self, x):
        return self.decoder(self.encode(x))


class FiLMLayer(nn.Module):
    def __init__(self, n_organisms, n_features):
        super().__init__()
        self.gamma = nn.Embedding(n_organisms, n_features)
        self.beta  = nn.Embedding(n_organisms, n_features)
        nn.init.ones_(self.gamma.weight)
        nn.init.zeros_(self.beta.weight)

    def forward(self, x, org_id):
        return self.gamma(org_id) * x + self.beta(org_id)


class DrugConditionedFiLMLayer(nn.Module):
    """
    FiLM layer conditioned on organism/task id plus optional drug chemistry.
    Used for E.coli panel experiments where a shared resistance head can score
    held-out drugs from a fingerprint rather than a trained per-drug head.
    """
    def __init__(self, n_tasks, n_features, drug_features,
                 conditioning="morgan",
                 id_embed_dim=DRUG_CONDITIONING_EMBED_DIM,
                 fp_embed_dim=DRUG_FINGERPRINT_EMBED_DIM):
        super().__init__()
        if conditioning not in DRUG_CONDITIONING_CHOICES or conditioning == "task_id":
            raise ValueError(f"DrugConditionedFiLMLayer needs a drug-aware mode, got {conditioning!r}")
        self.conditioning = conditioning
        self.task_embed = nn.Embedding(n_tasks, id_embed_dim)

        drug_features = torch.as_tensor(drug_features, dtype=torch.float32)
        if drug_features.ndim != 2 or drug_features.shape[0] != n_tasks:
            raise ValueError("drug_features must have shape [n_tasks, feature_dim]")
        self.register_buffer("drug_features", drug_features)

        cond_dim = id_embed_dim
        if "drug_id" in conditioning:
            self.drug_id_embed = nn.Embedding(n_tasks, id_embed_dim)
            cond_dim += id_embed_dim
        else:
            self.drug_id_embed = None

        if "morgan" in conditioning:
            self.fingerprint_proj = nn.Sequential(
                nn.Linear(drug_features.shape[1], fp_embed_dim),
                nn.ReLU(inplace=True),
                nn.Linear(fp_embed_dim, fp_embed_dim),
                nn.ReLU(inplace=True),
            )
            cond_dim += fp_embed_dim
        else:
            self.fingerprint_proj = None

        self.to_gamma = nn.Linear(cond_dim, n_features)
        self.to_beta  = nn.Linear(cond_dim, n_features)
        nn.init.zeros_(self.to_gamma.weight)
        nn.init.ones_(self.to_gamma.bias)
        nn.init.zeros_(self.to_beta.weight)
        nn.init.zeros_(self.to_beta.bias)

    def forward(self, x, org_id):
        cond = [self.task_embed(org_id)]
        if self.drug_id_embed is not None:
            cond.append(self.drug_id_embed(org_id))
        if self.fingerprint_proj is not None:
            cond.append(self.fingerprint_proj(self.drug_features[org_id]))
        z = torch.cat(cond, dim=1)
        return self.to_gamma(z) * x + self.to_beta(z)


class MALDICNN(nn.Module):
    def __init__(self, n_sites=N_SITES, n_organisms=None):
        super().__init__()
        n_organisms = n_organisms or N_ORGANISMS
        self.stem = nn.Sequential(
            nn.Conv1d(1, 32, 15, stride=2, padding=7, bias=False),
            nn.BatchNorm1d(32), nn.ReLU(inplace=True),
            nn.MaxPool1d(3, stride=2, padding=1),
        )
        self.layers = nn.Sequential(
            ResBlock(32,  64,  stride=2),
            ResBlock(64,  128, stride=2),
            ResBlock(128, 256, stride=2),
        )
        self.gap       = nn.AdaptiveAvgPool1d(1)
        self.drop      = nn.Dropout(0.3)
        self.film      = FiLMLayer(n_organisms, 256)
        self.res_heads = nn.ModuleList([nn.Linear(256, 1) for _ in range(n_organisms)])
        self.dom_fc    = nn.Linear(256, n_sites)

    def encode(self, x):
        return self.gap(self.layers(self.stem(x))).flatten(1)

    def forward(self, x, org_id):
        f      = self.drop(self.encode(x))
        f_film = self.film(f, org_id)
        all_logits = torch.stack([h(f_film) for h in self.res_heads], dim=1).squeeze(-1)
        return all_logits.gather(1, org_id.unsqueeze(1)).squeeze(1)

    def forward_dann(self, x, org_id, lam):
        f      = self.encode(x)
        f_drop = self.drop(f)
        f_film = self.film(f_drop, org_id)
        res    = torch.stack([h(f_film) for h in self.res_heads], dim=1).squeeze(-1)
        res    = res.gather(1, org_id.unsqueeze(1)).squeeze(1)
        dom    = self.dom_fc(grad_reverse(f, lam))
        return res, dom

    def forward_domain_only(self, x, lam):
        return self.dom_fc(grad_reverse(self.encode(x), lam))


class MALDICNNDrugConditioned(nn.Module):
    """Drug-aware FiLM model with a shared resistance head for zero-shot drugs."""
    def __init__(self, n_sites=N_SITES, n_organisms=None,
                 drug_features=None, drug_conditioning="morgan"):
        super().__init__()
        n_organisms = n_organisms or N_ORGANISMS
        if drug_features is None:
            drug_features = build_drug_feature_matrix(
                ORGANISM_DRUG_PAIRS, conditioning=drug_conditioning)
        self.stem = nn.Sequential(
            nn.Conv1d(1, 32, 15, stride=2, padding=7, bias=False),
            nn.BatchNorm1d(32), nn.ReLU(inplace=True),
            nn.MaxPool1d(3, stride=2, padding=1),
        )
        self.layers = nn.Sequential(
            ResBlock(32,  64,  stride=2),
            ResBlock(64,  128, stride=2),
            ResBlock(128, 256, stride=2),
        )
        self.gap      = nn.AdaptiveAvgPool1d(1)
        self.drop     = nn.Dropout(0.3)
        self.film     = DrugConditionedFiLMLayer(
            n_organisms, 256, drug_features, conditioning=drug_conditioning)
        self.res_head = nn.Linear(256, 1)
        self.dom_fc   = nn.Linear(256, n_sites)

    def encode(self, x):
        return self.gap(self.layers(self.stem(x))).flatten(1)

    def forward(self, x, org_id):
        f = self.drop(self.encode(x))
        return self.res_head(self.film(f, org_id)).squeeze(1)

    def forward_dann(self, x, org_id, lam):
        f = self.encode(x)
        res = self.res_head(self.film(self.drop(f), org_id)).squeeze(1)
        dom = self.dom_fc(grad_reverse(f, lam))
        return res, dom

    def forward_domain_only(self, x, lam):
        return self.dom_fc(grad_reverse(self.encode(x), lam))


def create_maldi_model(n_sites=N_SITES, n_organisms=None,
                       drug_conditioning="task_id"):
    if drug_conditioning == "task_id":
        return MALDICNN(n_sites=n_sites, n_organisms=n_organisms)
    drug_features = build_drug_feature_matrix(
        ORGANISM_DRUG_PAIRS, conditioning=drug_conditioning)
    return MALDICNNDrugConditioned(
        n_sites=n_sites,
        n_organisms=n_organisms,
        drug_features=drug_features,
        drug_conditioning=drug_conditioning)


class MALDICNNNoFiLM(nn.Module):
    """Ablation: organism one-hot concatenated instead of FiLM modulation."""
    def __init__(self, n_sites=N_SITES, n_organisms=None):
        super().__init__()
        n_organisms = n_organisms or N_ORGANISMS
        self.n_organisms = n_organisms
        self.stem = nn.Sequential(
            nn.Conv1d(1, 32, 15, stride=2, padding=7, bias=False),
            nn.BatchNorm1d(32), nn.ReLU(inplace=True),
            nn.MaxPool1d(3, stride=2, padding=1),
        )
        self.layers = nn.Sequential(
            ResBlock(32,  64,  stride=2),
            ResBlock(64,  128, stride=2),
            ResBlock(128, 256, stride=2),
        )
        self.gap      = nn.AdaptiveAvgPool1d(1)
        self.drop     = nn.Dropout(0.3)
        self.res_head = nn.Linear(256 + n_organisms, 1)
        self.dom_fc   = nn.Linear(256, n_sites)

    def encode(self, x):
        return self.gap(self.layers(self.stem(x))).flatten(1)

    def forward(self, x, org_id):
        f  = self.drop(self.encode(x))
        oh = F.one_hot(org_id, self.n_organisms).float()
        return self.res_head(torch.cat([f, oh], dim=1)).squeeze(1)

    def forward_dann(self, x, org_id, lam):
        f  = self.encode(x)
        fd = self.drop(f)
        oh = F.one_hot(org_id, self.n_organisms).float()
        res = self.res_head(torch.cat([fd, oh], dim=1)).squeeze(1)
        dom = self.dom_fc(grad_reverse(f, lam))
        return res, dom

    def forward_domain_only(self, x, lam):
        return self.dom_fc(grad_reverse(self.encode(x), lam))


# ═══════════════════════════════════════════════════════════════════════════════
# §9  MAE PRETRAINING
# ═══════════════════════════════════════════════════════════════════════════════

def pretrain_mae(mae_model, mae_paths):
    n_patches = N_BINS // MAE_PATCH_SIZE
    loader    = DataLoader(MAEDataset(mae_paths), batch_size=128, shuffle=True,
                           num_workers=2, pin_memory=True, drop_last=True)
    optimizer = torch.optim.AdamW(mae_model.parameters(), lr=MAE_LR, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=MAE_EPOCHS)
    mae_model.train()
    for epoch in range(1, MAE_EPOCHS + 1):
        loss_sum, n = 0.0, 0
        for x in loader:
            x = x.to(DEVICE)
            B = x.shape[0]
            n_masked = int(n_patches * MAE_MASK_RATIO)
            mask = torch.zeros(B, n_patches, dtype=torch.bool, device=DEVICE)
            masked_ids = torch.rand(B, n_patches, device=DEVICE).argsort(dim=1)[:, :n_masked]
            mask.scatter_(1, masked_ids, True)
            mask_exp = mask.repeat_interleave(MAE_PATCH_SIZE, dim=1).unsqueeze(1)
            recon    = mae_model(x * (~mask_exp).float())
            loss     = F.mse_loss(recon[mask_exp.squeeze(1)], x.squeeze(1)[mask_exp.squeeze(1)])
            optimizer.zero_grad(); loss.backward()
            nn.utils.clip_grad_norm_(mae_model.parameters(), 1.0)
            optimizer.step()
            loss_sum += loss.item() * B; n += B
        scheduler.step()
        if epoch % 5 == 0 or epoch == 1:
            print(f"  MAE Ep{epoch:>3}/{MAE_EPOCHS}  loss={loss_sum/n:.6f}", flush=True)


def load_mae_weights(mae_model, maldicnn):
    maldicnn.stem.load_state_dict(mae_model.stem.state_dict())
    maldicnn.layers.load_state_dict(mae_model.layers.state_dict())


# ═══════════════════════════════════════════════════════════════════════════════
# §10  TRAINING UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def run_epoch_dann(model, src_loader, dom_loader, optimizer,
                   res_criterion, dom_criterion, lam):
    model.train()
    loss_sum, correct, n = 0.0, 0, 0
    dom_iter = iter(dom_loader)
    for x_s, y_s, org_s in src_loader:
        try:
            x_d, y_d = next(dom_iter)
        except StopIteration:
            dom_iter = iter(dom_loader)
            x_d, y_d = next(dom_iter)
        x_s, y_s, org_s = x_s.to(DEVICE), y_s.to(DEVICE), org_s.to(DEVICE)
        x_d, y_d        = x_d.to(DEVICE), y_d.to(DEVICE)
        x_s, y_s = intraclass_mixup(x_s, y_s, org_s)
        optimizer.zero_grad()
        res_out, dom_out_s = model.forward_dann(x_s, org_s, lam)
        y_smooth           = y_s * (1 - LABEL_SMOOTH) + (1 - y_s) * LABEL_SMOOTH
        res_loss           = res_criterion(res_out, y_smooth)
        src_dom_labels     = torch.zeros(len(x_s), dtype=torch.long, device=DEVICE)
        dom_loss_s         = dom_criterion(dom_out_s, src_dom_labels)
        dom_loss_t         = dom_criterion(model.forward_domain_only(x_d, lam), y_d)
        loss = res_loss + dom_loss_s + dom_loss_t
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        loss_sum += res_loss.item() * len(y_s)
        correct  += ((res_out > 0) == y_s.bool()).sum().item()
        n        += len(y_s)
    return loss_sum / n, correct / n


def _run_epoch_no_dann(model, loader, optimizer, res_criterion):
    model.train()
    loss_sum, correct, n = 0.0, 0, 0
    for x_s, y_s, org_s in loader:
        x_s, y_s, org_s = x_s.to(DEVICE), y_s.to(DEVICE), org_s.to(DEVICE)
        x_s, y_s = intraclass_mixup(x_s, y_s, org_s)
        optimizer.zero_grad()
        out      = model(x_s, org_s)
        y_smooth = y_s * (1 - LABEL_SMOOTH) + (1 - y_s) * LABEL_SMOOTH
        loss     = res_criterion(out, y_smooth)
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        loss_sum += loss.item() * len(y_s)
        correct  += ((out > 0) == y_s.bool()).sum().item()
        n        += len(y_s)
    return loss_sum / n, correct / n


@torch.no_grad()
def evaluate_val(model, loader):
    model.eval()
    logits_all, labels_all, orgs_all = [], [], []
    for x, y, org in loader:
        logits_all.append(model(x.to(DEVICE), org.to(DEVICE)).cpu())
        labels_all.append(y); orgs_all.append(org)
    return torch.cat(logits_all), torch.cat(labels_all), torch.cat(orgs_all)


def compute_per_pair_auc(logits, labels, orgs, active_pairs, min_r=5):
    aucs = {}
    for org_id, organism, drug in active_pairs:
        mask = (orgs == org_id)
        if mask.sum() == 0:
            continue
        l = labels[mask].numpy(); p = logits[mask].numpy()
        n_r = int(l.sum())
        if n_r < min_r or n_r == len(l):
            continue
        try:
            aucs[org_id] = float(roc_auc_score(l, p))
        except ValueError:
            pass
    return aucs


def macro_mean_auc(auc_dict):
    vals = [v for v in auc_dict.values() if not math.isnan(v)]
    return float(np.mean(vals)) if vals else float("nan")


def checkpoint_breakdown(model, val_loader, epoch, tr_loss, tr_acc, active_pairs):
    logits, labels, orgs = evaluate_val(model, val_loader)
    auc_dict = compute_per_pair_auc(logits, labels, orgs, active_pairs)
    print(f"\n  ── Epoch {epoch} checkpoint ──")
    for org_id, organism, drug in active_pairs:
        mask = (orgs == org_id)
        n_r  = int(labels[mask].sum()) if mask.sum() > 0 else 0
        auc  = auc_dict.get(org_id, float("nan"))
        note = ("✓" if auc >= 0.80 else "→" if auc >= 0.70 else "✗") if not math.isnan(auc) else "⚠"
        print(f"  [{org_id}] {organism[:22]:22s}/{drug[:14]:14s}  "
              f"AUC={auc:.3f}  R={n_r}  {note}")
    print(f"  loss={tr_loss:.4f}  acc={tr_acc:.3f}  "
          f"macro={macro_mean_auc(auc_dict):.3f}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# §11  CALIBRATION  (source val only)
# ═══════════════════════════════════════════════════════════════════════════════

def learn_ensemble_temperature(models, val_samples):
    """Fit scalar T on ensemble-averaged val logits — fixes single-model bug."""
    assert_source_only("learn_ensemble_temperature")
    loader = make_loader(val_samples, training=False)
    all_logits, all_labels = [], []
    with torch.no_grad():
        for x, y, org in loader:
            seed_lgt = [m(x.to(DEVICE), org.to(DEVICE)).cpu().numpy() for m in models]
            all_logits.append(np.mean(seed_lgt, axis=0))
            all_labels.append(y.numpy())
    logits = np.concatenate(all_logits)
    labels = np.concatenate(all_labels)

    def nll(log_T):
        T     = float(np.exp(log_T[0]))
        probs = np.clip(1.0 / (1.0 + np.exp(-logits / T)), 1e-7, 1 - 1e-7)
        return -np.mean(labels * np.log(probs) + (1 - labels) * np.log(1 - probs))

    if scipy_optimize is None:
        log_T_grid = np.linspace(math.log(0.25), math.log(4.0), 121)
        losses = [nll([log_T]) for log_T in log_T_grid]
        T_opt = float(np.exp(log_T_grid[int(np.argmin(losses))]))
        print("  scipy not available; used grid-search temperature fallback")
    else:
        T_opt = float(np.exp(scipy_optimize.minimize(nll, [0.0], method="Nelder-Mead").x[0]))
    print(f"  Ensemble temperature T = {T_opt:.4f}")
    return T_opt


def _clip_prob(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return 0.5
    return float(np.clip(v if math.isfinite(v) else 0.5, PROB_EPS, 1.0 - PROB_EPS))


def _best_youden_threshold(labels, probs, default=0.5):
    fpr, tpr, ths = roc_curve(labels, probs)
    j = tpr + (1 - fpr) - 1
    finite = np.flatnonzero(np.isfinite(ths))
    if finite.size == 0:
        return _clip_prob(default), float("nan"), float("nan"), float("nan")
    best = int(finite[np.argmax(j[finite])])
    t = _clip_prob(ths[best])
    return t, float(j[best]), float(tpr[best]), float(1 - fpr[best])


def find_thresholds_per_pair(models, val_samples, temperature, active_pairs):
    """Per-pair Youden thresholds on ensemble-averaged val probs."""
    assert_source_only("find_thresholds_per_pair")
    loader = make_loader(val_samples, training=False)
    all_probs, all_labels, all_orgs = [], [], []
    with torch.no_grad():
        for x, y, org in loader:
            seed_p = []
            for m in models:
                m.eval()
                seed_p.append(torch.sigmoid(m(x.to(DEVICE), org.to(DEVICE)).cpu() / temperature).numpy())
            all_probs.append(np.mean(seed_p, axis=0))
            all_labels.append(y.numpy()); all_orgs.append(org.numpy())
    probs  = np.concatenate(all_probs)
    labels = np.concatenate(all_labels)
    orgs   = np.concatenate(all_orgs)

    thresholds = {}
    for org_id, organism, drug in active_pairs:
        mask = (orgs == org_id)
        if mask.sum() == 0:
            thresholds[org_id] = 0.5; continue
        p = probs[mask]; l = labels[mask]
        n_r = int(l.sum())
        if n_r < MIN_R_VAL or n_r == len(l):
            thresholds[org_id] = 0.5
            print(f"  [{org_id}] {organism[:22]:22s}/{drug}: n_r={n_r} → 0.5")
            continue
        t, youden, sens, spec = _best_youden_threshold(l, p)
        thresholds[org_id] = t
        print(f"  [{org_id}] {organism[:22]:22s}/{drug:14s}  "
              f"thresh={t:.4f}  J={youden:.3f}  Sens={sens:.3f}  Spec={spec:.3f}")
    return thresholds


def capped_logit_shift(base_prev, target_prev,
                       max_odds_shift=MAX_PREVALENCE_ODDS_SHIFT):
    base_prev   = _clip_prob(base_prev)
    target_prev = _clip_prob(target_prev)
    shift = (math.log(target_prev / (1.0 - target_prev))
             - math.log(base_prev / (1.0 - base_prev)))
    if max_odds_shift is None:
        return shift
    max_odds_shift = float(max_odds_shift)
    if not math.isfinite(max_odds_shift):
        return shift
    cap = math.log(max(max_odds_shift, 1.0))
    return float(np.clip(shift, -cap, cap))


def _prevalence_shift(base_thresh, base_prev, target_prev,
                      max_odds_shift=MAX_PREVALENCE_ODDS_SHIFT):
    """Log-odds threshold shift. Returns base_thresh if either prev is None."""
    base_thresh = _clip_prob(base_thresh)
    if target_prev is None or base_prev is None:
        return base_thresh
    logit_base = math.log(base_thresh / (1.0 - base_thresh))
    logit_shift = capped_logit_shift(base_prev, target_prev, max_odds_shift)
    return _clip_prob(1.0 / (1.0 + math.exp(-(logit_base - logit_shift))))


def site_thresholds_for(base_thresholds, site, active_pairs, prevalence_shift_mode="capped"):
    """Apply per-pair prevalence shift for a given external site."""
    out = {}
    for org_id, organism, drug in active_pairs:
        base_t      = base_thresholds.get(org_id, 0.5)
        if prevalence_shift_mode == "none":
            out[org_id] = base_t
            continue
        pair_prev   = ALL_SITE_PREVALENCE.get((organism, drug), {})
        train_prev  = pair_prev.get(TRAIN_SITE)
        target_prev = pair_prev.get(site)
        if (train_prev is not None and target_prev is not None
                and abs(target_prev - train_prev) > 0.005):
            shifted = _prevalence_shift(base_t, train_prev, target_prev)
            print(f"  [{org_id}] {organism[:20]:20s}/{drug:14s}  "
                  f"prev {train_prev:.2f}→{target_prev:.2f}  "
                  f"thresh {base_t:.3f}→{shifted:.3f}")
            out[org_id] = shifted
        else:
            out[org_id] = base_t
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# §12  BN ADAPTATION  (Validation Protocol §3.3 fix)
# ═══════════════════════════════════════════════════════════════════════════════

@contextlib.contextmanager
def adapted_batchnorm(models, adapt_samples, use_adapt=True, n_batches=20):
    """
    Forward-pass unlabeled spectra through encode() to update BN stats,
    then restore original stats on exit. Labels never accessed.
    """
    if not use_adapt or not adapt_samples:
        yield; return

    saved = []
    for m in models:
        stats = {}
        for name, mod in m.named_modules():
            if isinstance(mod, nn.BatchNorm1d):
                stats[name] = (mod.running_mean.clone(),
                               mod.running_var.clone(),
                               mod.num_batches_tracked.clone())
        saved.append(stats)

    paths   = [s[0] for s in adapt_samples]
    loader  = DataLoader(MAEDataset(paths), batch_size=BATCH_SIZE, shuffle=True,
                         num_workers=2, pin_memory=True)
    for m in models:
        m.train()
        with torch.no_grad():
            for i, x in enumerate(loader):
                if i >= n_batches:
                    break
                m.encode(x.to(DEVICE))
        m.eval()

    try:
        yield
    finally:
        for m, stats in zip(models, saved):
            for name, mod in m.named_modules():
                if isinstance(mod, nn.BatchNorm1d) and name in stats:
                    mod.running_mean.copy_(stats[name][0])
                    mod.running_var.copy_(stats[name][1])
                    mod.num_batches_tracked.copy_(stats[name][2])


# ═══════════════════════════════════════════════════════════════════════════════
# §13  ENSEMBLE INFERENCE
# ═══════════════════════════════════════════════════════════════════════════════

@torch.no_grad()
def ensemble_predict(models, samples, temperature, n_passes=TTA_PASSES):
    """Returns (probs, labels, orgs) averaged across all seeds × TTA passes."""
    if not samples:
        return np.array([]), np.array([]), np.array([])
    clean_loader = make_loader(samples, training=False)
    labels_np = torch.cat([y for _, y, _ in clean_loader]).numpy()
    orgs_np   = torch.cat([o for _, _, o in clean_loader]).numpy()

    all_model_probs = []
    for m in models:
        m.eval()
        pass_probs = []
        for _ in range(n_passes):
            run_p = []
            for x, _, org in make_tta_loader(samples):
                logits = m(x.to(DEVICE), org.to(DEVICE)).cpu() / temperature
                run_p.append(torch.sigmoid(logits).numpy())
            pass_probs.append(np.concatenate(run_p))
        all_model_probs.append(np.mean(pass_probs, axis=0))

    return np.mean(all_model_probs, axis=0), labels_np, orgs_np


# ═══════════════════════════════════════════════════════════════════════════════
# §14  METRICS
# ═══════════════════════════════════════════════════════════════════════════════

def _safe_auc(labels, probs):
    n_r = int(labels.sum())
    if n_r == 0 or n_r == len(labels):
        return float("nan")
    try:
        return float(roc_auc_score(labels, probs))
    except ValueError:
        return float("nan")


def _safe_aupr(labels, probs):
    if int(labels.sum()) == 0:
        return float("nan")
    try:
        return float(average_precision_score(labels, probs))
    except ValueError:
        return float("nan")


def compute_multi_metrics(probs, labels, orgs, thresholds, site_label, active_pairs):
    """
    thresholds: dict {org_id: float} or single float.
    Returns (rows_list, macro_auc).
    """
    rows, per_org_aucs, pooled_p_, pooled_l_ = [], [], [], []

    for org_id, organism, drug in active_pairs:
        mask = (orgs == org_id)
        if not mask.any():
            continue
        p = probs[mask]; l = labels[mask]
        n = len(l); n_r = int(l.sum()); n_s = n - n_r
        thresh = thresholds[org_id] if isinstance(thresholds, dict) else float(thresholds)
        auc  = _safe_auc(l, p)
        aupr = _safe_aupr(l, p)
        pred = (p >= thresh).astype(float)
        tp   = float(((pred == 1) & (l == 1)).sum())
        tn   = float(((pred == 0) & (l == 0)).sum())
        rows.append(dict(
            site=site_label, organism=organism, drug=drug,
            auc=round(auc, 4), aupr=round(aupr, 4),
            sens=round(tp / n_r if n_r > 0 else float("nan"), 4),
            spec=round(tn / n_s if n_s > 0 else float("nan"), 4),
            n=n, n_r=n_r, threshold=round(thresh, 5),
        ))
        if not np.isnan(auc):
            per_org_aucs.append(auc)
        if n_r > 0:
            pooled_p_.append(p); pooled_l_.append(l)

    if pooled_p_:
        pp = np.concatenate(pooled_p_); pl = np.concatenate(pooled_l_)
        n_tot = len(pl); n_r_tot = int(pl.sum())
        pred  = (pp >= 0.5).astype(float)
        tp    = float(((pred == 1) & (pl == 1)).sum())
        tn    = float(((pred == 0) & (pl == 0)).sum())
        rows.append(dict(
            site=site_label, organism="__micro__", drug="pooled",
            auc=round(_safe_auc(pl, pp), 4), aupr=round(_safe_aupr(pl, pp), 4),
            sens=round(tp / n_r_tot if n_r_tot > 0 else float("nan"), 4),
            spec=round(tn / (n_tot - n_r_tot) if (n_tot - n_r_tot) > 0 else float("nan"), 4),
            n=n_tot, n_r=n_r_tot, threshold=0.5,
        ))

    return rows, float(np.mean(per_org_aucs)) if per_org_aucs else float("nan")


def compute_subset_metrics(probs, labels, orgs, thresholds, site_label, active_pairs,
                           subset_pairs=RUN14_OVERLAP_PAIRS, subset_name="run14_overlap"):
    """Summarize a named subset, mainly for direct run14 overlap comparison."""
    pair_to_id = {(organism, drug): org_id for org_id, organism, drug in active_pairs}
    subset_ids = [pair_to_id[pair] for pair in subset_pairs if pair in pair_to_id]
    if not subset_ids or len(labels) == 0:
        return dict(site=site_label, subset=subset_name, n=0, n_r=0,
                    macro_auc=float("nan"), pooled_auc=float("nan"),
                    mean_sens=float("nan"), mean_spec=float("nan"),
                    included_pairs=0)

    pair_aucs, pair_sens, pair_spec, pooled_p, pooled_l = [], [], [], [], []
    for org_id in subset_ids:
        mask = (orgs == org_id)
        if not mask.any():
            continue
        p = probs[mask]; l = labels[mask]
        n_r = int(l.sum()); n_s = len(l) - n_r
        auc = _safe_auc(l, p)
        if not np.isnan(auc):
            pair_aucs.append(auc)
        thresh = thresholds.get(org_id, 0.5) if isinstance(thresholds, dict) else float(thresholds)
        pred = (p >= thresh).astype(float)
        tp = float(((pred == 1) & (l == 1)).sum())
        tn = float(((pred == 0) & (l == 0)).sum())
        if n_r > 0:
            pair_sens.append(tp / n_r)
        if n_s > 0:
            pair_spec.append(tn / n_s)
        if n_r > 0:
            pooled_p.append(p); pooled_l.append(l)

    if pooled_p:
        pp = np.concatenate(pooled_p)
        pl = np.concatenate(pooled_l)
        n_total = len(pl)
        n_r_total = int(pl.sum())
        pooled_auc = _safe_auc(pl, pp)
    else:
        n_total = 0
        n_r_total = 0
        pooled_auc = float("nan")

    return dict(
        site=site_label, subset=subset_name,
        n=n_total, n_r=n_r_total, included_pairs=len(subset_ids),
        macro_auc=round(float(np.mean(pair_aucs)), 4) if pair_aucs else float("nan"),
        pooled_auc=round(pooled_auc, 4) if not np.isnan(pooled_auc) else float("nan"),
        mean_sens=round(float(np.mean(pair_sens)), 4) if pair_sens else float("nan"),
        mean_spec=round(float(np.mean(pair_spec)), 4) if pair_spec else float("nan"),
    )


def compute_per_organism_summary(metric_rows):
    """Aggregate pair rows into per-site, per-organism mean metrics."""
    rows = [r for r in metric_rows if r.get("organism") != "__micro__"]
    if not rows:
        return []
    df = pd.DataFrame(rows)
    out = []
    for (site, organism), group in df.groupby(["site", "organism"], dropna=False):
        out.append(dict(
            site=site, organism=organism, n_pairs=int(len(group)),
            mean_auc=round(float(group["auc"].mean()), 4),
            mean_aupr=round(float(group["aupr"].mean()), 4),
            mean_sens=round(float(group["sens"].mean()), 4),
            mean_spec=round(float(group["spec"].mean()), 4),
            n=int(group["n"].sum()), n_r=int(group["n_r"].sum()),
        ))
    return out


def compute_threshold_sanity_report(metric_rows, site_threshold_maps=None,
                                    base_thresholds=None, active_pairs=None,
                                    sens_floor=0.10, spec_floor=0.10,
                                    max_threshold_move=0.25):
    """Flag extreme sensitivity/specificity or large prevalence-threshold moves."""
    report = []
    for row in metric_rows:
        if row.get("organism") == "__micro__":
            continue
        sens = float(row.get("sens", float("nan")))
        spec = float(row.get("spec", float("nan")))
        if math.isfinite(sens) and sens < sens_floor:
            report.append(dict(
                site=row["site"], organism=row["organism"], drug=row["drug"],
                issue="low_sensitivity", value=round(sens, 4),
                threshold=row.get("threshold"), n=row.get("n"), n_r=row.get("n_r"),
            ))
        if math.isfinite(spec) and spec < spec_floor:
            report.append(dict(
                site=row["site"], organism=row["organism"], drug=row["drug"],
                issue="low_specificity", value=round(spec, 4),
                threshold=row.get("threshold"), n=row.get("n"), n_r=row.get("n_r"),
            ))

    if site_threshold_maps and base_thresholds and active_pairs:
        for site, thresh_map in site_threshold_maps.items():
            for org_id, organism, drug in active_pairs:
                base_t = base_thresholds.get(org_id, 0.5)
                site_t = thresh_map.get(org_id, base_t)
                move = abs(site_t - base_t)
                if move > max_threshold_move:
                    report.append(dict(
                        site=site, organism=organism, drug=drug,
                        issue="large_threshold_move", value=round(move, 4),
                        threshold=round(site_t, 5), base_threshold=round(base_t, 5),
                    ))
    return report


def _bootstrap_metric_ci(labels, probs, metric_fn, n_boot=STAT_BOOTSTRAP_N,
                         seed=STAT_BOOTSTRAP_SEED):
    labels = np.asarray(labels).astype(int)
    probs = np.asarray(probs, dtype=float)
    if len(labels) == 0 or labels.sum() == 0 or labels.sum() == len(labels):
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(labels), len(labels))
        y = labels[idx]
        if y.sum() == 0 or y.sum() == len(y):
            continue
        vals.append(metric_fn(y, probs[idx]))
    if not vals:
        return float("nan"), float("nan")
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def _resampled_pair_metrics(probs, labels, orgs, active_pairs, rng):
    """Bootstrap each pair within a site and return pair-level AUC/AUPR rows."""
    rows = []
    for org_id, organism, drug in active_pairs:
        mask = (orgs == org_id)
        if not mask.any():
            continue
        p = probs[mask]
        l = labels[mask]
        if len(l) == 0:
            continue
        idx = rng.integers(0, len(l), len(l))
        y = l[idx]
        pp = p[idx]
        auc = _safe_auc(y, pp)
        aupr = _safe_aupr(y, pp)
        rows.append(dict(
            org_id=org_id, organism=organism, drug=drug,
            expected_detectability=_mechanism_metadata(
                organism, drug)["expected_detectability"],
            auc=auc, aupr=aupr,
        ))
    return rows


def _mechanism_contrast_from_pair_rows(pair_rows, metric="auc"):
    more = [
        _finite_float(r.get(metric))
        for r in pair_rows
        if r.get("expected_detectability") == "more_detectable"
    ]
    less = [
        _finite_float(r.get(metric))
        for r in pair_rows
        if r.get("expected_detectability") == "less_detectable"
    ]
    more = [v for v in more if math.isfinite(v)]
    less = [v for v in less if math.isfinite(v)]
    if not more or not less:
        return float("nan")
    return float(np.mean(more) - np.mean(less))


def _mechanism_permutation_pvalue(pair_rows):
    observed = _mechanism_contrast_from_pair_rows(pair_rows, "auc")
    if not math.isfinite(observed):
        return observed, float("nan"), 0
    vals = []
    aucs = [_finite_float(r.get("auc")) for r in pair_rows]
    labels = [r.get("expected_detectability") for r in pair_rows]
    valid = [
        (auc, label) for auc, label in zip(aucs, labels)
        if math.isfinite(auc) and label in ("more_detectable", "less_detectable")
    ]
    if len(valid) < 3:
        return observed, float("nan"), 0
    aucs = [v[0] for v in valid]
    n_more = sum(1 for _, label in valid if label == "more_detectable")
    from itertools import combinations
    for more_idx in combinations(range(len(aucs)), n_more):
        more_idx = set(more_idx)
        more = [auc for i, auc in enumerate(aucs) if i in more_idx]
        less = [auc for i, auc in enumerate(aucs) if i not in more_idx]
        vals.append(float(np.mean(more) - np.mean(less)))
    if not vals:
        return observed, float("nan"), 0
    # One-sided hypothesis: more-detectable mechanisms transfer better.
    permutation_p = (sum(v >= observed - 1e-12 for v in vals) + 1) / (len(vals) + 1)
    return observed, float(permutation_p), len(vals)


def compute_prediction_statistical_report(cnn_cache, active_pairs,
                                          n_boot=STAT_BOOTSTRAP_N):
    """
    Bootstrap 95% CIs for CNN AUC/AUPR, macro AUC/AUPR, and the
    pre-specified mechanism-group contrast. Also run an exact pair-label
    permutation test for the mechanism contrast at each site.
    """
    pair_ci_rows, macro_ci_rows = [], []
    contrast_ci_rows, permutation_rows = [], []

    for site_label, (probs, labels, orgs) in cnn_cache.items():
        site_pair_rows = []
        for org_id, organism, drug in active_pairs:
            mask = (orgs == org_id)
            if not mask.any():
                continue
            p = probs[mask]
            l = labels[mask]
            auc = _safe_auc(l, p)
            aupr = _safe_aupr(l, p)
            auc_lo, auc_hi = _bootstrap_metric_ci(
                l, p, _safe_auc, n_boot=n_boot,
                seed=STAT_BOOTSTRAP_SEED + org_id + 101 * len(site_label))
            aupr_lo, aupr_hi = _bootstrap_metric_ci(
                l, p, _safe_aupr, n_boot=n_boot,
                seed=STAT_BOOTSTRAP_SEED + org_id + 211 * len(site_label))
            meta = _mechanism_metadata(organism, drug)
            row = dict(
                site=site_label, org_id=org_id,
                organism=organism, drug=drug,
                expected_detectability=meta["expected_detectability"],
                mechanism_class=meta["mechanism_class"],
                auc=round(auc, 4) if math.isfinite(auc) else float("nan"),
                auc_bootstrap_ci_low=round(auc_lo, 4) if math.isfinite(auc_lo) else float("nan"),
                auc_bootstrap_ci_high=round(auc_hi, 4) if math.isfinite(auc_hi) else float("nan"),
                aupr=round(aupr, 4) if math.isfinite(aupr) else float("nan"),
                aupr_bootstrap_ci_low=round(aupr_lo, 4) if math.isfinite(aupr_lo) else float("nan"),
                aupr_bootstrap_ci_high=round(aupr_hi, 4) if math.isfinite(aupr_hi) else float("nan"),
                n=int(len(l)), n_r=int(l.sum()),
            )
            pair_ci_rows.append(row)
            site_pair_rows.append(row)

        observed_macro_auc = _mean_finite([r["auc"] for r in site_pair_rows])
        observed_macro_aupr = _mean_finite([r["aupr"] for r in site_pair_rows])
        rng = np.random.default_rng(STAT_BOOTSTRAP_SEED + len(site_label))
        macro_auc_vals, macro_aupr_vals = [], []
        contrast_auc_vals, contrast_aupr_vals = [], []
        for _ in range(n_boot):
            boot_rows = _resampled_pair_metrics(probs, labels, orgs, active_pairs, rng)
            auc_vals = [_finite_float(r["auc"]) for r in boot_rows]
            aupr_vals = [_finite_float(r["aupr"]) for r in boot_rows]
            auc_vals = [v for v in auc_vals if math.isfinite(v)]
            aupr_vals = [v for v in aupr_vals if math.isfinite(v)]
            if auc_vals:
                macro_auc_vals.append(float(np.mean(auc_vals)))
            if aupr_vals:
                macro_aupr_vals.append(float(np.mean(aupr_vals)))
            contrast_auc = _mechanism_contrast_from_pair_rows(boot_rows, "auc")
            contrast_aupr = _mechanism_contrast_from_pair_rows(boot_rows, "aupr")
            if math.isfinite(contrast_auc):
                contrast_auc_vals.append(contrast_auc)
            if math.isfinite(contrast_aupr):
                contrast_aupr_vals.append(contrast_aupr)

        macro_auc_lo, macro_auc_hi = (
            np.percentile(macro_auc_vals, [2.5, 97.5])
            if macro_auc_vals else (float("nan"), float("nan"))
        )
        macro_aupr_lo, macro_aupr_hi = (
            np.percentile(macro_aupr_vals, [2.5, 97.5])
            if macro_aupr_vals else (float("nan"), float("nan"))
        )
        macro_ci_rows.append(dict(
            site=site_label, row_type="macro_bootstrap_ci",
            macro_auc=round(observed_macro_auc, 4) if math.isfinite(observed_macro_auc) else float("nan"),
            macro_auc_ci_low=round(float(macro_auc_lo), 4) if math.isfinite(float(macro_auc_lo)) else float("nan"),
            macro_auc_ci_high=round(float(macro_auc_hi), 4) if math.isfinite(float(macro_auc_hi)) else float("nan"),
            macro_aupr=round(observed_macro_aupr, 4) if math.isfinite(observed_macro_aupr) else float("nan"),
            macro_aupr_ci_low=round(float(macro_aupr_lo), 4) if math.isfinite(float(macro_aupr_lo)) else float("nan"),
            macro_aupr_ci_high=round(float(macro_aupr_hi), 4) if math.isfinite(float(macro_aupr_hi)) else float("nan"),
            n_boot=len(macro_auc_vals),
        ))

        observed_contrast, permutation_p, n_perm = _mechanism_permutation_pvalue(site_pair_rows)
        contrast_auc_lo, contrast_auc_hi = (
            np.percentile(contrast_auc_vals, [2.5, 97.5])
            if contrast_auc_vals else (float("nan"), float("nan"))
        )
        contrast_aupr_lo, contrast_aupr_hi = (
            np.percentile(contrast_aupr_vals, [2.5, 97.5])
            if contrast_aupr_vals else (float("nan"), float("nan"))
        )
        contrast_ci_rows.append(dict(
            site=site_label,
            contrast="more_detectable_minus_less_detectable",
            auc_contrast=round(observed_contrast, 4) if math.isfinite(observed_contrast) else float("nan"),
            auc_contrast_ci_low=round(float(contrast_auc_lo), 4) if math.isfinite(float(contrast_auc_lo)) else float("nan"),
            auc_contrast_ci_high=round(float(contrast_auc_hi), 4) if math.isfinite(float(contrast_auc_hi)) else float("nan"),
            aupr_contrast_ci_low=round(float(contrast_aupr_lo), 4) if math.isfinite(float(contrast_aupr_lo)) else float("nan"),
            aupr_contrast_ci_high=round(float(contrast_aupr_hi), 4) if math.isfinite(float(contrast_aupr_hi)) else float("nan"),
            n_boot=len(contrast_auc_vals),
        ))
        permutation_rows.append(dict(
            site=site_label,
            contrast="more_detectable_minus_less_detectable",
            observed_auc_contrast=round(observed_contrast, 4) if math.isfinite(observed_contrast) else float("nan"),
            permutation_p=round(permutation_p, 4) if math.isfinite(permutation_p) else float("nan"),
            n_permutations=n_perm,
            hypothesis="more_detectable_auc_greater",
        ))

    return dict(
        pair_ci_rows=pair_ci_rows,
        macro_ci_rows=macro_ci_rows,
        contrast_ci_rows=contrast_ci_rows,
        permutation_rows=permutation_rows,
    )


def build_statistical_reporting_markdown(pair_ci_rows, macro_ci_rows,
                                         contrast_ci_rows, permutation_rows):
    lines = [
        "# Statistical Reporting\n\n",
        "Bootstrap intervals are computed from held-out prediction scores and labels. "
        "Mechanism permutation tests shuffle the pre-specified mechanism labels across "
        "available organism-drug pairs within each site. These tests are exploratory; "
        "include a multiple-comparison caution because clinical4 has only two pairs per "
        "mechanism group when all pairs are present.\n\n",
    ]
    if macro_ci_rows:
        lines.append("## Macro AUC/AUPR Bootstrap CIs\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(macro_ci_rows)))
        lines.append("\n")
    if contrast_ci_rows:
        lines.append("## Mechanism-Group Contrast Bootstrap CIs\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(contrast_ci_rows)))
        lines.append("\n")
    if permutation_rows:
        lines.append("## Mechanism Label Permutation Tests\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(permutation_rows)))
        lines.append("\n")
    if pair_ci_rows:
        lines.append("## Pair-Level AUC/AUPR Bootstrap CIs\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(pair_ci_rows), columns=[
            "site", "organism", "drug", "auc", "auc_bootstrap_ci_low",
            "auc_bootstrap_ci_high", "aupr", "aupr_bootstrap_ci_low",
            "aupr_bootstrap_ci_high", "n", "n_r",
        ]))
    return "".join(lines)


def _finite_float(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return value if math.isfinite(value) else float("nan")


def _mean_finite(values):
    vals = [_finite_float(v) for v in values]
    vals = [v for v in vals if math.isfinite(v)]
    return float(np.mean(vals)) if vals else float("nan")


def _site_sort_key(site):
    order = {f"A-{TEST_YEAR}": 0, "DRIAMS-B": 1, "DRIAMS-C": 2, "DRIAMS-D": 3}
    return order.get(site, 99), site


def _format_markdown_value(value):
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.3f}" if math.isfinite(float(value)) else "-"
    if value is None:
        return "-"
    text = str(value)
    return text if text.lower() != "nan" else "-"


def dataframe_to_markdown(df, columns=None):
    """Small dependency-free markdown renderer for paper artifacts."""
    if df.empty:
        return "_No rows._\n"
    if columns is None:
        columns = list(df.columns)
    header = "| " + " | ".join(columns) + " |\n"
    sep = "| " + " | ".join("---" for _ in columns) + " |\n"
    lines = [header, sep]
    for _, row in df.iterrows():
        vals = [_format_markdown_value(row.get(c, "")) for c in columns]
        lines.append("| " + " | ".join(vals) + " |\n")
    return "".join(lines)


def write_markdown_table(df, path, title, note=None, columns=None):
    lines = [f"# {title}\n\n"]
    if note:
        lines.append(note.rstrip() + "\n\n")
    lines.append(dataframe_to_markdown(df, columns=columns))
    with open(path, "w") as f:
        f.writelines(lines)


def _matching_metric(df, site, organism, drug, model=None):
    if df.empty:
        return {}
    mask = (
        (df["site"] == site)
        & (df["organism"] == organism)
        & (df["drug"] == drug)
    )
    if model is not None and "model" in df.columns:
        mask = mask & (df["model"] == model)
    rows = df[mask]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def build_core_results_table(cnn_rows, lgbm_rows, active_pairs):
    """
    Build the paper's core model-comparison table at pair and site level.
    CNN rows use calibrated pair thresholds; LightGBM rows use its 0.5 baseline.
    """
    cnn_df = pd.DataFrame([r for r in cnn_rows if r.get("organism") != "__micro__"])
    if cnn_df.empty:
        return []
    lgbm_df = pd.DataFrame(lgbm_rows) if lgbm_rows else pd.DataFrame()

    pair_order = {(org, drug): i for i, org, drug in active_pairs}
    rows = []
    cnn_df = cnn_df.sort_values(
        by=["site", "organism", "drug"],
        key=lambda col: col.map(_site_sort_key) if col.name == "site" else col,
    )

    for _, cnn in cnn_df.iterrows():
        site = cnn["site"]
        organism = cnn["organism"]
        drug = cnn["drug"]
        single = _matching_metric(lgbm_df, site, organism, drug, "lgbm_single")
        multi = _matching_metric(lgbm_df, site, organism, drug, "lgbm_multi")
        lgbm_single_auc = _finite_float(single.get("auc", float("nan")))
        lgbm_multi_auc = _finite_float(multi.get("auc", float("nan")))
        best_lgbm_vals = [v for v in (lgbm_single_auc, lgbm_multi_auc) if math.isfinite(v)]
        best_lgbm_auc = max(best_lgbm_vals) if best_lgbm_vals else float("nan")
        cnn_auc = _finite_float(cnn.get("auc", float("nan")))
        delta = cnn_auc - best_lgbm_auc if math.isfinite(cnn_auc) and math.isfinite(best_lgbm_auc) else float("nan")

        rows.append(dict(
            row_type="pair",
            site=site,
            pair_index=pair_order.get((organism, drug), -1),
            organism=organism,
            drug=drug,
            n=int(cnn.get("n", 0)),
            n_r=int(cnn.get("n_r", 0)),
            cnn_auc=round(cnn_auc, 4) if math.isfinite(cnn_auc) else float("nan"),
            cnn_aupr=round(_finite_float(cnn.get("aupr", float("nan"))), 4),
            cnn_sens=round(_finite_float(cnn.get("sens", float("nan"))), 4),
            cnn_spec=round(_finite_float(cnn.get("spec", float("nan"))), 4),
            lgbm_single_auc=round(lgbm_single_auc, 4) if math.isfinite(lgbm_single_auc) else float("nan"),
            lgbm_single_sens=round(_finite_float(single.get("sens", float("nan"))), 4),
            lgbm_single_spec=round(_finite_float(single.get("spec", float("nan"))), 4),
            lgbm_multi_auc=round(lgbm_multi_auc, 4) if math.isfinite(lgbm_multi_auc) else float("nan"),
            lgbm_multi_sens=round(_finite_float(multi.get("sens", float("nan"))), 4),
            lgbm_multi_spec=round(_finite_float(multi.get("spec", float("nan"))), 4),
            best_lgbm_auc=round(best_lgbm_auc, 4) if math.isfinite(best_lgbm_auc) else float("nan"),
            cnn_minus_best_lgbm_auc=round(delta, 4) if math.isfinite(delta) else float("nan"),
        ))

    for site, group in cnn_df.groupby("site", dropna=False):
        site_lgbm = lgbm_df[lgbm_df["site"] == site] if not lgbm_df.empty else pd.DataFrame()
        single_group = site_lgbm[site_lgbm["model"] == "lgbm_single"] if not site_lgbm.empty else pd.DataFrame()
        multi_group = site_lgbm[site_lgbm["model"] == "lgbm_multi"] if not site_lgbm.empty else pd.DataFrame()
        cnn_auc = _mean_finite(group["auc"])
        lgbm_single_auc = _mean_finite(single_group["auc"]) if not single_group.empty else float("nan")
        lgbm_multi_auc = _mean_finite(multi_group["auc"]) if not multi_group.empty else float("nan")
        best_lgbm_vals = [v for v in (lgbm_single_auc, lgbm_multi_auc) if math.isfinite(v)]
        best_lgbm_auc = max(best_lgbm_vals) if best_lgbm_vals else float("nan")
        delta = cnn_auc - best_lgbm_auc if math.isfinite(cnn_auc) and math.isfinite(best_lgbm_auc) else float("nan")
        rows.append(dict(
            row_type="macro",
            site=site,
            pair_index=-1,
            organism="Macro mean",
            drug="active pairs",
            n=int(group["n"].sum()),
            n_r=int(group["n_r"].sum()),
            cnn_auc=round(cnn_auc, 4) if math.isfinite(cnn_auc) else float("nan"),
            cnn_aupr=round(_mean_finite(group["aupr"]), 4),
            cnn_sens=round(_mean_finite(group["sens"]), 4),
            cnn_spec=round(_mean_finite(group["spec"]), 4),
            lgbm_single_auc=round(lgbm_single_auc, 4) if math.isfinite(lgbm_single_auc) else float("nan"),
            lgbm_single_sens=round(_mean_finite(single_group["sens"]), 4) if not single_group.empty else float("nan"),
            lgbm_single_spec=round(_mean_finite(single_group["spec"]), 4) if not single_group.empty else float("nan"),
            lgbm_multi_auc=round(lgbm_multi_auc, 4) if math.isfinite(lgbm_multi_auc) else float("nan"),
            lgbm_multi_sens=round(_mean_finite(multi_group["sens"]), 4) if not multi_group.empty else float("nan"),
            lgbm_multi_spec=round(_mean_finite(multi_group["spec"]), 4) if not multi_group.empty else float("nan"),
            best_lgbm_auc=round(best_lgbm_auc, 4) if math.isfinite(best_lgbm_auc) else float("nan"),
            cnn_minus_best_lgbm_auc=round(delta, 4) if math.isfinite(delta) else float("nan"),
        ))

    return sorted(rows, key=lambda r: (_site_sort_key(r["site"]), r["row_type"] != "macro", r["pair_index"]))


def _mechanism_metadata(organism, drug):
    return MECHANISM_DETECTABILITY.get((organism, drug), {
        "expected_detectability": "unknown",
        "mechanism_class": "unknown",
        "mechanism_label": "Unknown / not pre-specified",
        "mechanism_rationale": "No mechanism detectability hypothesis was pre-specified for this pair.",
    })


def compute_mechanism_detectability_summary(metric_rows, active_pairs):
    """
    Group clinical pair performance by the paper's mechanism hypothesis:
    more_detectable chromosomal_structural pairs versus less_detectable
    mobile_heterogeneous pairs.
    """
    active_set = {(organism, drug) for _, organism, drug in active_pairs}
    pair_rows = []
    for row in metric_rows:
        if row.get("organism") == "__micro__":
            continue
        pair = (row.get("organism"), row.get("drug"))
        if pair not in active_set:
            continue
        meta = _mechanism_metadata(*pair)
        enriched = dict(row)
        enriched.update(meta)
        enriched["mechanism_group_label"] = MECHANISM_GROUP_LABELS.get(
            meta["expected_detectability"], MECHANISM_GROUP_LABELS["unknown"])
        pair_rows.append(enriched)

    if not pair_rows:
        return [], []

    df = pd.DataFrame(pair_rows)
    summary_rows = []
    for (site, expected), group in df.groupby(["site", "expected_detectability"], dropna=False):
        summary_rows.append(dict(
            row_type="mechanism_group",
            site=site,
            expected_detectability=expected,
            mechanism_group_label=MECHANISM_GROUP_LABELS.get(
                expected, MECHANISM_GROUP_LABELS["unknown"]),
            n_pairs=int(len(group)),
            n=int(group["n"].sum()),
            n_r=int(group["n_r"].sum()),
            mean_auc=round(_mean_finite(group["auc"]), 4),
            mean_aupr=round(_mean_finite(group["aupr"]), 4),
            mean_sens=round(_mean_finite(group["sens"]), 4),
            mean_spec=round(_mean_finite(group["spec"]), 4),
            pairs="; ".join(f"{r.organism}/{r.drug}" for r in group.itertuples()),
        ))

    for site, group in df.groupby("site", dropna=False):
        more = group[group["expected_detectability"] == "more_detectable"]
        less = group[group["expected_detectability"] == "less_detectable"]
        if more.empty or less.empty:
            continue
        more_auc = _mean_finite(more["auc"])
        less_auc = _mean_finite(less["auc"])
        summary_rows.append(dict(
            row_type="mechanism_contrast",
            site=site,
            expected_detectability="more_minus_less",
            mechanism_group_label="More detectable minus less detectable",
            n_pairs=int(len(more) + len(less)),
            n=int(more["n"].sum() + less["n"].sum()),
            n_r=int(more["n_r"].sum() + less["n_r"].sum()),
            mean_auc=round(more_auc - less_auc, 4)
                     if math.isfinite(more_auc) and math.isfinite(less_auc)
                     else float("nan"),
            mean_aupr=float("nan"),
            mean_sens=float("nan"),
            mean_spec=float("nan"),
            pairs="contrast",
        ))

    return pair_rows, sorted(summary_rows, key=lambda r: (_site_sort_key(r["site"]), r["row_type"], r["expected_detectability"]))


def build_mechanism_framing_markdown(pair_rows, summary_rows, active_pairs=None):
    active_set = (
        {(organism, drug) for _, organism, drug in active_pairs}
        if active_pairs is not None else set(MECHANISM_DETECTABILITY)
    )
    lines = [
        "# Mechanism Framing\n\n",
        "Hypothesis: MALDI-TOF AMR prediction should transfer better when resistance "
        "creates a stable chromosomal or structural phenotype, and worse when resistance "
        "is mobile, plasmid-mediated, efflux-mediated, or otherwise heterogeneous across strains.\n\n",
        "## Pre-specified Groups\n\n",
    ]
    for expected in ("more_detectable", "less_detectable"):
        label = MECHANISM_GROUP_LABELS[expected]
        lines.append(f"### {label}\n\n")
        for (organism, drug), meta in MECHANISM_DETECTABILITY.items():
            if (organism, drug) not in active_set:
                continue
            if meta["expected_detectability"] != expected:
                continue
            lines.append(
                f"- **{organism} / {drug}**: {meta['mechanism_label']}. "
                f"{meta['mechanism_rationale']}\n"
            )
        lines.append("\n")

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        lines.extend([
            "## Group-Level Results\n\n",
            dataframe_to_markdown(summary_df, columns=[
                "row_type", "site", "mechanism_group_label", "n_pairs",
                "mean_auc", "mean_aupr", "mean_sens", "mean_spec",
            ]),
            "\n\n",
        ])

    if pair_rows:
        pair_df = pd.DataFrame(pair_rows)
        lines.extend([
            "## Pair-Level Mechanism Annotation\n\n",
            dataframe_to_markdown(pair_df, columns=[
                "site", "organism", "drug", "expected_detectability",
                "mechanism_class", "auc", "aupr", "sens", "spec",
            ]),
            "\n",
        ])
    return "".join(lines)


def _select_balanced_saliency_samples(samples, org_id,
                                      max_samples=SALIENCY_MAX_SAMPLES_PER_PAIR):
    pair_samples = sorted([s for s in samples if s[2] == org_id], key=lambda s: s[0])
    pos = [s for s in pair_samples if s[1] == 1]
    neg = [s for s in pair_samples if s[1] == 0]
    half = max_samples // 2
    selected = pos[:half] + neg[:max_samples - min(len(pos), half)]
    seen = {s[0] for s in selected}
    for sample in pair_samples:
        if len(selected) >= max_samples:
            break
        if sample[0] not in seen:
            selected.append(sample)
            seen.add(sample[0])
    return selected[:max_samples]


def _ensemble_batch_probs(models, x, org, temperature):
    model_probs = []
    for model in models:
        model.eval()
        logits = model(x, org) / temperature
        model_probs.append(torch.sigmoid(logits).detach().cpu().numpy())
    return np.mean(model_probs, axis=0)


def compute_occlusion_saliency(models, samples, temperature, active_pairs,
                               max_samples_per_pair=SALIENCY_MAX_SAMPLES_PER_PAIR,
                               top_k=SALIENCY_TOP_K,
                               site_label="A-val"):
    """
    Lightweight interpretability: zero one m/z window at a time and rank windows by
    mean absolute change in ensemble probability on source validation spectra.
    """
    rows = []
    if not samples:
        return rows

    for org_id, organism, drug in active_pairs:
        selected = _select_balanced_saliency_samples(
            samples, org_id, max_samples=max_samples_per_pair)
        if not selected:
            continue

        mz_axis = load_mz_axis(selected[0][0])
        xs, ys = [], []
        for path, label, _ in selected:
            try:
                xs.append(load_spectrum(path))
                ys.append(int(label))
            except Exception as exc:
                print(f"  Saliency skip {path}: {exc}")
        if not xs:
            continue

        x = torch.from_numpy(np.stack(xs).astype(np.float32)).unsqueeze(1).to(DEVICE)
        org = torch.full((x.size(0),), org_id, dtype=torch.long, device=DEVICE)

        scores = []
        with torch.no_grad():
            base_probs = _ensemble_batch_probs(models, x, org, temperature)
            for start in range(0, N_BINS, SALIENCY_STRIDE):
                end = min(start + SALIENCY_WINDOW, N_BINS)
                if start >= end:
                    continue
                x_occ = x.clone()
                x_occ[:, :, start:end] = 0.0
                occ_probs = _ensemble_batch_probs(models, x_occ, org, temperature)
                delta = base_probs - occ_probs
                scores.append(dict(
                    bin_start=start,
                    bin_end=end - 1,
                    importance=float(np.mean(np.abs(delta))),
                    signed_delta=float(np.mean(delta)),
                ))

        scores = sorted(scores, key=lambda r: r["importance"], reverse=True)[:top_k]
        for rank, score in enumerate(scores, start=1):
            mz_start_idx = min(score["bin_start"], len(mz_axis) - 1)
            mz_end_idx = min(score["bin_end"], len(mz_axis) - 1)
            rows.append(dict(
                site=site_label,
                org_id=org_id,
                organism=organism,
                drug=drug,
                rank=rank,
                bin_start=score["bin_start"],
                bin_end=score["bin_end"],
                mz_start=round(float(mz_axis[mz_start_idx]), 4),
                mz_end=round(float(mz_axis[mz_end_idx]), 4),
                importance=round(score["importance"], 6),
                signed_delta=round(score["signed_delta"], 6),
                n_samples=len(xs),
                n_r=int(sum(ys)),
            ))
    return rows


def _top_bin_set(group, top_k=SALIENCY_STABILITY_TOP_K):
    if group is None or len(group) == 0:
        return set()
    ordered = group.sort_values(["rank", "importance"], ascending=[True, False])
    return set(
        (int(r.bin_start), int(r.bin_end))
        for r in ordered.head(top_k).itertuples()
    )


def compute_saliency_stability_summary(saliency_rows, active_pairs,
                                       reference_site=None,
                                       top_k=SALIENCY_STABILITY_TOP_K):
    """Compare top occlusion windows across sites with Jaccard overlap."""
    if not saliency_rows:
        return []
    df = pd.DataFrame(saliency_rows)
    if df.empty or "site" not in df.columns:
        return []
    if reference_site is None:
        reference_site = f"A-{TEST_YEAR}" if (df["site"] == f"A-{TEST_YEAR}").any() else "A-val"

    rows = []
    active_lookup = {(org, drug): oid for oid, org, drug in active_pairs}
    for (organism, drug), org_id in active_lookup.items():
        pair_df = df[(df["organism"] == organism) & (df["drug"] == drug)]
        if pair_df.empty:
            continue
        ref_bins = _top_bin_set(pair_df[pair_df["site"] == reference_site], top_k=top_k)
        meta = _mechanism_metadata(organism, drug)
        jaccards = []
        for site in sorted(pair_df["site"].dropna().unique(), key=_site_sort_key):
            if site == reference_site:
                continue
            site_bins = _top_bin_set(pair_df[pair_df["site"] == site], top_k=top_k)
            union = ref_bins | site_bins
            jaccard = len(ref_bins & site_bins) / len(union) if union else float("nan")
            if math.isfinite(jaccard):
                jaccards.append(jaccard)
            rows.append(dict(
                row_type="pair_site",
                reference_site=reference_site,
                site=site,
                org_id=org_id,
                organism=organism,
                drug=drug,
                expected_detectability=meta["expected_detectability"],
                mechanism_class=meta["mechanism_class"],
                jaccard_top_bins=round(jaccard, 4) if math.isfinite(jaccard) else float("nan"),
                n_ref_bins=len(ref_bins),
                n_site_bins=len(site_bins),
            ))
        if jaccards:
            rows.append(dict(
                row_type="pair_mean",
                reference_site=reference_site,
                site="all_external",
                org_id=org_id,
                organism=organism,
                drug=drug,
                expected_detectability=meta["expected_detectability"],
                mechanism_class=meta["mechanism_class"],
                jaccard_top_bins=round(float(np.mean(jaccards)), 4),
                mean_jaccard=round(float(np.mean(jaccards)), 4),
                n_ref_bins=len(ref_bins),
                n_site_bins=top_k,
            ))

    pair_mean = [r for r in rows if r["row_type"] == "pair_mean"]
    if pair_mean:
        mean_df = pd.DataFrame(pair_mean)
        for expected, group in mean_df.groupby("expected_detectability", dropna=False):
            mean_val = _mean_finite(group["mean_jaccard"])
            rows.append(dict(
                row_type="mechanism_group",
                reference_site=reference_site,
                site="all_external",
                org_id=-1,
                organism=MECHANISM_GROUP_LABELS.get(expected, expected),
                drug="saliency stability",
                expected_detectability=expected,
                mechanism_class="group_summary",
                jaccard_top_bins=round(mean_val, 4) if math.isfinite(mean_val) else float("nan"),
                mean_jaccard=round(mean_val, 4) if math.isfinite(mean_val) else float("nan"),
                n_ref_bins=top_k,
                n_site_bins=top_k,
            ))
    return rows


def build_saliency_stability_markdown(stability_rows):
    lines = [
        "# Cross-Site Saliency Stability\n\n",
        "Top m/z windows are compared by Jaccard overlap between the reference A-site "
        "saliency profile and each external site. Higher overlap supports a more stable "
        "spectral mechanism across sites.\n\n",
    ]
    df = pd.DataFrame(stability_rows)
    if df.empty:
        lines.append("_No saliency stability rows were generated._\n")
        return "".join(lines)
    lines.append(dataframe_to_markdown(df, columns=[
        "row_type", "reference_site", "site", "organism", "drug",
        "expected_detectability", "jaccard_top_bins", "mean_jaccard",
    ]))
    lines.append(
        "\n\nInterpretation target: Ec/Cipro and Sa/Oxa should show higher cross-site "
        "top-window overlap than Ec/Amox-Clav and S.epi/Ery if the chromosomal/"
        "structural phenotype hypothesis is correct.\n"
    )
    return "".join(lines)


def _rank_correlation(ref_group, site_group, top_k=SALIENCY_TOP_K):
    if ref_group is None or site_group is None or len(ref_group) == 0 or len(site_group) == 0:
        return float("nan")
    ref_rank = {
        (int(r.bin_start), int(r.bin_end)): int(r.rank)
        for r in ref_group.sort_values("rank").head(top_k).itertuples()
    }
    site_rank = {
        (int(r.bin_start), int(r.bin_end)): int(r.rank)
        for r in site_group.sort_values("rank").head(top_k).itertuples()
    }
    keys = sorted(set(ref_rank) | set(site_rank))
    if len(keys) < 2:
        return float("nan")
    missing = top_k + 1
    x = np.array([ref_rank.get(k, missing) for k in keys], dtype=float)
    y = np.array([site_rank.get(k, missing) for k in keys], dtype=float)
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _broad_bin_set(group, top_k, broad_window=SALIENCY_BROAD_WINDOW):
    bins = _top_bin_set(group, top_k=top_k)
    return {start // broad_window for start, _ in bins}


def _jaccard(a, b):
    union = set(a) | set(b)
    return len(set(a) & set(b)) / len(union) if union else float("nan")


def compute_saliency_robustness_summary(saliency_rows, metric_rows, active_pairs,
                                        reference_site=None,
                                        top_ks=SALIENCY_STABILITY_TOP_KS,
                                        exclude_low_bins_below=SALIENCY_EXCLUDE_LOW_BINS_BELOW):
    """
    Robustness view of saliency stability: top-5 and top-10 Jaccard, rank
    correlation over top windows, broad-window overlap, and a random-window
    null control. This is a saliency-overlap negative control; a fully random
    label retrain remains the stronger follow-up control.
    """
    if not saliency_rows:
        return [], [], []
    df = pd.DataFrame(saliency_rows)
    if df.empty:
        return [], [], []
    if reference_site is None:
        reference_site = f"A-{TEST_YEAR}" if (df["site"] == f"A-{TEST_YEAR}").any() else "A-val"

    metric_df = pd.DataFrame([r for r in metric_rows if r.get("organism") != "__micro__"])
    rng = np.random.default_rng(STAT_BOOTSTRAP_SEED)
    robustness_rows, null_rows = [], []
    active_lookup = {(org, drug): oid for oid, org, drug in active_pairs}

    for (organism, drug), org_id in active_lookup.items():
        pair_df = df[(df["organism"] == organism) & (df["drug"] == drug)]
        if pair_df.empty:
            continue
        meta = _mechanism_metadata(organism, drug)
        filtered_pair_df = pair_df[pair_df["bin_start"].astype(int) >= exclude_low_bins_below]
        ref = filtered_pair_df[filtered_pair_df["site"] == reference_site]
        if ref.empty:
            ref = pair_df[pair_df["site"] == reference_site]
        all_bins = sorted({
            (int(r.bin_start), int(r.bin_end))
            for r in filtered_pair_df.itertuples()
        })
        for site in sorted(pair_df["site"].dropna().unique(), key=_site_sort_key):
            if site == reference_site:
                continue
            site_group = filtered_pair_df[filtered_pair_df["site"] == site]
            if site_group.empty:
                site_group = pair_df[pair_df["site"] == site]
            row = dict(
                row_type="pair_site",
                reference_site=reference_site,
                site=site,
                org_id=org_id,
                organism=organism,
                drug=drug,
                expected_detectability=meta["expected_detectability"],
                mechanism_class=meta["mechanism_class"],
                low_bins_excluded_below=exclude_low_bins_below,
            )
            for top_k in top_ks:
                ref_bins = _top_bin_set(ref, top_k=top_k)
                site_bins = _top_bin_set(site_group, top_k=top_k)
                row[f"top{top_k}_jaccard"] = round(_jaccard(ref_bins, site_bins), 4)
            row["rank_correlation"] = round(
                _rank_correlation(ref, site_group, top_k=max(top_ks)), 4)
            row["broad_window_jaccard"] = round(
                _jaccard(
                    _broad_bin_set(ref, top_k=max(top_ks)),
                    _broad_bin_set(site_group, top_k=max(top_ks)),
                ), 4)
            robustness_rows.append(row)

            if all_bins:
                observed = row.get(f"top{min(top_ks)}_jaccard", float("nan"))
                null_vals = []
                draw_n = min(min(top_ks), len(all_bins))
                for _ in range(SALIENCY_NULL_N):
                    a_idx = rng.choice(len(all_bins), draw_n, replace=False)
                    b_idx = rng.choice(len(all_bins), draw_n, replace=False)
                    a = {all_bins[i] for i in a_idx}
                    b = {all_bins[i] for i in b_idx}
                    null_vals.append(_jaccard(a, b))
                null_mean = float(np.mean(null_vals)) if null_vals else float("nan")
                null_p = (
                    (sum(v >= observed for v in null_vals) + 1) / (len(null_vals) + 1)
                    if null_vals and math.isfinite(observed) else float("nan")
                )
                null_rows.append(dict(
                    site=site,
                    organism=organism,
                    drug=drug,
                    expected_detectability=meta["expected_detectability"],
                    null_control="random_window_control",
                    observed_top5_jaccard=observed,
                    null_mean_jaccard=round(null_mean, 4) if math.isfinite(null_mean) else float("nan"),
                    null_p=round(null_p, 4) if math.isfinite(null_p) else float("nan"),
                    n_random=SALIENCY_NULL_N,
                ))

    alignment_rows = []
    robust_df = pd.DataFrame(robustness_rows)
    if not robust_df.empty and not metric_df.empty:
        for _, organism, drug in active_pairs:
            rpair = robust_df[(robust_df["organism"] == organism) & (robust_df["drug"] == drug)]
            mpair = metric_df[(metric_df["organism"] == organism) & (metric_df["drug"] == drug)]
            if rpair.empty or mpair.empty:
                continue
            external_auc = _mean_finite(
                mpair[mpair["site"].isin(TEST_SITES)]["auc"]
                if "site" in mpair.columns else []
            )
            meta = _mechanism_metadata(organism, drug)
            alignment_rows.append(dict(
                organism=organism,
                drug=drug,
                expected_detectability=meta["expected_detectability"],
                mean_external_auc=round(external_auc, 4) if math.isfinite(external_auc) else float("nan"),
                mean_top5_jaccard=round(_mean_finite(rpair["top5_jaccard"]), 4),
                mean_top10_jaccard=round(_mean_finite(rpair["top10_jaccard"]), 4)
                                   if "top10_jaccard" in rpair.columns else float("nan"),
                mean_rank_correlation=round(_mean_finite(rpair["rank_correlation"]), 4),
                mean_broad_window_jaccard=round(_mean_finite(rpair["broad_window_jaccard"]), 4),
            ))
    return robustness_rows, alignment_rows, null_rows


def build_saliency_robustness_markdown(robustness_rows, alignment_rows, null_rows):
    lines = [
        "# Saliency Robustness Checks\n\n",
        "This table repeats the cross-site saliency comparison with top-5, top-10, "
        "rank-correlation, broader m/z bands, and low-bin sensitivity analyses. "
        "Report both including-low-bin and excluding-low-bin views side by side, "
        "because low m/z windows may be biologically meaningful or preprocessing-sensitive. "
        "The null table is a random-window overlap control; a random-label retrain "
        "is still the stronger negative-control experiment.\n\n",
    ]
    if robustness_rows:
        lines.append("## Cross-Site Robustness Metrics\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(robustness_rows)))
        lines.append("\n")
    if alignment_rows:
        lines.append("## Saliency Stability Versus External AUC\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(alignment_rows)))
        lines.append("\n")
    if null_rows:
        lines.append("## Random-Window Negative Control\n\n")
        lines.append(dataframe_to_markdown(pd.DataFrame(null_rows)))
    return "".join(lines)


def compute_d_site_shift_analysis(metric_rows, active_pairs,
                                  baseline_site=None,
                                  d_site="DRIAMS-D"):
    """Make DRIAMS-D a named community/inpatient distribution-shift finding."""
    baseline_site = baseline_site or f"A-{TEST_YEAR}"
    df = pd.DataFrame([r for r in metric_rows if r.get("organism") != "__micro__"])
    if df.empty:
        return []
    rows = []
    for _, organism, drug in active_pairs:
        base = _matching_metric(df, baseline_site, organism, drug)
        drow = _matching_metric(df, d_site, organism, drug)
        if not base or not drow:
            continue
        base_auc = _finite_float(base.get("auc", float("nan")))
        d_auc = _finite_float(drow.get("auc", float("nan")))
        base_aupr = _finite_float(base.get("aupr", float("nan")))
        d_aupr = _finite_float(drow.get("aupr", float("nan")))
        meta = _mechanism_metadata(organism, drug)
        context = SITE_CONTEXT.get(d_site, {})
        rows.append(dict(
            row_type="pair",
            baseline_site=baseline_site,
            site=d_site,
            site_context=context.get("site_context", "unknown"),
            organism=organism,
            drug=drug,
            expected_detectability=meta["expected_detectability"],
            mechanism_class=meta["mechanism_class"],
            baseline_auc=round(base_auc, 4) if math.isfinite(base_auc) else float("nan"),
            d_site_auc=round(d_auc, 4) if math.isfinite(d_auc) else float("nan"),
            delta_auc=round(d_auc - base_auc, 4)
                      if math.isfinite(base_auc) and math.isfinite(d_auc)
                      else float("nan"),
            baseline_aupr=round(base_aupr, 4) if math.isfinite(base_aupr) else float("nan"),
            d_site_aupr=round(d_aupr, 4) if math.isfinite(d_aupr) else float("nan"),
            delta_aupr=round(d_aupr - base_aupr, 4)
                       if math.isfinite(base_aupr) and math.isfinite(d_aupr)
                       else float("nan"),
            interpretation=context.get("interpretation", ""),
        ))

    if rows:
        pair_df = pd.DataFrame(rows)
        rows.append(dict(
            row_type="summary",
            baseline_site=baseline_site,
            site=d_site,
            site_context=SITE_CONTEXT.get(d_site, {}).get("site_context", "unknown"),
            organism="All active DRIAMS-D pairs",
            drug="community_acquired_shift_summary",
            expected_detectability="all",
            mechanism_class="site_context",
            baseline_auc=round(_mean_finite(pair_df["baseline_auc"]), 4),
            d_site_auc=round(_mean_finite(pair_df["d_site_auc"]), 4),
            delta_auc=round(_mean_finite(pair_df["delta_auc"]), 4),
            baseline_aupr=round(_mean_finite(pair_df["baseline_aupr"]), 4),
            d_site_aupr=round(_mean_finite(pair_df["d_site_aupr"]), 4),
            delta_aupr=round(_mean_finite(pair_df["delta_aupr"]), 4),
            interpretation=SITE_CONTEXT.get(d_site, {}).get("interpretation", ""),
        ))
    return rows


def build_d_site_shift_markdown(d_rows):
    lines = [
        "# DRIAMS-D Community/External Shift Analysis\n\n",
        "DRIAMS-D is reported as a named site-level finding because its drop is "
        "consistent with a community_acquired distribution shift. This section "
        "separates that biological/site-context interpretation from generic model failure.\n\n",
    ]
    df = pd.DataFrame(d_rows)
    if df.empty:
        lines.append("_No DRIAMS-D rows were generated._\n")
        return "".join(lines)
    lines.append(dataframe_to_markdown(df, columns=[
        "row_type", "baseline_site", "site", "site_context", "organism", "drug",
        "expected_detectability", "baseline_auc", "d_site_auc", "delta_auc",
        "baseline_aupr", "d_site_aupr", "delta_aupr",
    ]))
    lines.append(
        "\n\nSuggested wording: The largest degradation occurred at DRIAMS-D, "
        "consistent with a community-acquired distribution shift rather than an "
        "undifferentiated failure of the MALDI signal.\n"
    )
    return "".join(lines)


# ═══════════════════════════════════════════════════════════════════════════════
# §15  LIGHTGBM BASELINES  (Task 5)
# ═══════════════════════════════════════════════════════════════════════════════

def _load_spectra_matrix(samples):
    X_, y_, org_ = [], [], []
    for path, label, org_id in samples:
        try:
            X_.append(load_spectrum(path))
            y_.append(float(label)); org_.append(int(org_id))
        except Exception as e:
            print(f"  WARNING: {path}: {e}")
    return np.stack(X_), np.array(y_, np.float32), np.array(org_, np.int32)


def train_lgbm_singletask(train_s, val_s, active_pairs):
    try:
        import lightgbm as lgb
    except ImportError:
        print("  lightgbm not installed — skipping"); return {}
    X_tr, y_tr, org_tr = _load_spectra_matrix(train_s)
    X_v,  y_v,  org_v  = _load_spectra_matrix(val_s)
    params = dict(objective="binary", metric="auc", learning_rate=0.05,
                  num_leaves=63, min_child_samples=20, subsample=0.8,
                  colsample_bytree=0.1, reg_lambda=1.0, random_state=42,
                  verbose=-1, n_jobs=-1)
    models = {}
    for org_id, organism, drug in active_pairs:
        tr_m = (org_tr == org_id); v_m = (org_v == org_id)
        if not tr_m.any():
            continue
        cbs   = [lgb.early_stopping(50, verbose=False), lgb.log_evaluation(-1)]
        model = lgb.train(params, lgb.Dataset(X_tr[tr_m], label=y_tr[tr_m]),
                          num_boost_round=500,
                          valid_sets=[lgb.Dataset(X_v[v_m], label=y_v[v_m])],
                          callbacks=cbs)
        models[org_id] = model
        print(f"  [{org_id}] {organism[:22]:22s}/{drug}: "
              f"n={int(tr_m.sum())}  R={int(y_tr[tr_m].sum())}  "
              f"iter={getattr(model,'best_iteration',500)}")
    return models


def train_lgbm_multitask(train_s, val_s):
    try:
        import lightgbm as lgb
    except ImportError:
        return None
    X_tr, y_tr, org_tr = _load_spectra_matrix(train_s)
    X_v,  y_v,  org_v  = _load_spectra_matrix(val_s)
    X_tr_a = np.hstack([X_tr, org_tr.reshape(-1, 1).astype(np.float32)])
    X_v_a  = np.hstack([X_v,  org_v.reshape(-1, 1).astype(np.float32)])
    params = dict(objective="binary", metric="auc", learning_rate=0.05,
                  num_leaves=63, min_child_samples=20, subsample=0.8,
                  colsample_bytree=0.1, reg_lambda=1.0, random_state=42,
                  verbose=-1, n_jobs=-1)
    cbs   = [lgb.early_stopping(50, verbose=False), lgb.log_evaluation(-1)]
    model = lgb.train(params, lgb.Dataset(X_tr_a, label=y_tr),
                      num_boost_round=500,
                      valid_sets=[lgb.Dataset(X_v_a, label=y_v)],
                      callbacks=cbs)
    print(f"  Multi-task LGBM: n={len(y_tr)}  iter={getattr(model,'best_iteration',500)}")
    return model


def evaluate_lgbm_site(lgbm_single, lgbm_multi, site_samples, active_pairs, threshold=0.5):
    X, y, orgs = _load_spectra_matrix(site_samples)
    rows = []
    for variant in ("lgbm_single", "lgbm_multi"):
        if variant == "lgbm_single" and not lgbm_single:
            continue
        if variant == "lgbm_multi" and lgbm_multi is None:
            continue
        X_aug = np.hstack([X, orgs.reshape(-1, 1).astype(np.float32)]) if variant == "lgbm_multi" else None
        for org_id, organism, drug in active_pairs:
            mask = (orgs == org_id)
            if not mask.any():
                continue
            l = y[mask]; n = len(l); n_r = int(l.sum()); n_s = n - n_r
            if variant == "lgbm_single":
                if org_id not in lgbm_single:
                    continue
                p = lgbm_single[org_id].predict(X[mask])
            else:
                p = lgbm_multi.predict(X_aug[mask])
            auc  = _safe_auc(l, p); aupr = _safe_aupr(l, p)
            pred = (p >= threshold).astype(float)
            tp   = float(((pred == 1) & (l == 1)).sum())
            tn   = float(((pred == 0) & (l == 0)).sum())
            rows.append(dict(
                model=variant, organism=organism, drug=drug,
                auc=round(auc, 4), aupr=round(aupr, 4),
                sens=round(tp / n_r if n_r > 0 else float("nan"), 4),
                spec=round(tn / n_s if n_s > 0 else float("nan"), 4),
                n=n, n_r=n_r,
            ))
    return rows


def make_random_cv_split(all_samples, holdout=RANDOM_CV_HOLDOUT,
                         random_state=RANDOM_CV_SEED):
    """Random pair-label split across all source years, used only as an inflation diagnostic."""
    source_samples = [(p, l, o) for p, l, o, _ in all_samples]
    return split_train_val_by_pair_label(
        source_samples, val_holdout=holdout, random_state=random_state)


def run_random_cv_inflation_analysis(all_samples, temporal_lgbm_rows, active_pairs):
    """
    Train LightGBM on a random source split and compare against temporal A-2018.
    This is a diagnostic for literature-style random-CV inflation, not the main protocol.
    """
    if not temporal_lgbm_rows:
        return []
    rand_train, rand_holdout = make_random_cv_split(all_samples)
    random_single = train_lgbm_singletask(rand_train, rand_holdout, active_pairs)
    random_multi = train_lgbm_multitask(rand_train, rand_holdout)
    random_rows = evaluate_lgbm_site(
        random_single, random_multi, rand_holdout, active_pairs)
    for row in random_rows:
        row["site"] = "DRIAMS-A-random-holdout"

    temporal_df = pd.DataFrame(temporal_lgbm_rows)
    out = []
    for row in random_rows:
        temporal = _matching_metric(
            temporal_df, f"A-{TEST_YEAR}", row["organism"], row["drug"], row["model"])
        if not temporal:
            continue
        random_auc = _finite_float(row.get("auc", float("nan")))
        temporal_auc = _finite_float(temporal.get("auc", float("nan")))
        random_aupr = _finite_float(row.get("aupr", float("nan")))
        temporal_aupr = _finite_float(temporal.get("aupr", float("nan")))
        out.append(dict(
            row_type="pair",
            model=row["model"],
            organism=row["organism"],
            drug=row["drug"],
            random_site="DRIAMS-A-random-holdout",
            temporal_site=f"A-{TEST_YEAR}",
            random_auc=round(random_auc, 4) if math.isfinite(random_auc) else float("nan"),
            temporal_auc=round(temporal_auc, 4) if math.isfinite(temporal_auc) else float("nan"),
            auc_inflation=round(random_auc - temporal_auc, 4)
                          if math.isfinite(random_auc) and math.isfinite(temporal_auc)
                          else float("nan"),
            random_aupr=round(random_aupr, 4) if math.isfinite(random_aupr) else float("nan"),
            temporal_aupr=round(temporal_aupr, 4) if math.isfinite(temporal_aupr) else float("nan"),
            aupr_inflation=round(random_aupr - temporal_aupr, 4)
                           if math.isfinite(random_aupr) and math.isfinite(temporal_aupr)
                           else float("nan"),
            n_random=int(row.get("n", 0)),
            n_r_random=int(row.get("n_r", 0)),
        ))

    if out:
        out_df = pd.DataFrame(out)
        for model, group in out_df.groupby("model", dropna=False):
            out.append(dict(
                row_type="macro",
                model=model,
                organism="Macro mean",
                drug="active pairs",
                random_site="DRIAMS-A-random-holdout",
                temporal_site=f"A-{TEST_YEAR}",
                random_auc=round(_mean_finite(group["random_auc"]), 4),
                temporal_auc=round(_mean_finite(group["temporal_auc"]), 4),
                auc_inflation=round(_mean_finite(group["auc_inflation"]), 4),
                random_aupr=round(_mean_finite(group["random_aupr"]), 4),
                temporal_aupr=round(_mean_finite(group["temporal_aupr"]), 4),
                aupr_inflation=round(_mean_finite(group["aupr_inflation"]), 4),
                n_random=int(group["n_random"].sum()),
                n_r_random=int(group["n_r_random"].sum()),
            ))
    return out


def build_temporal_vs_random_cv_markdown(rows):
    lines = [
        "# Temporal Holdout vs Random CV Diagnostic\n\n",
        "This LightGBM diagnostic estimates how much a literature-style random split "
        "can inflate apparent source-site performance relative to the stricter "
        f"temporal holdout on A-{TEST_YEAR}. It is not used for calibration or model selection.\n\n",
    ]
    df = pd.DataFrame(rows)
    if df.empty:
        lines.append("_Random-CV diagnostic was skipped or produced no rows._\n")
        return "".join(lines)
    lines.append(dataframe_to_markdown(df, columns=[
        "row_type", "model", "organism", "drug",
        "random_auc", "temporal_auc", "auc_inflation",
        "random_aupr", "temporal_aupr", "aupr_inflation",
    ]))
    lines.append(
        "\n\nPaper framing: random cross-validation should be presented as an optimistic "
        "upper-bound diagnostic; temporal and external-site evaluation remain the main evidence.\n"
    )
    return "".join(lines)


def build_evaluation_critique_table(cnn_macro_rows, lgbm_rows,
                                    random_cv_rows=None):
    """
    One table for the evaluation critique: random source split vs temporal
    A-2018 vs external B/C/D. This keeps the paper wording at "random CV is
    optimistic" rather than overclaiming that all prior work is invalid.
    """
    rows = []
    random_cv_rows = random_cv_rows or []
    for r in cnn_macro_rows:
        site = r.get("site")
        rows.append(dict(
            model="CNN",
            evaluation=(
                "temporal_A-2018" if site == f"A-{TEST_YEAR}"
                else "external_site"
            ),
            site=site,
            macro_auc=round(_finite_float(r.get("macro_auc")), 4),
            random_auc=float("nan"),
            temporal_auc=float("nan"),
            auc_inflation=float("nan"),
            note="Primary evaluation uses temporal and external-site splits.",
        ))

    lgbm_df = pd.DataFrame(lgbm_rows) if lgbm_rows else pd.DataFrame()
    if not lgbm_df.empty:
        for (model, site), group in lgbm_df.groupby(["model", "site"], dropna=False):
            rows.append(dict(
                model=model,
                evaluation=(
                    "temporal_A-2018" if site == f"A-{TEST_YEAR}"
                    else "external_site"
                ),
                site=site,
                macro_auc=round(_mean_finite(group["auc"]), 4),
                random_auc=float("nan"),
                temporal_auc=float("nan"),
                auc_inflation=float("nan"),
                note="Classical baseline evaluated with the same temporal/external protocol.",
            ))

    for r in random_cv_rows:
        if r.get("row_type") != "macro":
            continue
        rows.append(dict(
            model=r.get("model"),
            evaluation="random_cv_inflation_diagnostic",
            site="DRIAMS-A-random-holdout_vs_A-2018",
            macro_auc=float("nan"),
            random_auc=round(_finite_float(r.get("random_auc")), 4),
            temporal_auc=round(_finite_float(r.get("temporal_auc")), 4),
            auc_inflation=round(_finite_float(r.get("auc_inflation")), 4),
            note="Random CV is optimistic; temporal/external results are the main evidence.",
        ))
    return rows


def build_evaluation_critique_markdown(rows):
    lines = [
        "# Evaluation Critique Table\n\n",
        "This table places random-CV diagnostics next to temporal and external-site "
        "evaluation. The intended claim is conservative: random CV is optimistic, "
        "while temporal and external sites are the deployment-facing evidence.\n\n",
    ]
    lines.append(dataframe_to_markdown(pd.DataFrame(rows)))
    return "".join(lines)


def compute_mechanism_confound_checks(metric_rows, active_pairs):
    """Small reviewer-facing checks for sample size, prevalence, organism, site, and pair effects."""
    rows = []
    df = pd.DataFrame([r for r in metric_rows if r.get("organism") != "__micro__"])
    if df.empty:
        return rows
    df["prevalence"] = df["n_r"].astype(float) / df["n"].astype(float)
    df["sample_size"] = df["n"].astype(float)
    df["expected_detectability"] = [
        _mechanism_metadata(r.organism, r.drug)["expected_detectability"]
        for r in df.itertuples()
    ]
    df["mechanism_class"] = [
        _mechanism_metadata(r.organism, r.drug)["mechanism_class"]
        for r in df.itertuples()
    ]

    for site, group in df.groupby("site", dropna=False):
        for col, label in [("sample_size", "sample_size"), ("prevalence", "prevalence")]:
            vals = group[col].astype(float).values
            aucs = group["auc"].astype(float).values
            if len(vals) >= 2 and np.std(vals) > 0 and np.std(aucs) > 0:
                corr = float(np.corrcoef(vals, aucs)[0, 1])
            else:
                corr = float("nan")
            rows.append(dict(
                check_type=f"{label}_auc_correlation",
                site=site,
                excluded="none",
                grouping=label,
                value=round(corr, 4) if math.isfinite(corr) else float("nan"),
                n_pairs=int(len(group)),
                interpretation=f"Checks whether mechanism contrast could be explained by {label}.",
            ))

        for organism, org_group in group.groupby("organism", dropna=False):
            rows.append(dict(
                check_type="organism_identity",
                site=site,
                excluded="none",
                grouping=organism,
                value=round(_mean_finite(org_group["auc"]), 4),
                n_pairs=int(len(org_group)),
                interpretation="Checks whether results are dominated by organism identity.",
            ))

        more = group[group["expected_detectability"] == "more_detectable"]
        less = group[group["expected_detectability"] == "less_detectable"]
        if not more.empty and not less.empty:
            contrast = _mean_finite(more["auc"]) - _mean_finite(less["auc"])
            rows.append(dict(
                check_type="mechanism_contrast",
                site=site,
                excluded="none",
                grouping="mechanism",
                value=round(contrast, 4),
                n_pairs=int(len(more) + len(less)),
                interpretation="Pre-specified more-minus-less detectable AUC contrast.",
            ))

    non_d = df[df["site"] != "DRIAMS-D"]
    if not non_d.empty:
        more = non_d[non_d["expected_detectability"] == "more_detectable"]
        less = non_d[non_d["expected_detectability"] == "less_detectable"]
        if not more.empty and not less.empty:
            rows.append(dict(
                check_type="excluding_DRIAMS-D",
                site="all_non_D",
                excluded="DRIAMS-D",
                grouping="mechanism",
                value=round(_mean_finite(more["auc"]) - _mean_finite(less["auc"]), 4),
                n_pairs=int(len(more) + len(less)),
                interpretation="Checks whether the mechanism result is one-site dominated.",
            ))

    active_pairs_only = [(organism, drug) for _, organism, drug in active_pairs]
    for organism, drug in active_pairs_only:
        loo = df[~((df["organism"] == organism) & (df["drug"] == drug))]
        more = loo[loo["expected_detectability"] == "more_detectable"]
        less = loo[loo["expected_detectability"] == "less_detectable"]
        if more.empty or less.empty:
            continue
        rows.append(dict(
            check_type="leave_one_pair_out",
            site="all_sites",
            excluded=f"{organism}/{drug}",
            grouping="mechanism",
            value=round(_mean_finite(more["auc"]) - _mean_finite(less["auc"]), 4),
            n_pairs=int(len(more) + len(less)),
            interpretation="Checks whether a single pair drives the mechanism contrast.",
        ))
    return rows


def build_mechanism_confound_markdown(rows):
    lines = [
        "# Mechanism Confound Checks\n\n",
        "These checks are designed for reviewer questions: whether the mechanism "
        "contrast is mostly sample size, prevalence, organism identity, one site, "
        "or one organism-drug pair. They are sensitivity analyses, not causal proof.\n\n",
    ]
    lines.append(dataframe_to_markdown(pd.DataFrame(rows)))
    return "".join(lines)


# ═══════════════════════════════════════════════════════════════════════════════
# §16  ABLATION TRAINING  (Task 4) — calibrated BEFORE safety gate
# ═══════════════════════════════════════════════════════════════════════════════

def _train_ablation_variant(name, train_s, val_s, target_site_data, mae_model,
                             use_mae=True, use_dann=True, use_film=True,
                             active_pairs=None, early_stop="macro"):
    assert_source_only(f"_train_ablation_variant({name})")
    print(f"\n[Ablation: {name}]  MAE={use_mae}  DANN={use_dann}  "
          f"FiLM={use_film}  early_stop={early_stop}")
    if active_pairs is None:
        active_pairs = [(i, org, drug) for i, (org, drug) in enumerate(ORGANISM_DRUG_PAIRS)]
    train_loader  = make_loader(train_s, training=True)
    val_loader    = make_loader(val_s,   training=False)
    dom_loader    = make_domain_loader(target_site_data) if use_dann else None
    ModelClass    = MALDICNN if use_film else MALDICNNNoFiLM
    res_criterion = nn.BCEWithLogitsLoss()
    dom_criterion = nn.CrossEntropyLoss()
    trained       = []

    for seed in range(N_SEEDS_ABL):
        torch.manual_seed(seed + 200); np.random.seed(seed + 200)
        model = ModelClass(n_sites=N_SITES, n_organisms=N_ORGANISMS).to(DEVICE)
        if use_mae and mae_model is not None:
            load_mae_weights(mae_model, model)
        optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-4)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS_ABL)
        best_auc, patience_ct, best_state = 0.0, 0, None

        for epoch in range(1, EPOCHS_ABL + 1):
            if use_dann and dom_loader is not None:
                p   = (epoch - 1) / max(EPOCHS_ABL - 1, 1)
                lam = GRL_LAMBDA_MAX * (2.0 / (1.0 + math.exp(-10.0 * p)) - 1.0)
                tr_loss, _ = run_epoch_dann(model, train_loader, dom_loader,
                                             optimizer, res_criterion, dom_criterion, lam)
            else:
                tr_loss, _ = _run_epoch_no_dann(model, train_loader, optimizer, res_criterion)
            scheduler.step()
            val_logits, val_labels, val_orgs = evaluate_val(model, val_loader)
            auc_dict = compute_per_pair_auc(val_logits, val_labels, val_orgs, active_pairs)
            macro = macro_mean_auc(auc_dict)
            primary = auc_dict.get(PRIMARY_PAIR_IDX, macro)
            if early_stop == "primary":
                stop_val = primary
            else:
                stop_val = macro
            if math.isnan(stop_val):
                stop_val = macro if early_stop == "primary" else primary
            if stop_val > best_auc:
                best_auc    = stop_val
                best_state  = {k: v.cpu().clone() for k, v in model.state_dict().items()}
                patience_ct = 0
            else:
                patience_ct += 1
            if patience_ct >= PATIENCE_ABL:
                break

        m = ModelClass(n_sites=N_SITES, n_organisms=N_ORGANISMS).to(DEVICE)
        m.load_state_dict(best_state or {k: v.cpu() for k, v in model.state_dict().items()})
        m.eval(); trained.append(m)
        print(f"  [{name}] seed {seed+1}/{N_SEEDS_ABL}  best_auc={best_auc:.3f}")
    return trained


# ═══════════════════════════════════════════════════════════════════════════════
# §17  MODEL TRAINING  (load or train)
# ═══════════════════════════════════════════════════════════════════════════════

def select_seed_ensemble(trained_models, seed_scores, policy=SEED_POLICY,
                         top_k=TOP_K_SEEDS, min_macro_auc=MIN_SEED_MACRO_AUC):
    """Optionally keep only the strongest seeds by validation macro AUC."""
    if policy == "all" or len(trained_models) <= 1:
        return trained_models, [s.get("index", i) for i, s in enumerate(seed_scores)]

    finite_scores = [
        s for s in seed_scores
        if math.isfinite(float(s.get("val_macro", float("nan"))))
    ]
    if not finite_scores:
        print("  Seed QC: no finite macro-AUC scores; keeping all seeds")
        return trained_models, [s.get("index", i) for i, s in enumerate(seed_scores)]

    for score in finite_scores:
        if score["val_macro"] < NEAR_CHANCE_MACRO_AUC:
            print(f"  WARNING: seed {score['index']} is near chance "
                  f"(val_macro={score['val_macro']:.3f}); consider excluding it.")

    if policy == "threshold":
        chosen = [s for s in finite_scores if s["val_macro"] >= min_macro_auc]
        if not chosen:
            chosen = finite_scores
    elif policy == "topk":
        eligible = [s for s in finite_scores if s["val_macro"] >= min_macro_auc]
        pool = eligible if eligible else finite_scores
        chosen = sorted(pool, key=lambda s: s["val_macro"], reverse=True)[:max(1, min(top_k, len(pool)))]
    else:
        print(f"  Seed QC: unknown SEED_POLICY={policy!r}; keeping all seeds")
        return trained_models, [s.get("index", i) for i, s in enumerate(seed_scores)]

    chosen_indices = [int(s["index"]) for s in chosen]
    model_by_index = {
        int(score.get("index", i)): trained_models[i]
        for i, score in enumerate(seed_scores)
        if i < len(trained_models)
    }
    selected = [model_by_index[i] for i in chosen_indices if i in model_by_index]
    if not selected:
        print("  Seed QC: selected indices did not map to models; keeping all seeds")
        return trained_models, [s.get("index", i) for i, s in enumerate(seed_scores)]

    pretty = ", ".join(f"seed {s['index']} macro={s['val_macro']:.3f}" for s in chosen)
    print(f"  Seed QC ({policy}): keeping {len(selected)}/{len(trained_models)} seeds: {pretty}")
    return selected, chosen_indices


def train_all_seeds(train_s, val_s, target_site_data, ckpt_dir, mae_ckpt_path=None,
                    early_stop="macro", seed_policy=SEED_POLICY,
                    top_k_seeds=TOP_K_SEEDS,
                    min_seed_macro_auc=MIN_SEED_MACRO_AUC,
                    strict_mae_source_only=False,
                    drug_conditioning="task_id"):
    """
    Load checkpoints if they exist in ckpt_dir; otherwise train from scratch.
    Returns (trained_models, mae_model, seed_scores, selected_seed_indices).
    """
    ckpt_dir = Path(ckpt_dir)
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    active_pairs = [(i, org, drug) for i, (org, drug) in enumerate(ORGANISM_DRUG_PAIRS)]

    # ── MAE ──────────────────────────────────────────────────────────────────
    mae_model   = MAEWrapper().to(DEVICE)
    cached_mae  = ckpt_dir / "maldi_mae_pretrained.pt"

    if mae_ckpt_path and Path(mae_ckpt_path).exists():
        mae_model.load_state_dict(torch.load(mae_ckpt_path, map_location=DEVICE))
        print(f"MAE loaded from {mae_ckpt_path}")
    elif cached_mae.exists():
        mae_model.load_state_dict(torch.load(cached_mae, map_location=DEVICE))
        print(f"MAE loaded from cache: {cached_mae}")
    else:
        # Exclude test_s — bug fix vs multidrug_expansion.py line 1101
        mae_paths = list({p for p, _, _ in train_s + val_s})
        if not strict_mae_source_only:
            for sp in target_site_data.values():
                mae_paths.extend(p for p, _, _ in sp)
        mae_paths = list(set(mae_paths))
        print(f"\nMAE pre-training: {len(mae_paths)} spectra, {MAE_EPOCHS} epochs")
        pretrain_mae(mae_model, mae_paths)
        torch.save(mae_model.state_dict(), cached_mae)
        print(f"  MAE saved → {cached_mae}")

    # ── DANN seeds ────────────────────────────────────────────────────────────
    train_loader = make_loader(train_s, training=True)
    val_loader   = make_loader(val_s,   training=False)
    dom_loader   = make_domain_loader(target_site_data)
    if dom_loader is None:
        print("  WARNING: DANN disabled because no external target-site spectra were loaded.")
    res_crit     = nn.BCEWithLogitsLoss()
    dom_crit     = nn.CrossEntropyLoss()
    trained      = []
    seed_scores  = []

    for seed in range(N_SEEDS):
        seed_ckpt = ckpt_dir / f"maldi_amr_seed{seed}.pt"
        if seed_ckpt.exists():
            state = torch.load(seed_ckpt, map_location=DEVICE)
            if any(k.startswith("module.") for k in state):
                state = {k.removeprefix("module."): v for k, v in state.items()}
            m = create_maldi_model(
                n_sites=N_SITES, n_organisms=N_ORGANISMS,
                drug_conditioning=drug_conditioning).to(DEVICE)
            try:
                m.load_state_dict(state)
            except RuntimeError as exc:
                print(f"  Seed {seed}: {seed_ckpt.name} incompatible with current "
                      f"pair profile; retraining ({exc})")
            else:
                m.eval(); trained.append(m)
                val_logits, val_labels, val_orgs = evaluate_val(m, val_loader)
                auc_dict = compute_per_pair_auc(val_logits, val_labels, val_orgs, active_pairs)
                val_macro = macro_mean_auc(auc_dict)
                seed_scores.append(dict(index=seed, val_macro=val_macro, source="loaded"))
                print(f"  Seed {seed}: loaded {seed_ckpt.name}  val_macro={val_macro:.3f}")
                continue

        print(f"\n─── Seed {seed+1}/{N_SEEDS} ───")
        torch.manual_seed(seed); np.random.seed(seed)
        model     = create_maldi_model(
            n_sites=N_SITES, n_organisms=N_ORGANISMS,
            drug_conditioning=drug_conditioning).to(DEVICE)
        load_mae_weights(mae_model, model)
        optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=1e-4)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=EPOCHS)
        best_auc, patience_ct, best_state = 0.0, 0, None

        for epoch in range(1, EPOCHS + 1):
            p   = (epoch - 1) / max(EPOCHS - 1, 1)
            lam = GRL_LAMBDA_MAX * (2.0 / (1.0 + math.exp(-10.0 * p)) - 1.0)
            if dom_loader is not None:
                tr_loss, tr_acc = run_epoch_dann(model, train_loader, dom_loader,
                                                  optimizer, res_crit, dom_crit, lam)
            else:
                tr_loss, tr_acc = _run_epoch_no_dann(model, train_loader, optimizer, res_crit)
            scheduler.step()

            val_logits, val_labels, val_orgs = evaluate_val(model, val_loader)
            auc_dict = compute_per_pair_auc(val_logits, val_labels, val_orgs, active_pairs)
            macro    = macro_mean_auc(auc_dict)
            primary  = auc_dict.get(PRIMARY_PAIR_IDX, macro)
            if early_stop == "primary":
                stop_val = primary
            else:
                stop_val = macro
            if math.isnan(stop_val):
                stop_val = macro if early_stop == "primary" else primary

            marker = ""
            if stop_val > best_auc:
                best_auc    = stop_val; patience_ct = 0; marker = " *"
                best_state  = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            else:
                patience_ct += 1

            if epoch % CHECK_INTERVAL == 0:
                checkpoint_breakdown(model, val_loader, epoch, tr_loss, tr_acc, active_pairs)
            else:
                print(f"  Ep{epoch:3d}  loss={tr_loss:.4f}  primary={primary:.3f}  "
                      f"macro={macro:.3f}  λ={lam:.3f}{marker}", flush=True)

            if patience_ct >= PATIENCE:
                print(f"  Early stop at epoch {epoch}"); break

        if best_state is None:
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
        m = create_maldi_model(
            n_sites=N_SITES, n_organisms=N_ORGANISMS,
            drug_conditioning=drug_conditioning).to(DEVICE)
        m.load_state_dict(best_state); m.eval(); trained.append(m)
        seed_scores.append(dict(index=seed, val_macro=best_auc, source="trained"))
        torch.save(best_state, seed_ckpt)
        print(f"  Seed {seed+1} best_macro_auc={best_auc:.3f} → {seed_ckpt.name}")

    trained, selected_seed_indices = select_seed_ensemble(
        trained, seed_scores, policy=seed_policy,
        top_k=top_k_seeds, min_macro_auc=min_seed_macro_auc)
    return trained, mae_model, seed_scores, selected_seed_indices


def _slugify(text):
    return "".join(c.lower() if c.isalnum() else "_" for c in str(text)).strip("_")


def _binary_metrics_at_threshold(probs, labels, threshold=0.5):
    n = len(labels)
    n_r = int(labels.sum()) if n else 0
    n_s = n - n_r
    pred = (probs >= threshold).astype(float)
    tp = float(((pred == 1) & (labels == 1)).sum())
    tn = float(((pred == 0) & (labels == 0)).sum())
    return dict(
        auc=round(_safe_auc(labels, probs), 4),
        aupr=round(_safe_aupr(labels, probs), 4),
        sens=round(tp / n_r if n_r > 0 else float("nan"), 4),
        spec=round(tn / n_s if n_s > 0 else float("nan"), 4),
        n=n,
        n_r=n_r,
    )


def run_leave_one_drug_out_experiment(args, train_s, val_s, test_s,
                                      target_site_data, active_pairs,
                                      out_dir, base_ckpt_dir):
    """
    Optional zero-shot drug diagnostic for the E.coli panel.
    Held-out drug labels are excluded from source train/val and target UDA inputs.
    Metrics use raw T=1.0 probabilities and threshold=0.5 because no held-out
    validation labels are available for calibration.
    """
    rows = []
    eval_sites = {}
    if test_s:
        eval_sites[f"A-{TEST_YEAR}"] = test_s
    eval_sites.update(target_site_data)

    mae_ckpt = Path(base_ckpt_dir) / "maldi_mae_pretrained.pt"
    for heldout_id, organism, drug in active_pairs:
        loo_train = [s for s in train_s if s[2] != heldout_id]
        loo_val = [s for s in val_s if s[2] != heldout_id]
        if not loo_train or not loo_val:
            rows.append(dict(
                heldout_org_id=heldout_id, organism=organism, drug=drug,
                site="all", skipped=True,
                reason="empty train or validation split after holding drug out",
            ))
            continue

        loo_target = {
            site: [s for s in samples if s[2] != heldout_id]
            for site, samples in target_site_data.items()
        }
        slug = _slugify(f"{heldout_id}_{drug}_{args.drug_conditioning}")
        loo_ckpt_dir = Path(base_ckpt_dir) / f"leave_one_drug_out_{slug}"
        models, _, seed_scores, selected = train_all_seeds(
            loo_train, loo_val, loo_target, loo_ckpt_dir,
            mae_ckpt_path=str(mae_ckpt) if mae_ckpt.exists() else None,
            early_stop=args.early_stop,
            seed_policy=args.seed_policy,
            top_k_seeds=args.top_k_seeds,
            min_seed_macro_auc=args.min_seed_macro_auc,
            strict_mae_source_only=args.strict_mae_source_only,
            drug_conditioning=args.drug_conditioning)

        for site_label, site_samples in eval_sites.items():
            heldout_samples = [s for s in site_samples if s[2] == heldout_id]
            if not heldout_samples:
                continue
            probs, labels, _ = ensemble_predict(
                models, heldout_samples, temperature=1.0, n_passes=TTA_PASSES)
            row = dict(
                heldout_org_id=heldout_id,
                organism=organism,
                drug=drug,
                site=site_label,
                drug_conditioning=args.drug_conditioning,
                skipped=False,
                selected_seed_indices=";".join(map(str, selected)),
                seed_val_macro=";".join(
                    f"{float(s.get('val_macro', float('nan'))):.4f}"
                    for s in seed_scores),
                threshold=0.5,
            )
            row.update(_binary_metrics_at_threshold(probs, labels, threshold=0.5))
            rows.append(row)
    return rows


# ═══════════════════════════════════════════════════════════════════════════════
# §18  MAIN EXPERIMENT RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

def run_experiment(args):
    global MAE_EPOCHS
    global _EXTERNAL_EVAL_OPEN
    _EXTERNAL_EVAL_OPEN = False

    pair_profile = getattr(args, "pair_profile", None)
    if pair_profile is None:
        legacy_mode = getattr(args, "mode", None)
        pair_profile = "run14" if legacy_mode == "run14" else (
            "gram_negative6" if legacy_mode == "multidrug" else NOTEBOOK_PAIR_PROFILE)
    args.pair_profile = pair_profile

    init_config(pair_profile)
    data_root = args.data_root or DATA_ROOT
    args.early_stop = getattr(args, "early_stop", NOTEBOOK_EARLY_STOP)
    args.seed_policy = getattr(args, "seed_policy", SEED_POLICY)
    args.top_k_seeds = getattr(args, "top_k_seeds", TOP_K_SEEDS)
    args.min_seed_macro_auc = getattr(args, "min_seed_macro_auc", MIN_SEED_MACRO_AUC)
    args.strict_mae_source_only = getattr(
        args, "strict_mae_source_only", NOTEBOOK_STRICT_MAE_SOURCE_ONLY)
    args.prevalence_shift = "none" if getattr(args, "no_prevalence_shift", False) else getattr(
        args, "prevalence_shift", NOTEBOOK_PREVALENCE_SHIFT)
    args.mae_epochs = int(getattr(args, "mae_epochs", NOTEBOOK_MAE_EPOCHS))
    args.with_random_cv = bool(getattr(args, "with_random_cv", NOTEBOOK_WITH_RANDOM_CV))
    args.drug_conditioning = getattr(
        args, "drug_conditioning", NOTEBOOK_DRUG_CONDITIONING)
    args.with_leave_one_drug_out = bool(getattr(
        args, "with_leave_one_drug_out", NOTEBOOK_WITH_LEAVE_ONE_DRUG_OUT))
    if args.mae_epochs <= 0:
        raise ValueError("--mae-epochs must be a positive integer")
    if args.drug_conditioning not in DRUG_CONDITIONING_CHOICES:
        raise ValueError(f"--drug-conditioning must be one of {DRUG_CONDITIONING_CHOICES}")
    MAE_EPOCHS = args.mae_epochs

    # ── Directory + logging setup ─────────────────────────────────────────────
    out_dir = Path(args.output_dir) / args.experiment
    for sub in ("metrics", "paper", "models", "logs"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "logs" / "run.log"
    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s  %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger()
    log.info(f"=== Experiment: {args.experiment}  pair_profile={args.pair_profile} ===")
    log.info(f"MAE epochs: {MAE_EPOCHS}")
    log.info(f"Drug conditioning: {args.drug_conditioning}")

    # ── Active pairs (MIN_R_TRAIN screening in multidrug mode) ────────────────
    active_pairs = screen_active_pairs(data_root, args.pair_profile)

    # ── Data loading ──────────────────────────────────────────────────────────
    log.info(f"\nLoading {TRAIN_SITE} ...")
    spectrum_index_A = build_spectrum_index(data_root, TRAIN_SITE)
    all_samples = load_all_organisms(data_root, TRAIN_SITE, active_pairs, spectrum_index_A)
    if not all_samples:
        raise RuntimeError("No samples found.")
    if len(all_samples) < MIN_TOTAL_SOURCE_SAMPLES_WARN:
        log.warning(f"Suspiciously small source sample count: {len(all_samples)} "
                    f"(< {MIN_TOTAL_SOURCE_SAMPLES_WARN}). Check DATA_ROOT, "
                    f"{TRAIN_SITE}/binned_6000, and label CSV paths.")
    train_s, val_s, test_s = make_split(all_samples, active_pairs)

    log.info("Loading external sites ...")
    log.info("External target spectra are used only as unlabeled unsupervised domain adaptation inputs "
             "for DANN/BN adaptation; labels are "
             "reserved for evaluation metrics.")
    target_site_data = {}
    for site in TEST_SITES:
        if not (Path(data_root) / site).exists():
            log.info(f"  {site}: not found, skipping")
            continue
        sp = load_all_organisms(data_root, site, active_pairs)
        if sp:
            target_site_data[site] = [(p, l, o) for p, l, o, _ in sp]

    # ── Load or train models ──────────────────────────────────────────────────
    ckpt_dir = args.ckpt_dir or str(out_dir / "models")
    log.info(f"\nCheckpoint dir: {ckpt_dir}")
    models, mae_model, seed_scores, selected_seed_indices = train_all_seeds(
        train_s, val_s, target_site_data, ckpt_dir,
        mae_ckpt_path=None,
        early_stop=args.early_stop,
        seed_policy=args.seed_policy,
        top_k_seeds=args.top_k_seeds,
        min_seed_macro_auc=args.min_seed_macro_auc,
        strict_mae_source_only=args.strict_mae_source_only,
        drug_conditioning=args.drug_conditioning)

    # ── Calibration (MUST happen before open_external_eval) ───────────────────
    log.info("\nCalibration (source val only) ...")

    T_opt = None
    base_thresholds = None

    if args.calib_json and Path(args.calib_json).exists():
        with open(args.calib_json) as f:
            c = json.load(f)
        T_opt = float(c["temperature"])
        log.info(f"  Loaded T={T_opt:.4f} from {args.calib_json}")
        if "thresholds_per_pair" in c:
            base_thresholds = {int(k): float(v) for k, v in c["thresholds_per_pair"].items()}
            log.info(f"  Loaded per-pair thresholds from JSON")

    if T_opt is None:
        T_opt = learn_ensemble_temperature(models, val_s)
    if base_thresholds is None:
        log.info("  Computing per-pair Youden thresholds ...")
        base_thresholds = find_thresholds_per_pair(models, val_s, T_opt, active_pairs)

    # Pre-compute per-site threshold maps
    site_threshold_maps = {f"A-{TEST_YEAR}": base_thresholds}
    for site in TEST_SITES:
        if site in target_site_data:
            log.info(f"  Prevalence shift for {site}:")
            site_threshold_maps[site] = site_thresholds_for(
                base_thresholds, site, active_pairs,
                prevalence_shift_mode=args.prevalence_shift)

    # ── Ablations (trained and calibrated BEFORE gate) ────────────────────────
    ablation_data: dict = {}
    if args.with_ablation:
        log.info("\n=== Pre-gate: Ablation Training ===")
        ablation_specs = [
            ("no_mae",  False, True,  True),
            ("no_dann", True,  False, True),
            ("no_film", True,  True,  False),
        ]
        for abl_name, use_mae, use_dann, use_film in ablation_specs:
            try:
                abl_models = _train_ablation_variant(
                    abl_name, train_s, val_s, target_site_data, mae_model,
                    use_mae=use_mae, use_dann=use_dann, use_film=use_film,
                    active_pairs=active_pairs, early_stop=args.early_stop)
                abl_T      = learn_ensemble_temperature(abl_models, val_s)
                abl_thresh = find_thresholds_per_pair(abl_models, val_s, abl_T, active_pairs)
                # Pre-compute site threshold maps for ablation
                abl_site_thresh = {f"A-{TEST_YEAR}": abl_thresh}
                for site in TEST_SITES:
                    if site in target_site_data:
                        abl_site_thresh[site] = site_thresholds_for(
                            abl_thresh, site, active_pairs,
                            prevalence_shift_mode=args.prevalence_shift)
                ablation_data[abl_name] = (abl_models, abl_T, abl_site_thresh)
                # Save ablation checkpoints
                for i, m in enumerate(abl_models):
                    torch.save(m.state_dict(),
                               out_dir / "models" / f"ablation_{abl_name}_seed{i}.pt")
            except Exception:
                log.warning(f"  Ablation {abl_name} failed:\n" + traceback.format_exc())

    # ── Open safety gate — no calibration after this ──────────────────────────
    open_external_eval()
    use_bn = not args.no_bn_adapt

    # ── Save config ───────────────────────────────────────────────────────────
    with open(out_dir / "config.json", "w") as f:
        json.dump(dict(
            experiment=args.experiment, pair_profile=args.pair_profile,
            temperature=T_opt,
            thresholds_per_pair={str(k): round(v, 5) for k, v in base_thresholds.items()},
            n_seeds=N_SEEDS, tta_passes=TTA_PASSES,
            seed_policy=args.seed_policy,
            top_k_seeds=args.top_k_seeds,
            min_seed_macro_auc=args.min_seed_macro_auc,
            early_stop=args.early_stop,
            mae_epochs=args.mae_epochs,
            drug_conditioning=args.drug_conditioning,
            with_leave_one_drug_out=args.with_leave_one_drug_out,
            with_random_cv=args.with_random_cv,
            prevalence_shift=args.prevalence_shift,
            strict_mae_source_only=args.strict_mae_source_only,
            no_saliency=getattr(args, "no_saliency", NOTEBOOK_NO_SALIENCY),
            seed_scores=seed_scores,
            selected_seed_indices=selected_seed_indices,
            use_bn_adapt=use_bn, data_root=data_root, ckpt_dir=ckpt_dir,
            pairs=[(o, d) for o, d in ORGANISM_DRUG_PAIRS],
            active_pairs=[(oid, o, d) for oid, o, d in active_pairs],
        ), f, indent=2)

    seed_summary_rows = []
    selected_seed_set = set(selected_seed_indices)
    for score in seed_scores:
        row = dict(score)
        row["selected"] = int(row.get("index", -1)) in selected_seed_set
        row["near_chance_warning"] = (
            math.isfinite(float(row.get("val_macro", float("nan"))))
            and float(row.get("val_macro")) < NEAR_CHANCE_MACRO_AUC
        )
        seed_summary_rows.append(row)
    pd.DataFrame(seed_summary_rows).to_csv(out_dir / "metrics" / "seed_summary.csv", index=False)

    # Evaluation sites
    eval_sites: dict = {}
    if test_s:
        eval_sites[f"A-{TEST_YEAR}"] = (test_s, False)
    else:
        log.info(f"  A-{TEST_YEAR}: no same-site test samples, skipping")
    for site, sp in target_site_data.items():
        eval_sites[site] = (sp, True)

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 1 — Multi-organism CNN evaluation
    # ═════════════════════════════════════════════════════════════════════════
    log.info("\n=== Task 1: Multi-Organism Evaluation (CNN) ===")
    all_multi_rows, macro_rows, micro_rows = [], [], []
    cnn_cache: dict = {}

    for site_label, (site_sp, is_ext) in eval_sites.items():
        thresh_map = site_threshold_maps.get(site_label, base_thresholds)
        log.info(f"  {site_label}  n={len(site_sp)}")

        with adapted_batchnorm(models, site_sp, use_adapt=is_ext and use_bn):
            probs, labels, orgs = ensemble_predict(models, site_sp, T_opt, TTA_PASSES)

        cnn_cache[site_label] = (probs, labels, orgs)
        rows, macro_auc = compute_multi_metrics(probs, labels, orgs, thresh_map, site_label, active_pairs)
        all_multi_rows.extend(rows)

        for r in rows:
            if r["organism"] != "__micro__":
                log.info(f"    {r['organism'][:22]:22s}/{r['drug']:14s}  "
                         f"AUC={r['auc']:.3f}  AUPR={r['aupr']:.3f}  "
                         f"Sens={r['sens']:.3f}  Spec={r['spec']:.3f}")
        log.info(f"    Macro AUC={macro_auc:.3f}")

        macro_rows.append(dict(site=site_label, macro_auc=round(macro_auc, 4),
                               n_pairs=len(active_pairs)))
        micro_r = next((r for r in rows if r["organism"] == "__micro__"), None)
        if micro_r:
            micro_rows.append(dict(site=site_label, **{k: micro_r[k] for k in
                ("auc", "aupr", "sens", "spec", "n", "n_r")}))

    pd.DataFrame(all_multi_rows).to_csv(out_dir / "metrics" / "multi_metrics_summary.csv", index=False)
    pd.DataFrame(macro_rows).to_csv(out_dir / "metrics" / "multi_macro_summary.csv", index=False)
    pd.DataFrame(micro_rows).to_csv(out_dir / "metrics" / "multi_micro_summary.csv", index=False)
    log.info("  Saved multi_metrics_summary / macro / micro CSVs")

    per_organism_rows = compute_per_organism_summary(all_multi_rows)
    pd.DataFrame(per_organism_rows).to_csv(
        out_dir / "metrics" / "per_organism_summary.csv", index=False)
    log.info("  Saved per_organism_summary.csv")

    run14_overlap_rows = []
    for site_label, (probs, labels, orgs) in cnn_cache.items():
        thresh_map = site_threshold_maps.get(site_label, base_thresholds)
        run14_overlap_rows.append(compute_subset_metrics(
            probs, labels, orgs, thresh_map, site_label, active_pairs,
            RUN14_OVERLAP_PAIRS, subset_name="run14_overlap"))
    pd.DataFrame(run14_overlap_rows).to_csv(
        out_dir / "metrics" / "run14_overlap_summary.csv", index=False)
    log.info("  Saved run14_overlap_summary.csv")

    threshold_sanity_rows = compute_threshold_sanity_report(
        all_multi_rows, site_threshold_maps, base_thresholds, active_pairs)
    pd.DataFrame(threshold_sanity_rows).to_csv(
        out_dir / "metrics" / "threshold_sanity_report.csv", index=False)
    log.info(f"  Saved threshold_sanity_report.csv ({len(threshold_sanity_rows)} flags)")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 5 — LightGBM baselines
    # ═════════════════════════════════════════════════════════════════════════
    lgbm_single, lgbm_multi, lgbm_all_rows = {}, None, []
    if not args.no_lgbm:
        log.info("\n=== Task 5: LightGBM Baselines ===")
        try:
            log.info("  Training single-task LightGBM ...")
            lgbm_single = train_lgbm_singletask(train_s, val_s, active_pairs)
            log.info("  Training multi-task LightGBM ...")
            lgbm_multi  = train_lgbm_multitask(train_s, val_s)
            for oid, m in lgbm_single.items():
                with open(out_dir / "models" / f"lgbm_single_org{oid}.pkl", "wb") as f:
                    pickle.dump(m, f)
            if lgbm_multi:
                with open(out_dir / "models" / "lgbm_multi.pkl", "wb") as f:
                    pickle.dump(lgbm_multi, f)
            for site_label, (site_sp, _) in eval_sites.items():
                rows = evaluate_lgbm_site(lgbm_single, lgbm_multi, site_sp, active_pairs)
                for r in rows:
                    r["site"] = site_label
                lgbm_all_rows.extend(rows)
            pd.DataFrame(lgbm_all_rows).to_csv(out_dir / "metrics" / "lgbm_results.csv", index=False)
            log.info("  Saved lgbm_results.csv")
        except Exception:
            log.warning("  LightGBM failed:\n" + traceback.format_exc())

    random_cv_rows = []
    if args.with_random_cv and lgbm_all_rows:
        log.info("\n=== Temporal vs Random CV Diagnostic ===")
        try:
            random_cv_rows = run_random_cv_inflation_analysis(
                all_samples, lgbm_all_rows, active_pairs)
            pd.DataFrame(random_cv_rows).to_csv(
                out_dir / "metrics" / "temporal_vs_random_cv.csv", index=False)
            with open(out_dir / "paper" / "temporal_vs_random_cv.md", "w") as f:
                f.write(build_temporal_vs_random_cv_markdown(random_cv_rows))
            log.info("  Saved temporal_vs_random_cv.csv and temporal_vs_random_cv.md")
        except Exception:
            log.warning("  Temporal-vs-random CV diagnostic failed:\n" + traceback.format_exc())
    else:
        log.info("  Temporal-vs-random CV diagnostic skipped")

    log.info("\n=== Evaluation Critique Table ===")
    evaluation_critique_rows = build_evaluation_critique_table(
        macro_rows, lgbm_all_rows, random_cv_rows)
    pd.DataFrame(evaluation_critique_rows).to_csv(
        out_dir / "metrics" / "evaluation_critique_table.csv", index=False)
    with open(out_dir / "paper" / "evaluation_critique_table.md", "w") as f:
        f.write(build_evaluation_critique_markdown(evaluation_critique_rows))
    log.info("  Saved evaluation_critique_table.csv and evaluation_critique_table.md")

    log.info("\n=== Paper Core Results Table ===")
    core_results_rows = build_core_results_table(all_multi_rows, lgbm_all_rows, active_pairs)
    core_results_df = pd.DataFrame(core_results_rows)
    core_results_df.to_csv(out_dir / "metrics" / "core_results_table.csv", index=False)
    write_markdown_table(
        core_results_df,
        out_dir / "paper" / "core_results_table.md",
        title="Core Results Table",
        note=(
            "CNN is the clinical multi-head MAE+FiLM+DANN model. "
            "LightGBM columns are included as classical baselines; "
            "delta is CNN AUC minus the stronger LightGBM AUC for that row."
        ),
        columns=[
            "row_type", "site", "organism", "drug", "n", "n_r",
            "cnn_auc", "cnn_aupr", "cnn_sens", "cnn_spec",
            "lgbm_single_auc", "lgbm_multi_auc",
            "best_lgbm_auc", "cnn_minus_best_lgbm_auc",
        ],
    )
    log.info("  Saved core_results_table.csv and core_results_table.md")

    log.info("\n=== Mechanism Detectability Framing ===")
    mechanism_pair_rows, mechanism_summary_rows = compute_mechanism_detectability_summary(
        all_multi_rows, active_pairs)
    pd.DataFrame(mechanism_pair_rows).to_csv(
        out_dir / "metrics" / "mechanism_pair_metrics.csv", index=False)
    pd.DataFrame(mechanism_summary_rows).to_csv(
        out_dir / "metrics" / "mechanism_detectability_summary.csv", index=False)
    with open(out_dir / "paper" / "mechanism_framing.md", "w") as f:
        f.write(build_mechanism_framing_markdown(
            mechanism_pair_rows, mechanism_summary_rows, active_pairs))
    log.info("  Saved mechanism_pair_metrics.csv, mechanism_detectability_summary.csv, and mechanism_framing.md")

    log.info("\n=== Statistical Reporting ===")
    stat_report = compute_prediction_statistical_report(cnn_cache, active_pairs)
    pd.DataFrame(stat_report["pair_ci_rows"]).to_csv(
        out_dir / "metrics" / "auc_bootstrap_ci.csv", index=False)
    pd.DataFrame(stat_report["macro_ci_rows"]).to_csv(
        out_dir / "metrics" / "macro_bootstrap_ci.csv", index=False)
    pd.DataFrame(stat_report["contrast_ci_rows"]).to_csv(
        out_dir / "metrics" / "mechanism_contrast_bootstrap_ci.csv", index=False)
    pd.DataFrame(stat_report["permutation_rows"]).to_csv(
        out_dir / "metrics" / "mechanism_permutation_tests.csv", index=False)
    with open(out_dir / "paper" / "statistical_reporting.md", "w") as f:
        f.write(build_statistical_reporting_markdown(
            stat_report["pair_ci_rows"],
            stat_report["macro_ci_rows"],
            stat_report["contrast_ci_rows"],
            stat_report["permutation_rows"],
        ))
    log.info("  Saved bootstrap CI, mechanism permutation, and statistical_reporting artifacts")

    log.info("\n=== Mechanism Confound Checks ===")
    confound_rows = compute_mechanism_confound_checks(all_multi_rows, active_pairs)
    pd.DataFrame(confound_rows).to_csv(
        out_dir / "metrics" / "mechanism_confound_checks.csv", index=False)
    with open(out_dir / "paper" / "mechanism_confound_checks.md", "w") as f:
        f.write(build_mechanism_confound_markdown(confound_rows))
    log.info("  Saved mechanism_confound_checks.csv and mechanism_confound_checks.md")

    log.info("\n=== DRIAMS-D Community Shift Analysis ===")
    d_site_rows = compute_d_site_shift_analysis(all_multi_rows, active_pairs)
    pd.DataFrame(d_site_rows).to_csv(
        out_dir / "metrics" / "d_site_shift_analysis.csv", index=False)
    with open(out_dir / "paper" / "d_site_shift_analysis.md", "w") as f:
        f.write(build_d_site_shift_markdown(d_site_rows))
    log.info("  Saved d_site_shift_analysis.csv and d_site_shift_analysis.md")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 3 — Domain shift analysis
    # ═════════════════════════════════════════════════════════════════════════
    log.info("\n=== Task 3: Domain Shift Analysis ===")
    shift_rows = []
    site_a     = f"A-{TEST_YEAR}"
    if site_a in cnn_cache:
        p_a, l_a, o_a = cnn_cache[site_a]
        for org_id, organism, drug in active_pairs:
            mask_a = (o_a == org_id)
            auc_a  = _safe_auc(l_a[mask_a], p_a[mask_a]) if mask_a.any() else float("nan")
            for ext_site in TEST_SITES:
                if ext_site not in cnn_cache:
                    continue
                p_s, l_s, o_s = cnn_cache[ext_site]
                mask_s = (o_s == org_id)
                auc_s  = _safe_auc(l_s[mask_s], p_s[mask_s]) if mask_s.any() else float("nan")
                delta  = (auc_s - auc_a) if not (np.isnan(auc_s) or np.isnan(auc_a)) else float("nan")
                shift_rows.append(dict(
                    model="CNN", organism=organism, drug=drug,
                    site_a=site_a, site=ext_site,
                    auc_a=round(auc_a, 4), auc_site=round(auc_s, 4),
                    delta_auc=round(delta, 4),
                ))
    if lgbm_all_rows:
        lgbm_df = pd.DataFrame(lgbm_all_rows)
        for variant in ("lgbm_single", "lgbm_multi"):
            vdf  = lgbm_df[lgbm_df["model"] == variant]
            a_df = vdf[vdf["site"] == site_a]
            for org_id, organism, drug in active_pairs:
                ar = a_df[a_df["organism"] == organism]
                if ar.empty:
                    continue
                auc_a = float(ar["auc"].iloc[0])
                for ext_site in TEST_SITES:
                    sr = vdf[(vdf["site"] == ext_site) & (vdf["organism"] == organism)]
                    if sr.empty:
                        continue
                    shift_rows.append(dict(
                        model=variant, organism=organism, drug=drug,
                        site_a=site_a, site=ext_site,
                        auc_a=round(auc_a, 4), auc_site=round(float(sr["auc"].iloc[0]), 4),
                        delta_auc=round(float(sr["auc"].iloc[0]) - auc_a, 4),
                    ))
    pd.DataFrame(shift_rows).to_csv(out_dir / "metrics" / "domain_shift_analysis.csv", index=False)
    log.info("  Saved domain_shift_analysis.csv")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 4 — Ablation evaluation
    # ═════════════════════════════════════════════════════════════════════════
    if ablation_data:
        log.info("\n=== Task 4: Ablation Evaluation ===")
        ablation_rows = []
        for variant, (abl_models, abl_T, abl_site_thresh) in ablation_data.items():
            for site_label, (site_sp, is_ext) in eval_sites.items():
                thresh_map = abl_site_thresh.get(site_label, base_thresholds)
                try:
                    with adapted_batchnorm(abl_models, site_sp, use_adapt=is_ext and use_bn):
                        probs, labels, orgs = ensemble_predict(abl_models, site_sp, abl_T, n_passes=5)
                    rows, macro_auc = compute_multi_metrics(
                        probs, labels, orgs, thresh_map, site_label, active_pairs)
                    for r in rows:
                        ablation_rows.append(dict(ablation=variant, **r))
                    log.info(f"  [{variant}] {site_label}: macro={macro_auc:.3f}")
                except Exception:
                    log.warning(f"  [{variant}] {site_label} failed:\n" + traceback.format_exc())

        # Add full model as reference row
        for site_label in eval_sites:
            if site_label in cnn_cache:
                thresh_map = site_threshold_maps.get(site_label, base_thresholds)
                rows, _ = compute_multi_metrics(
                    *cnn_cache[site_label], thresh_map, site_label, active_pairs)
                for r in rows:
                    ablation_rows.append(dict(ablation="full_model", **r))

        pd.DataFrame(ablation_rows).to_csv(out_dir / "metrics" / "ablation_multi.csv", index=False)
        log.info("  Saved ablation_multi.csv")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 6 — Threshold tradeoff analysis
    # ═════════════════════════════════════════════════════════════════════════
    log.info("\n=== Task 6: Threshold Analysis ===")

    # Derive three strategies from source val probs (no gate violation — already open,
    # but this uses only val data to define strategies, not to calibrate the production model)
    val_probs_all, val_labels_all, val_orgs_all = ensemble_predict(
        models, val_s, T_opt, n_passes=1)

    # Strategy 1: per-pair Youden (already computed — base_thresholds)
    # Strategy 2: global Youden on E.coli val (ECOLI org_id=0 in both modes)
    ec_mask = (val_orgs_all == 0)
    if ec_mask.any() and val_labels_all[ec_mask].sum() > 0:
        t_youden_ec, _, _, _ = _best_youden_threshold(
            val_labels_all[ec_mask], val_probs_all[ec_mask])
        # Strategy 3: macro-mean val probs Youden
        t_youden_global, _, _, _ = _best_youden_threshold(val_labels_all, val_probs_all)
    else:
        t_youden_ec = t_youden_global = 0.5

    thresh_rows = []
    strategy_defs = {
        "per_pair_youden": base_thresholds,
        "global_ecoli_youden": {oid: t_youden_ec for oid, _, _ in active_pairs},
        "global_pooled_youden": {oid: t_youden_global for oid, _, _ in active_pairs},
    }

    for strat_name, base_t_map in strategy_defs.items():
        for site_label, (site_sp, is_ext) in eval_sites.items():
            if is_ext:
                t_map = site_thresholds_for(
                    base_t_map, site_label, active_pairs,
                    prevalence_shift_mode=args.prevalence_shift)
            else:
                t_map = base_t_map
            probs, labels, orgs = (
                cnn_cache[site_label] if site_label in cnn_cache
                else ensemble_predict(models, site_sp, T_opt, n_passes=5)
            )
            for org_id, organism, drug in active_pairs:
                mask = (orgs == org_id)
                if not mask.any():
                    continue
                p = probs[mask]; l = labels[mask]
                n = len(l); n_r = int(l.sum()); n_s = n - n_r
                thr  = t_map[org_id] if isinstance(t_map, dict) else float(t_map)
                pred = (p >= thr).astype(float)
                tp   = float(((pred == 1) & (l == 1)).sum())
                tn   = float(((pred == 0) & (l == 0)).sum())
                thresh_rows.append(dict(
                    strategy=strat_name, site=site_label,
                    organism=organism, drug=drug, threshold=round(thr, 5),
                    auc=round(_safe_auc(l, p), 4),
                    sens=round(tp / n_r if n_r > 0 else float("nan"), 4),
                    spec=round(tn / n_s if n_s > 0 else float("nan"), 4),
                    n=n, n_r=n_r,
                ))

    pd.DataFrame(thresh_rows).to_csv(
        out_dir / "metrics" / "threshold_tradeoff_multi.csv", index=False)
    log.info("  Saved threshold_tradeoff_multi.csv")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 8 — Saliency / interpretability
    # ═════════════════════════════════════════════════════════════════════════
    if getattr(args, "no_saliency", NOTEBOOK_NO_SALIENCY):
        log.info("\n=== Task 8: Saliency / Interpretability skipped ===")
    else:
        log.info("\n=== Task 8: Saliency / Interpretability ===")
        saliency_rows = compute_occlusion_saliency(
            models, val_s, T_opt, active_pairs, site_label="A-val")
        saliency_df = pd.DataFrame(saliency_rows)
        saliency_df.to_csv(out_dir / "metrics" / "saliency_top_bins.csv", index=False)
        write_markdown_table(
            saliency_df,
            out_dir / "paper" / "saliency_top_bins.md",
            title="Occlusion Saliency Top Spectral Windows",
            note=(
                "Each row occludes one fixed-width m/z window on source validation spectra "
                "and ranks windows by mean absolute change in ensemble resistance probability."
            ),
            columns=[
                "organism", "drug", "rank", "bin_start", "bin_end",
                "mz_start", "mz_end", "importance", "signed_delta",
                "n_samples", "n_r",
            ],
        )
        log.info("  Saved saliency_top_bins.csv and saliency_top_bins.md")

        saliency_by_site_rows = list(saliency_rows)
        for site_label, (site_sp, is_ext) in eval_sites.items():
            with adapted_batchnorm(models, site_sp, use_adapt=is_ext and use_bn):
                saliency_by_site_rows.extend(compute_occlusion_saliency(
                    models, site_sp, T_opt, active_pairs, site_label=site_label))
        saliency_by_site_df = pd.DataFrame(saliency_by_site_rows)
        saliency_by_site_df.to_csv(
            out_dir / "metrics" / "saliency_top_bins_by_site.csv", index=False)
        reference_site = f"A-{TEST_YEAR}" if f"A-{TEST_YEAR}" in eval_sites else "A-val"
        saliency_stability_rows = compute_saliency_stability_summary(
            saliency_by_site_rows, active_pairs, reference_site=reference_site)
        pd.DataFrame(saliency_stability_rows).to_csv(
            out_dir / "metrics" / "saliency_stability_summary.csv", index=False)
        with open(out_dir / "paper" / "saliency_stability.md", "w") as f:
            f.write(build_saliency_stability_markdown(saliency_stability_rows))
        log.info("  Saved saliency_top_bins_by_site.csv, saliency_stability_summary.csv, and saliency_stability.md")

        saliency_robust_including_rows, saliency_alignment_including_rows, saliency_null_including_rows = (
            compute_saliency_robustness_summary(
                saliency_by_site_rows, all_multi_rows, active_pairs,
                reference_site=reference_site, exclude_low_bins_below=0)
        )
        saliency_robust_excluding_rows, saliency_alignment_excluding_rows, saliency_null_excluding_rows = (
            compute_saliency_robustness_summary(
                saliency_by_site_rows, all_multi_rows, active_pairs,
                reference_site=reference_site,
                exclude_low_bins_below=SALIENCY_EXCLUDE_LOW_BINS_BELOW)
        )
        for row in saliency_robust_including_rows:
            row["low_bin_sensitivity"] = "including_low_bins"
        for row in saliency_robust_excluding_rows:
            row["low_bin_sensitivity"] = "excluding_low_bins"
        for row in saliency_alignment_including_rows:
            row["low_bin_sensitivity"] = "including_low_bins"
        for row in saliency_alignment_excluding_rows:
            row["low_bin_sensitivity"] = "excluding_low_bins"
        for row in saliency_null_including_rows:
            row["low_bin_sensitivity"] = "including_low_bins"
        for row in saliency_null_excluding_rows:
            row["low_bin_sensitivity"] = "excluding_low_bins"

        saliency_robust_rows = saliency_robust_including_rows + saliency_robust_excluding_rows
        saliency_alignment_rows = saliency_alignment_including_rows + saliency_alignment_excluding_rows
        saliency_null_rows = saliency_null_including_rows + saliency_null_excluding_rows

        pd.DataFrame(saliency_robust_rows).to_csv(
            out_dir / "metrics" / "saliency_stability_robustness.csv", index=False)
        pd.DataFrame(saliency_robust_including_rows).to_csv(
            out_dir / "metrics" / "saliency_stability_robustness_including_low_bins.csv", index=False)
        pd.DataFrame(saliency_robust_excluding_rows).to_csv(
            out_dir / "metrics" / "saliency_stability_robustness_excluding_low_bins.csv", index=False)
        pd.DataFrame(saliency_alignment_rows).to_csv(
            out_dir / "metrics" / "saliency_auc_alignment.csv", index=False)
        pd.DataFrame(saliency_null_rows).to_csv(
            out_dir / "metrics" / "saliency_null_control.csv", index=False)
        with open(out_dir / "paper" / "saliency_robustness.md", "w") as f:
            f.write(build_saliency_robustness_markdown(
                saliency_robust_rows, saliency_alignment_rows, saliency_null_rows))
        log.info("  Saved saliency robustness, AUC alignment, and null-control artifacts")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 8b — Leave-one-drug-out for drug fingerprint experiments
    # ═════════════════════════════════════════════════════════════════════════
    if args.with_leave_one_drug_out:
        log.info("\n=== Task 8b: Leave-One-Drug-Out Diagnostic ===")
        if args.drug_conditioning == "task_id":
            log.warning("  task_id conditioning has no chemical zero-shot prior; "
                        "use drug_id/morgan/drug_id+morgan comparisons carefully.")
        try:
            loodo_rows = run_leave_one_drug_out_experiment(
                args, train_s, val_s, test_s, target_site_data,
                active_pairs, out_dir, ckpt_dir)
            loodo_df = pd.DataFrame(loodo_rows)
            loodo_df.to_csv(out_dir / "metrics" / "leave_one_drug_out.csv", index=False)
            write_markdown_table(
                loodo_df,
                out_dir / "paper" / "leave_one_drug_out.md",
                title="Leave-One-Drug-Out Diagnostic",
                note=(
                    "Each held-out drug is excluded from source train/validation labels "
                    "and target UDA inputs. Metrics use uncalibrated T=1.0 probabilities "
                    "because held-out validation labels are intentionally unavailable."
                ),
                columns=[
                    "drug_conditioning", "site", "organism", "drug", "auc", "aupr",
                    "sens", "spec", "n", "n_r", "selected_seed_indices",
                ],
            )
            log.info("  Saved leave_one_drug_out.csv and leave_one_drug_out.md")
        except Exception:
            log.warning("  Leave-one-drug-out diagnostic failed:\n" + traceback.format_exc())

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 2 — Paper comparison table
    # ═════════════════════════════════════════════════════════════════════════
    log.info("\n=== Task 2: Multi vs Single Comparison ===")
    ext_sites = [s for s in TEST_SITES if s in target_site_data]

    def _get_org_auc_from_cache(site_label, org_id):
        if site_label not in cnn_cache:
            return float("nan")
        p, l, o = cnn_cache[site_label]
        mask = (o == org_id)
        return _safe_auc(l[mask], p[mask]) if mask.any() else float("nan")

    lines = [
        "# Multi-Organism vs Single-Task Comparison\n\n",
        "CNN uses joint multi-head training (DANN+FiLM+MAE). "
        "LightGBM baselines are single-site with no domain adaptation.\n\n",
    ]
    # Build per-organism columns for pair 0 (E.coli/Cipro) across sites
    all_cols = ["Model", f"Pair0(A-{TEST_YEAR})", "MacroAUC(A)",
                ] + [f"Pair0({s})" for s in ext_sites] + ["Cross-site(Pair0)"]
    lines.append("| " + " | ".join(all_cols) + " |\n")
    lines.append("|" + "|".join(["-"*(len(c)+2) for c in all_cols]) + "|\n")

    def _row(label, auc_a, mac_a, ext_aucs):
        cross = float(np.nanmean(ext_aucs)) if ext_aucs else float("nan")
        vals  = [label, f"{auc_a:.3f}", f"{mac_a:.3f}"]
        vals += [f"{v:.3f}" if not np.isnan(v) else "—" for v in ext_aucs]
        vals += [f"{cross:.3f}" if not np.isnan(cross) else "—"]
        return "| " + " | ".join(vals) + " |\n"

    macro_a = next((r["macro_auc"] for r in macro_rows if r["site"] == site_a), float("nan"))
    lines.append(_row(
        f"CNN ({args.pair_profile})",
        _get_org_auc_from_cache(site_a, 0), macro_a,
        [_get_org_auc_from_cache(s, 0) for s in ext_sites],
    ))

    if lgbm_all_rows:
        lgbm_df = pd.DataFrame(lgbm_all_rows)
        for variant, vlabel in [("lgbm_single", "LightGBM single-task"),
                                 ("lgbm_multi",  "LightGBM multi-task")]:
            vdf  = lgbm_df[lgbm_df["model"] == variant]
            a_df = vdf[vdf["site"] == site_a]
            ec_a = float(a_df[a_df["organism"].str.contains("coli", case=False)]["auc"].iloc[0]) \
                   if not a_df.empty else float("nan")
            mac  = float(a_df["auc"].mean()) if not a_df.empty else float("nan")
            ext  = []
            for s in ext_sites:
                sr = vdf[(vdf["site"] == s) & vdf["organism"].str.contains("coli", case=False)]
                ext.append(float(sr["auc"].iloc[0]) if not sr.empty else float("nan"))
            lines.append(_row(vlabel, ec_a, mac, ext))

    with open(out_dir / "paper" / "paper_table_multi_vs_single.md", "w") as f:
        f.writelines(lines)
    log.info("  Saved paper_table_multi_vs_single.md")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 7 — Artifact inventory
    # ═════════════════════════════════════════════════════════════════════════
    log.info(f"\n=== Task 7: Artifacts in {out_dir} ===")
    for fp in sorted(out_dir.rglob("*")):
        if fp.is_file():
            log.info(f"  {str(fp.relative_to(out_dir)):55s}  {fp.stat().st_size:>10,} B")

    # ═════════════════════════════════════════════════════════════════════════
    # TASK 9 — Summary
    # ═════════════════════════════════════════════════════════════════════════
    log.info("\n=== Task 9: Summary ===")
    log.info(f"Experiment : {args.experiment}")
    log.info(f"Pair profile: {args.pair_profile}")
    log.info(f"Temperature: {T_opt:.4f}")
    log.info(f"Thresholds : { {k: round(v,4) for k,v in base_thresholds.items()} }")
    log.info("\nPer-site macro AUC:")
    for r in macro_rows:
        log.info(f"  {r['site']:15s}  macro={r['macro_auc']:.3f}")
    log.info(f"\nOutput: {out_dir}")
    log.info("Done.")


# ═══════════════════════════════════════════════════════════════════════════════
# §19  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def _in_notebook():
    try:
        get_ipython()  # noqa: F821
        return True
    except NameError:
        return False


def main():
    if _in_notebook():
        # argparse can't read sys.argv in Jupyter/Colab/Kaggle — use NOTEBOOK_* config
        args = argparse.Namespace(
            pair_profile=NOTEBOOK_PAIR_PROFILE,
            mode=NOTEBOOK_MODE,
            experiment=NOTEBOOK_EXPERIMENT,
            ckpt_dir=NOTEBOOK_CKPT_DIR,
            calib_json=NOTEBOOK_CALIB_JSON,
            data_root=NOTEBOOK_DATA_ROOT,
            output_dir=NOTEBOOK_OUTPUT_DIR,
            with_ablation=NOTEBOOK_WITH_ABLATION,
            no_lgbm=NOTEBOOK_NO_LGBM,
            no_bn_adapt=NOTEBOOK_NO_BN_ADAPT,
            early_stop=NOTEBOOK_EARLY_STOP,
            seed_policy=NOTEBOOK_SEED_POLICY,
            top_k_seeds=NOTEBOOK_TOP_K_SEEDS,
            min_seed_macro_auc=NOTEBOOK_MIN_SEED_MACRO_AUC,
            prevalence_shift=NOTEBOOK_PREVALENCE_SHIFT,
            no_prevalence_shift=False,
            strict_mae_source_only=NOTEBOOK_STRICT_MAE_SOURCE_ONLY,
            no_saliency=NOTEBOOK_NO_SALIENCY,
            mae_epochs=NOTEBOOK_MAE_EPOCHS,
            with_random_cv=NOTEBOOK_WITH_RANDOM_CV,
            drug_conditioning=NOTEBOOK_DRUG_CONDITIONING,
            with_leave_one_drug_out=NOTEBOOK_WITH_LEAVE_ONE_DRUG_OUT,
        )
    else:
        p = argparse.ArgumentParser(description="Mega MALDI-TOF AMR pipeline")
        p.add_argument("--pair-profile", choices=sorted(PAIR_PROFILES), default=None,
                       help=f"Pair set to train/evaluate (default: {NOTEBOOK_PAIR_PROFILE})")
        p.add_argument("--mode", choices=["run14", "multidrug"], default=None,
                       help="Deprecated alias. run14 maps to pair-profile=run14; "
                            "multidrug maps to pair-profile=gram_negative6.")
        p.add_argument("--experiment", default=NOTEBOOK_EXPERIMENT,
                       help="Experiment name (used as output subdirectory)")
        p.add_argument("--ckpt-dir",   default=None,
                       help="Checkpoint directory (default: output_dir/experiment/models/)")
        p.add_argument("--calib-json", default=None,
                       help="Calibration JSON with saved T (and optionally per-pair thresholds)")
        p.add_argument("--data-root",  default=None,
                       help=f"DRIAMS data root (default: {DATA_ROOT})")
        p.add_argument("--output-dir", default="/kaggle/working/runs",
                       help="Output base directory")
        p.add_argument("--with-ablation", action="store_true",
                       help="Train and evaluate ablation variants (slow — trains new models)")
        p.add_argument("--no-lgbm",   action="store_true",
                       help="Skip LightGBM baselines")
        p.add_argument("--no-bn-adapt", action="store_true", default=NOTEBOOK_NO_BN_ADAPT,
                       help="Disable BN adaptation on external sites")
        p.add_argument("--bn-adapt", action="store_false", dest="no_bn_adapt",
                       help="Enable BN adaptation on external sites")
        p.add_argument("--early-stop", choices=["macro", "primary"], default=NOTEBOOK_EARLY_STOP,
                       help="Validation metric used for early stopping")
        p.add_argument("--seed-policy", choices=["all", "threshold", "topk"],
                       default=NOTEBOOK_SEED_POLICY,
                       help="How to select seeds for the final ensemble")
        p.add_argument("--top-k-seeds", type=int, default=TOP_K_SEEDS,
                       help="Number of seeds kept when --seed-policy topk")
        p.add_argument("--min-seed-macro-auc", type=float, default=MIN_SEED_MACRO_AUC,
                       help="Minimum source-val macro AUC for seed threshold/topk eligibility")
        p.add_argument("--prevalence-shift", choices=["capped", "none"],
                       default=NOTEBOOK_PREVALENCE_SHIFT,
                       help="External-site prevalence threshold adjustment mode")
        p.add_argument("--no-prevalence-shift", action="store_true",
                       help="Disable external-site prevalence threshold adjustment")
        p.add_argument("--strict-mae-source-only", action="store_true",
                       help="Pretrain MAE only on source train/val spectra")
        p.add_argument("--no-saliency", action="store_true",
                       default=NOTEBOOK_NO_SALIENCY,
                       help="Skip occlusion saliency interpretability export")
        p.add_argument("--mae-epochs", type=int, default=NOTEBOOK_MAE_EPOCHS,
                       help="MAE pretraining epochs; set 40 for the optional performance run")
        p.add_argument("--with-random-cv", action="store_true",
                       default=NOTEBOOK_WITH_RANDOM_CV,
                       help="Run LightGBM random-CV inflation diagnostic")
        p.add_argument("--no-random-cv", action="store_false",
                       dest="with_random_cv",
                       help="Skip LightGBM random-CV inflation diagnostic")
        p.add_argument("--drug-conditioning", choices=DRUG_CONDITIONING_CHOICES,
                       default=NOTEBOOK_DRUG_CONDITIONING,
                       help="Conditioning mode: task_id keeps legacy per-pair heads; "
                            "drug-aware modes use a shared head for drug panel experiments")
        p.add_argument("--with-leave-one-drug-out", action="store_true",
                       default=NOTEBOOK_WITH_LEAVE_ONE_DRUG_OUT,
                       help="Run opt-in leave-one-drug-out zero-shot diagnostic "
                            "(intended for --pair-profile ecoli_drug_panel)")
        args = p.parse_args()
        if args.pair_profile is None:
            args.pair_profile = "run14" if args.mode == "run14" else (
                "gram_negative6" if args.mode == "multidrug" else NOTEBOOK_PAIR_PROFILE)

    run_experiment(args)


if __name__ == "__main__":
    main()
