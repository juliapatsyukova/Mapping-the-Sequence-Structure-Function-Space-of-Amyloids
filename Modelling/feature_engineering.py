"""
feature_engineering.py
=======================
Data loading, deduplication, and feature matrix preparation for
the Amyloid-DB modelling pipeline.
"""

import hashlib
import json
import warnings

import numpy as np
import pandas as pd

from pipeline_config import (
    TSV_PATH, ALL_FEATURES, PHYS_BINARY,
    TARGET, WEIGHT_COL, RANDOM_STATE
)

warnings.filterwarnings("ignore")


# ── Helpers ────────────────────────────────────────────────────────────────

def _safe(v) -> str:
    """Coerce to stripped string; empty string for None / NaN."""
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none", "") else s


def _sha8(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()[:16]


def _dedup_key(row: pd.Series) -> str:
    """Reproduce graph_build.py seq_node_id key (without the hash prefix)."""
    parts = [_safe(row.get("sequence")),
             _safe(row.get("uniprot_id")),
             _safe(row.get("region_start")),
             _safe(row.get("region_end"))]
    return "|".join(parts)


# ── Data loading ───────────────────────────────────────────────────────────

def load_raw_tsv() -> pd.DataFrame:
    """Load consensus_unified.tsv and apply basic type conversions."""
    df = pd.read_csv(TSV_PATH, sep="\t", dtype=str, keep_default_na=False)
    df.replace("nan", "", inplace=True)
    df.replace("None", "", inplace=True)

    # Boolean target
    df[TARGET] = df[TARGET].str.lower().isin(["true", "1", "yes"])

    # Binary feature columns: must be converted BEFORE numeric coercion
    # (they contain "True"/"False" strings that pd.to_numeric turns to NaN)
    for col in PHYS_BINARY:
        if col in df.columns:
            df[col] = df[col].astype(str).str.lower().isin(["true", "1", "yes"]).astype(float)

    # Continuous numeric columns (excludes PHYS_BINARY which are already handled)
    from pipeline_config import PHYS_CONTINUOUS, COMP_AA, COMP_DP, COMP_SIZE
    non_binary_features = PHYS_CONTINUOUS + COMP_AA + COMP_DP + COMP_SIZE
    num_cols = [WEIGHT_COL, "confidence"] + non_binary_features
    for col in num_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Composition dipeptide columns: empty string means the dipeptide was absent → 0
    for col in COMP_DP:
        if col in df.columns:
            df[col] = df[col].fillna(0.0)

    # Default evidence weight when missing
    df[WEIGHT_COL] = df[WEIGHT_COL].fillna(0.5)

    return df


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse to one row per unique biological identity
    (sequence, uniprot_id, region_start, region_end).

    Priority rule:
      1. Prefer rows with pre-computed features (non-empty phys_length).
      2. Among equally-featured rows, keep the one with highest evidence_weight.

    After deduplication, the max evidence_weight across all rows in the group
    is stored in `evidence_weight_max` for use as a sample weight — this
    captures the best experimental evidence for any observation of that sequence,
    even if that observation came from a no-sequence entry.

    This reproduces the ~4,348 unique SequenceInstance nodes in the graph.
    """
    df = df.copy()
    df["_dedup_key"]    = df.apply(_dedup_key, axis=1)
    df["_has_features"] = df["phys_length"].astype(str).str.strip().ne("").astype(int)

    # Record max evidence_weight per group before collapsing
    max_ew = (df.groupby("_dedup_key")[WEIGHT_COL]
                .max()
                .rename("evidence_weight_max")
                .reset_index())

    # Sort: feature-complete rows first, then highest evidence_weight
    df = (
        df.sort_values(["_has_features", WEIGHT_COL], ascending=[False, False])
          .drop_duplicates(subset=["_dedup_key"], keep="first")
    )

    # Attach max evidence weight
    df = df.merge(max_ew, on="_dedup_key", how="left")
    df = df.drop(columns=["_dedup_key", "_has_features"])
    df = df.reset_index(drop=True)

    return df


# ── Feature matrix ─────────────────────────────────────────────────────────

def build_feature_matrix(df: pd.DataFrame):
    """
    Build X (feature matrix), y (binary target), w (sample weights),
    and a metadata frame for downstream analysis.

    Returns
    -------
    X_full : pd.DataFrame — all rows, including those with missing features
    y      : pd.Series   — binary target (int)
    w      : pd.Series   — evidence-based sample weights
    meta   : pd.DataFrame — non-feature columns for analysis
    """
    present = [c for c in ALL_FEATURES if c in df.columns]
    missing_cols = [c for c in ALL_FEATURES if c not in df.columns]
    if missing_cols:
        print(f"  [WARN] Feature columns absent from TSV (will be NaN): {missing_cols}")

    X = df[present].copy()
    # Add missing columns as NaN so downstream code always sees expected shape
    for c in missing_cols:
        X[c] = np.nan
    X = X[ALL_FEATURES]  # enforce canonical order

    # Dipeptide absence should be 0, not NaN
    from pipeline_config import COMP_DP
    for c in COMP_DP:
        if c in X.columns:
            X[c] = X[c].fillna(0.0)

    y = df[TARGET].astype(int)
    w = df[WEIGHT_COL].fillna(0.5)

    meta_cols = [c for c in df.columns if c not in ALL_FEATURES]
    meta = df[meta_cols].copy()

    return X, y, w, meta


def load_and_prepare():
    """
    Full preparation pipeline:
    load → deduplicate → split into with-sequence and without-sequence subsets.

    Returns
    -------
    dict with keys:
        'df_full'   — deduplicated DataFrame (4 348 rows)
        'df_seq'    — rows that have sequence text and valid features (~3 850)
        'X', 'y', 'w', 'meta'  — aligned feature matrix / target / weights / metadata
                                  for df_seq only (the modelling subset)
        'feature_names'         — list of feature names
    """
    print("Loading consensus_unified.tsv …")
    df_raw = load_raw_tsv()
    print(f"  Loaded {len(df_raw):,} rows")

    print("Deduplicating …")
    df = deduplicate(df_raw)
    print(f"  {len(df):,} unique SequenceInstances after deduplication")
    print(f"  Class distribution: amyloid={df[TARGET].sum():,}  "
          f"non-amyloid={(~df[TARGET].astype(bool)).sum():,}")

    # Subset with sequence text (for reliable feature computation)
    has_seq = df["sequence"].str.strip().ne("")
    df_seq = df[has_seq].copy()
    print(f"  With sequence text: {len(df_seq):,}  "
          f"(without: {(~has_seq).sum():,} — excluded from ML)")

    # Use max evidence weight across all observations for sample weighting
    # (this captures e.g. structural evidence even when the structural row
    #  lacked sequence text and was superseded by a feature-complete row)
    if "evidence_weight_max" in df_seq.columns:
        df_seq = df_seq.copy()
        df_seq[WEIGHT_COL] = df_seq["evidence_weight_max"]

    X, y, w, meta = build_feature_matrix(df_seq)

    # Keep only rows where ALL features are present (not NaN)
    valid_mask = X.notna().all(axis=1)
    X   = X[valid_mask]
    y   = y[valid_mask]
    w   = w[valid_mask]
    meta = meta[valid_mask]

    print(f"  Rows with complete features: {valid_mask.sum():,}")

    # Fix 1: exclude literature-tier entries (evidence_weight < 1.0).
    # These are AmyPro (169) and WALTZ-DB literature rows (227) whose labels
    # were assigned by computational tools (WALTZ/TANGO) rather than assay.
    # All remaining entries have experimental backing: ThT (CPAD, 1.0),
    # kinetic seeding (AmyloGraph, 2.0), or structural methods (3.0).
    exp_mask = w >= 1.0
    n_dropped = (~exp_mask).sum()
    X    = X[exp_mask]
    y    = y[exp_mask]
    w    = w[exp_mask]
    meta = meta[exp_mask]
    print(f"  Dropped {n_dropped:,} literature-tier rows (evidence_weight < 1.0)")

    print(f"  Final class distribution: amyloid={y.sum():,}  "
          f"non-amyloid={(1-y).sum():,}")

    return {
        "df_full": df,
        "df_seq":  df_seq,
        "X":       X,
        "y":       y,
        "w":       w,
        "meta":    meta,
        "feature_names": list(X.columns),
    }


# ── Feature name mapping (for readable plots) ──────────────────────────────

FEATURE_LABELS = {
    "phys_length":                     "Sequence length",
    "phys_molecular_weight":           "Molecular weight",
    "phys_hydrophobicity_mean":        "Hydrophobicity (mean)",
    "phys_hydrophobicity_std":         "Hydrophobicity (std)",
    "phys_hydrophobicity_max":         "Hydrophobicity (max)",
    "phys_hydrophobicity_min":         "Hydrophobicity (min)",
    "phys_net_charge":                 "Net charge",
    "phys_positive_residues":          "Positive residues",
    "phys_negative_residues":          "Negative residues",
    "phys_charge_density":             "Charge density",
    "phys_beta_propensity_mean":       "β-sheet propensity (mean)",
    "phys_beta_propensity_max":        "β-sheet propensity (max)",
    "phys_aggregation_propensity_mean":"Aggregation propensity (mean)",
    "phys_aggregation_propensity_max": "Aggregation propensity (max)",
    "phys_aromatic_fraction":          "Aromatic fraction",
    "phys_aliphatic_fraction":         "Aliphatic fraction",
    "phys_polar_fraction":             "Polar fraction",
    "phys_charged_fraction":           "Charged fraction",
    "phys_has_polyq":                  "Has polyQ run",
    "phys_has_glycine_rich":           "Glycine-rich",
    "phys_has_proline":                "Has proline",
    "comp_aa_A":  "Ala (A)",  "comp_aa_C": "Cys (C)",  "comp_aa_D": "Asp (D)",
    "comp_aa_E":  "Glu (E)",  "comp_aa_F": "Phe (F)",  "comp_aa_G": "Gly (G)",
    "comp_aa_H":  "His (H)",  "comp_aa_I": "Ile (I)",  "comp_aa_K": "Lys (K)",
    "comp_aa_L":  "Leu (L)",  "comp_aa_M": "Met (M)",  "comp_aa_N": "Asn (N)",
    "comp_aa_P":  "Pro (P)",  "comp_aa_Q": "Gln (Q)",  "comp_aa_R": "Arg (R)",
    "comp_aa_S":  "Ser (S)",  "comp_aa_T": "Thr (T)",  "comp_aa_V": "Val (V)",
    "comp_aa_W":  "Trp (W)",  "comp_aa_Y": "Tyr (Y)",
    "comp_dp_AA": "Dipeptide AA",
    "comp_dp_AQ": "Dipeptide AQ",
    "comp_dp_QA": "Dipeptide QA",
    "comp_tiny_fraction":  "Tiny residue fraction",
    "comp_small_fraction": "Small residue fraction",
    "comp_large_fraction": "Large residue fraction",
}


def readable_name(col: str) -> str:
    return FEATURE_LABELS.get(col, col)
