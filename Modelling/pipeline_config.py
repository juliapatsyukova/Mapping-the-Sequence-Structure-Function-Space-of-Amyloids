"""
pipeline_config.py
==================
Central configuration for the Amyloid-DB modelling pipeline.
"""

from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT  = Path(__file__).resolve().parent.parent          # v2/
DATA_DIR      = PROJECT_ROOT / "database"
OUTPUT_DIR    = PROJECT_ROOT / "output"
MODELING_DIR  = Path(__file__).resolve().parent                 # v2/modeling/

TSV_PATH      = DATA_DIR   / "consensus_unified.tsv"
METRICS_PATH  = OUTPUT_DIR / "sequence_metrics.csv"
PLOTS_DIR     = MODELING_DIR / "plots"
RESULTS_DIR   = MODELING_DIR / "results"
MODELS_DIR    = MODELING_DIR / "models"

# ── Feature columns ────────────────────────────────────────────────────────
# Continuous physicochemical descriptors
PHYS_CONTINUOUS = [
    "phys_length", "phys_molecular_weight",
    "phys_hydrophobicity_mean", "phys_hydrophobicity_std",
    "phys_hydrophobicity_max",  "phys_hydrophobicity_min",
    "phys_net_charge", "phys_positive_residues", "phys_negative_residues",
    "phys_charge_density",
    "phys_beta_propensity_mean", "phys_beta_propensity_max",
    "phys_aggregation_propensity_mean", "phys_aggregation_propensity_max",
    "phys_aromatic_fraction", "phys_aliphatic_fraction",
    "phys_polar_fraction", "phys_charged_fraction",
]

# Binary physicochemical indicators
PHYS_BINARY = [
    "phys_has_polyq", "phys_has_glycine_rich", "phys_has_proline",
]

# Amino acid composition (single-residue frequencies)
COMP_AA = [
    "comp_aa_A", "comp_aa_C", "comp_aa_D", "comp_aa_E", "comp_aa_F",
    "comp_aa_G", "comp_aa_H", "comp_aa_I", "comp_aa_K", "comp_aa_L",
    "comp_aa_M", "comp_aa_N", "comp_aa_P", "comp_aa_Q", "comp_aa_R",
    "comp_aa_S", "comp_aa_T", "comp_aa_V", "comp_aa_W", "comp_aa_Y",
]

# Dipeptide frequencies (only 3 non-trivially populated in this dataset)
COMP_DP = ["comp_dp_AA", "comp_dp_AQ", "comp_dp_QA"]

# Size-class fractions
COMP_SIZE = ["comp_tiny_fraction", "comp_small_fraction", "comp_large_fraction"]

# All feature columns used in modelling
ALL_FEATURES = PHYS_CONTINUOUS + PHYS_BINARY + COMP_AA + COMP_DP + COMP_SIZE

# Human-readable group labels for feature-importance plots
FEATURE_GROUPS = {
    "Physicochemical": PHYS_CONTINUOUS + PHYS_BINARY,
    "AA composition":  COMP_AA,
    "Dipeptide":       COMP_DP,
    "Size fractions":  COMP_SIZE,
}

# ── Target and weighting ───────────────────────────────────────────────────
TARGET      = "is_amyloid"          # binary: True / False
WEIGHT_COL  = "evidence_weight"     # 0.5 (literature) → 3.0 (structural)

# ── Split sizes ────────────────────────────────────────────────────────────
TEST_SIZE   = 0.15                  # held-out test set
VAL_SIZE    = 0.15                  # validation set (from remaining data)
N_CV_FOLDS  = 5
RANDOM_STATE = 42

# ── Evidence weight tiers (for confidence-stratified analysis) ─────────────
EVIDENCE_TIERS = {
    "Literature (0.5)":   (0.0,  0.6),
    "Spectroscopic (1.0)":(0.6,  1.5),
    "Kinetic (2.0)":      (1.5,  2.5),
    "Structural (3.0)":   (2.5, 10.0),
}

# ── Plot style ─────────────────────────────────────────────────────────────
PALETTE = {
    "amyloid":     "#E63946",   # red
    "non_amyloid": "#2A9D8F",   # teal
    "uncertain":   "#F4A261",   # amber
    "neutral":     "#457B9D",   # blue-grey
    "highlight":   "#F4A261",   # amber
}
FIGURE_DPI  = 300
FONT_SIZE   = 11
