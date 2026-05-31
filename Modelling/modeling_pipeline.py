#!/usr/bin/env python3
"""
modeling_pipeline.py
====================
Complete Amyloid-DB modelling pipeline.

Sections
--------
 1. Imports and setup
 2. Data loading
 3. Exploratory data analysis  → plots 01–05
 4. Feature analysis           → plots 06–07
 5. Train / validation / test split
 6. Model training (LR, RF, GB) with evidence-weighted sample weights
 7. Cross-validation
 8. Test-set evaluation        → plots 08–12
 9. Feature importance         → plots 13–14
10. Confidence-stratified analysis → plots 15–16
11. Export (models, predictions, tables, captions)

Run
---
    cd /Users/julia_patsiukova/Documents/GitHub/Amyloid-DB/v2/modeling
    python modeling_pipeline.py
"""

# ═══════════════════════════════════════════════════════════════════════════
# 1. IMPORTS AND SETUP
# ═══════════════════════════════════════════════════════════════════════════

import json
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for script execution
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    RocCurveDisplay, PrecisionRecallDisplay,
    auc, average_precision_score, balanced_accuracy_score,
    classification_report, confusion_matrix, f1_score,
    matthews_corrcoef, precision_recall_curve, roc_auc_score, roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── project imports ──────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from pipeline_config import (
    PLOTS_DIR, RESULTS_DIR, MODELS_DIR,
    ALL_FEATURES, PHYS_CONTINUOUS, PHYS_BINARY,
    COMP_AA, COMP_DP, COMP_SIZE, FEATURE_GROUPS,
    TARGET, WEIGHT_COL,
    TEST_SIZE, VAL_SIZE, N_CV_FOLDS, RANDOM_STATE,
    EVIDENCE_TIERS, PALETTE, FIGURE_DPI, FONT_SIZE,
)
from feature_engineering import load_and_prepare, readable_name, FEATURE_LABELS

np.random.seed(RANDOM_STATE)

# ── Global plot style ────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size":        FONT_SIZE,
    "axes.titlesize":   FONT_SIZE + 1,
    "axes.labelsize":   FONT_SIZE,
    "legend.fontsize":  FONT_SIZE - 1,
    "xtick.labelsize":  FONT_SIZE - 1,
    "ytick.labelsize":  FONT_SIZE - 1,
    "figure.dpi":       100,           # screen; save at 300
    "savefig.dpi":      FIGURE_DPI,
    "savefig.bbox":     "tight",
    "axes.spines.top":  False,
    "axes.spines.right":False,
})
sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)

COLOR_AMY  = PALETTE["amyloid"]     # red
COLOR_NON  = PALETTE["non_amyloid"] # teal
COLOR_UNC  = PALETTE["uncertain"]   # amber
COLOR_NEU  = PALETTE["neutral"]     # blue-grey

MODEL_COLORS = {
    "Logistic Regression":  "#457B9D",
    "Random Forest":        "#E63946",
    "Gradient Boosting":    "#2A9D8F",
}

CAPTIONS: dict[str, str] = {}   # filled in as plots are saved


def _savefig(fname: str, caption: str):
    path = PLOTS_DIR / fname
    plt.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close()
    CAPTIONS[fname] = caption
    print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# 2. DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════

print("=" * 70)
print("Amyloid-DB Modelling Pipeline")
print("=" * 70)

data = load_and_prepare()
df_full = data["df_full"]
df_seq  = data["df_seq"]
X_all   = data["X"]
y_all   = data["y"]
w_all   = data["w"]
meta    = data["meta"]
FEAT    = data["feature_names"]

# Use the filtered modelling subset for EDA so all plots are consistent
# with the training data (3,594 experimentally-supported instances only).
df_eda  = meta.reset_index(drop=True).copy()

# ═══════════════════════════════════════════════════════════════════════════
# 3. EXPLORATORY DATA ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
print("\n[EDA] Generating dataset overview plots …")

# ── 01 Class distribution ────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

counts = df_eda[TARGET].value_counts()
labels = ["Amyloid\n(positive)", "Non-amyloid\n(negative)"]
colors = [COLOR_AMY, COLOR_NON]
bars = axes[0].bar(labels,
                   [counts.get(True, 0), counts.get(False, 0)],
                   color=colors, width=0.5, edgecolor="white", linewidth=1.5)
axes[0].set_title("Class distribution")
axes[0].set_ylabel("Number of sequence instances")
for bar in bars:
    axes[0].text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 25,
                 f"{bar.get_height():,}",
                 ha="center", va="bottom", fontsize=FONT_SIZE, fontweight="bold")
axes[0].set_ylim(0, max(counts.values) * 1.18)
axes[0].grid(axis="y", alpha=0.4)

# Evidence weight by class
ev_amy  = df_eda.loc[df_eda[TARGET],     WEIGHT_COL].dropna()
ev_non  = df_eda.loc[~df_eda[TARGET],    WEIGHT_COL].dropna()
bins = np.arange(0.25, 3.76, 0.5)
axes[1].hist(ev_amy,  bins=bins, alpha=0.75, color=COLOR_AMY,  label="Amyloid",
             density=True, edgecolor="white")
axes[1].hist(ev_non,  bins=bins, alpha=0.75, color=COLOR_NON,  label="Non-amyloid",
             density=True, edgecolor="white")
axes[1].set_xlabel("Evidence weight")
axes[1].set_ylabel("Density")
axes[1].set_title("Evidence weight by class")
axes[1].legend()
axes[1].set_xticks([0.5, 1.0, 2.0, 3.0])
axes[1].set_xticklabels(["Lit.\n(0.5)", "Spectro.\n(1.0)", "Kinetic\n(2.0)", "Struct.\n(3.0)"])

plt.tight_layout()
_savefig("01_class_distribution.png",
         "Fig. 1. Class distribution of the Amyloid-DB consensus dataset. "
         "(Left) Total counts of amyloid-positive and non-amyloid SequenceInstances "
         "after database integration and deduplication (n = 4,348). "
         "(Right) Distribution of evidence weights by class, where higher weights "
         "reflect stronger experimental evidence (structural methods weight 3.0 vs. "
         "literature-curated 0.5).")

# ── 02 Sequence length distribution ─────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

lengths_amy = df_eda.loc[df_eda[TARGET],  "sequence"].str.len().dropna()
lengths_non = df_eda.loc[~df_eda[TARGET], "sequence"].str.len().dropna()
lengths_all = df_eda["sequence"].str.len().dropna()

# ≤60 aa panel
mask_short = lengths_all <= 60
axes[0].hist(lengths_amy[lengths_amy <= 60], bins=30, alpha=0.75,
             color=COLOR_AMY,  label="Amyloid",     density=True, edgecolor="white")
axes[0].hist(lengths_non[lengths_non <= 60], bins=30, alpha=0.75,
             color=COLOR_NON,  label="Non-amyloid", density=True, edgecolor="white")
axes[0].set_xlabel("Sequence length (aa)")
axes[0].set_ylabel("Density")
axes[0].set_title("Sequence length ≤ 60 aa")
axes[0].legend()

# Full range (log x)
axes[1].hist(lengths_amy, bins=50, alpha=0.75,
             color=COLOR_AMY,  label="Amyloid",     density=True, edgecolor="white")
axes[1].hist(lengths_non, bins=50, alpha=0.75,
             color=COLOR_NON,  label="Non-amyloid", density=True, edgecolor="white")
axes[1].set_xscale("log")
axes[1].set_xlabel("Sequence length (aa, log scale)")
axes[1].set_ylabel("Density")
axes[1].set_title("Full length range")
axes[1].legend()

plt.tight_layout()
_savefig("02_sequence_length_distribution.png",
         "Fig. 2. Sequence length distribution by amyloid class. "
         "(Left) Short sequences (≤60 aa), which include the large hexapeptide library "
         "from WALTZ-DB and CPAD. (Right) Full range on a log scale. "
         "The distribution is strongly bimodal: most entries are short synthetic "
         "peptides (median 6 aa), with a second peak at full-protein lengths (100–2700 aa).")

# ── 03 Source database composition ──────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

src_counts = (df_eda.groupby(["source_db", TARGET])
                    .size()
                    .reset_index(name="n"))
src_totals = src_counts.groupby("source_db")["n"].sum().sort_values(ascending=False)
order = src_totals.index.tolist()

src_pivot = src_counts.pivot(index="source_db", columns=TARGET, values="n").fillna(0)
src_pivot = src_pivot.loc[order]

src_pivot.plot(kind="bar", ax=axes[0],
               color=[COLOR_NON, COLOR_AMY], edgecolor="white", linewidth=0.8)
axes[0].set_xlabel("Source database")
axes[0].set_ylabel("Number of sequences")
axes[0].set_title("Sequences contributed by source database")
axes[0].set_xticklabels(order, rotation=35, ha="right")
handles = [mpatches.Patch(color=COLOR_NON, label="Non-amyloid"),
           mpatches.Patch(color=COLOR_AMY, label="Amyloid")]
axes[0].legend(handles=handles)

# Evidence type composition
ev_counts = df_eda["evidence_type"].value_counts()
ev_labels = ev_counts.index.tolist()
ev_colors = sns.color_palette("Set2", len(ev_labels))
axes[1].barh(ev_labels, ev_counts.values, color=ev_colors, edgecolor="white")
axes[1].set_xlabel("Number of sequences")
axes[1].set_title("Evidence type composition")
for i, v in enumerate(ev_counts.values):
    axes[1].text(v + 15, i, f"{v:,}", va="center", fontsize=FONT_SIZE - 1)
axes[1].set_xlim(0, ev_counts.max() * 1.15)

plt.tight_layout()
_savefig("03_source_database_composition.png",
         "Fig. 3. Composition of the integrated Amyloid-DB dataset by source. "
         "(Left) Amyloid-positive and non-amyloid sequences contributed by each of "
         "the eight integrated databases. WALTZ-DB and CPAD dominate through their "
         "hexapeptide scan libraries; AmyloidAtlas and Cross-Beta DB contribute "
         "structurally-validated entries. (Right) Distribution of evidence type "
         "across all sequences.")

# ── 04 Evidence method composition ──────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

meth_counts = df_eda["method_universal"].replace("", "Unknown").value_counts()
ev_w_counts = df_eda[WEIGHT_COL].value_counts().sort_index()

colors_meth = sns.color_palette("tab10", len(meth_counts))
axes[0].barh(meth_counts.index.tolist(), meth_counts.values,
             color=colors_meth, edgecolor="white")
axes[0].set_xlabel("Number of sequences")
axes[0].set_title("Normalised experimental method")
for i, v in enumerate(meth_counts.values):
    axes[0].text(v + 10, i, f"{v:,}", va="center", fontsize=FONT_SIZE - 2)
axes[0].set_xlim(0, meth_counts.max() * 1.18)

tier_labels = list(EVIDENCE_TIERS.keys())
tier_counts = []
for lbl, (lo, hi) in EVIDENCE_TIERS.items():
    mask = (df_eda[WEIGHT_COL] > lo) & (df_eda[WEIGHT_COL] <= hi)
    tier_counts.append(mask.sum())

bar_colors = ["#D4E6F1", "#AED6F1", "#5DADE2", "#1A5276"]
bars = axes[1].bar(tier_labels, tier_counts, color=bar_colors, edgecolor="white")
axes[1].set_ylabel("Number of sequences")
axes[1].set_title("Evidence strength tier")
axes[1].set_xticklabels(tier_labels, rotation=20, ha="right")
for bar in bars:
    axes[1].text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 15,
                 f"{bar.get_height():,}", ha="center", va="bottom", fontsize=FONT_SIZE - 1)
axes[1].set_ylim(0, max(tier_counts) * 1.18)

plt.tight_layout()
_savefig("04_evidence_method_composition.png",
         "Fig. 4. Distribution of experimental evidence across the dataset. "
         "(Left) Counts of sequences by normalised experimental method. ThT binding "
         "and electron microscopy dominate, reflecting the method preferences of the "
         "contributing databases. (Right) Evidence strength tier distribution, "
         "stratified into literature-curated (weight 0.5), spectroscopic (1.0), "
         "kinetic (2.0), and structural (3.0) observations.")

# ── 05 Confidence distribution ──────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

conf_amy = pd.to_numeric(df_eda.loc[df_eda[TARGET],  "confidence"], errors="coerce").dropna()
conf_non = pd.to_numeric(df_eda.loc[~df_eda[TARGET], "confidence"], errors="coerce").dropna()

axes[0].hist(conf_amy, bins=20, alpha=0.75, color=COLOR_AMY,  label="Amyloid",
             density=True, edgecolor="white")
axes[0].hist(conf_non, bins=20, alpha=0.75, color=COLOR_NON,  label="Non-amyloid",
             density=True, edgecolor="white")
axes[0].set_xlabel("Confidence score")
axes[0].set_ylabel("Density")
axes[0].set_title("Confidence score by class")
axes[0].legend()

ew_vals = df_eda[WEIGHT_COL].dropna()
axes[1].hist(ew_vals, bins=15, color=COLOR_NEU, edgecolor="white", alpha=0.8)
axes[1].set_xlabel("Evidence weight")
axes[1].set_ylabel("Count")
axes[1].set_title("Evidence weight distribution")
axes[1].axvline(ew_vals.mean(), color="black", linestyle="--",
                label=f"Mean = {ew_vals.mean():.2f}")
axes[1].legend()

plt.tight_layout()
_savefig("05_confidence_distribution.png",
         "Fig. 5. Distribution of confidence scores and evidence weights. "
         "(Left) Confidence scores (0–100) by amyloid class; amyloid-positive entries "
         "show higher average confidence reflecting preferential coverage by structural "
         "methods. (Right) Evidence weight distribution across all 4,348 instances; "
         "the dominant peak at 1.0 corresponds to spectroscopic/binding assays.")

# ═══════════════════════════════════════════════════════════════════════════
# 4. FEATURE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Features] Generating feature analysis plots …")

# ── 06 Key physicochemical features by class ─────────────────────────────
key_phys = [
    "phys_hydrophobicity_mean", "phys_beta_propensity_mean",
    "phys_aggregation_propensity_mean", "phys_aromatic_fraction",
    "phys_charged_fraction", "phys_length",
]
key_labels = [readable_name(f) for f in key_phys]

df_plot = X_all[key_phys].copy()
df_plot["Class"] = y_all.map({1: "Amyloid", 0: "Non-amyloid"})

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()

for i, (col, lbl) in enumerate(zip(key_phys, key_labels)):
    ax = axes[i]
    data_amy = df_plot.loc[df_plot["Class"] == "Amyloid",     col].dropna()
    data_non = df_plot.loc[df_plot["Class"] == "Non-amyloid", col].dropna()
    bp = ax.boxplot([data_non.values, data_amy.values],
               labels=["Non-amyloid", "Amyloid"],
               patch_artist=True,
               boxprops=dict(facecolor=COLOR_NON + "55"),
               medianprops=dict(color="black", linewidth=2),
               flierprops=dict(marker=".", markersize=2, alpha=0.3))
    if len(bp["boxes"]) >= 2:
        bp["boxes"][0].set_facecolor(COLOR_NON + "55")
        bp["boxes"][1].set_facecolor(COLOR_AMY + "55")
    ax.set_title(lbl, fontsize=FONT_SIZE)
    ax.grid(axis="y", alpha=0.3)

plt.suptitle("Key physicochemical features by amyloid class", fontsize=FONT_SIZE + 2, y=1.01)
plt.tight_layout()
_savefig("06_physicochemical_by_class.png",
         "Fig. 6. Distribution of key physicochemical sequence descriptors by amyloid "
         "classification. Box plots show median, interquartile range, and outliers. "
         "Amyloid-positive sequences exhibit systematically higher hydrophobicity, "
         "β-sheet propensity, and aggregation propensity, while non-amyloid sequences "
         "show higher charged fractions.")

# ── 07 Feature correlation heatmap (physicochemical only) ────────────────
phys_cols_corr = PHYS_CONTINUOUS
corr_df = X_all[phys_cols_corr].corr()
labels_corr = [readable_name(c) for c in phys_cols_corr]

fig, ax = plt.subplots(figsize=(11, 9))
mask = np.triu(np.ones_like(corr_df, dtype=bool), k=1)
sns.heatmap(corr_df, mask=mask, annot=True, fmt=".2f", cmap="coolwarm",
            center=0, linewidths=0.5, annot_kws={"size": 7},
            xticklabels=labels_corr, yticklabels=labels_corr, ax=ax,
            cbar_kws={"shrink": 0.8, "label": "Pearson r"})
ax.set_title("Physicochemical feature correlation matrix")
plt.xticks(rotation=35, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()
_savefig("07_feature_correlation.png",
         "Fig. 7. Pearson correlation matrix of physicochemical sequence descriptors. "
         "Strong positive correlations exist between sequence length, molecular weight, "
         "and residue counts (as expected). Hydrophobicity and charged fraction are "
         "negatively correlated, reflecting the compositional trade-off between "
         "hydrophobic and ionic residues.")

# ═══════════════════════════════════════════════════════════════════════════
# 5. TRAIN / VALIDATION / TEST SPLIT
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Split] Preparing train / val / test sets …")

X_temp, X_test, y_temp, y_test, w_temp, w_test, meta_temp, meta_test = train_test_split(
    X_all, y_all, w_all, meta,
    test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y_all
)
# val_size relative to remaining
val_relative = VAL_SIZE / (1.0 - TEST_SIZE)
X_train, X_val, y_train, y_val, w_train, w_val, meta_train, meta_val = train_test_split(
    X_temp, y_temp, w_temp, meta_temp,
    test_size=val_relative, random_state=RANDOM_STATE, stratify=y_temp
)

print(f"  Train : {len(X_train):,} | Val : {len(X_val):,} | Test : {len(X_test):,}")
print(f"  Train amyloid: {y_train.sum():,} / {len(y_train):,} "
      f"({100*y_train.mean():.1f}%)")

# ═══════════════════════════════════════════════════════════════════════════
# 6. MODEL TRAINING
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Training] Fitting models …")

# Imputer + scaler shared across pipelines
imputer  = SimpleImputer(strategy="median")
scaler   = StandardScaler()

# Build pipelines
MODELS = {
    "Logistic Regression": Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale",  StandardScaler()),
        ("clf",    LogisticRegression(C=0.5, max_iter=2000,
                                      solver="saga",   # saga avoids MKL/BLAS crash
                                      random_state=RANDOM_STATE,
                                      class_weight="balanced")),
    ]),
    "Random Forest": Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("clf",    RandomForestClassifier(
                       n_estimators=500, min_samples_leaf=5,
                       max_features="sqrt", n_jobs=1,
                       random_state=RANDOM_STATE)),
    ]),
    "Gradient Boosting": Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("clf",    GradientBoostingClassifier(
                       n_estimators=300, learning_rate=0.05,
                       max_depth=5, subsample=0.8,
                       random_state=RANDOM_STATE)),
    ]),
}

# Fit each model with evidence-based sample weights
fitted_models: dict[str, Pipeline] = {}
for name, pipe in MODELS.items():
    print(f"  Fitting {name} …", end=" ", flush=True)
    pipe.fit(X_train, y_train, clf__sample_weight=w_train)
    fitted_models[name] = pipe
    print("done")

# ═══════════════════════════════════════════════════════════════════════════
# 7. CROSS-VALIDATION
# ═══════════════════════════════════════════════════════════════════════════
print("\n[CV] Running 5-fold stratified cross-validation on training set …")

cv = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)
cv_results: dict[str, dict] = {}

for name, pipe in MODELS.items():
    scores = cross_validate(
        pipe, X_train, y_train,
        cv=cv,
        scoring={"roc_auc": "roc_auc", "avg_precision": "average_precision",
                 "f1": "f1", "bal_acc": "balanced_accuracy"},
        fit_params={"clf__sample_weight": w_train},
        n_jobs=1,    # n_jobs=-1 can trigger MKL issues on this macOS environment
    )
    cv_results[name] = scores
    print(f"  {name:25s}  "
          f"ROC-AUC={scores['test_roc_auc'].mean():.3f}±{scores['test_roc_auc'].std():.3f}  "
          f"PR-AUC={scores['test_avg_precision'].mean():.3f}±{scores['test_avg_precision'].std():.3f}")

# ── CV performance plot ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

metrics_to_plot = [("test_roc_auc", "ROC-AUC"), ("test_avg_precision", "PR-AUC")]
for ax, (metric_key, metric_label) in zip(axes, metrics_to_plot):
    data_for_box = [cv_results[n][metric_key] for n in MODELS.keys()]
    colors_box   = [MODEL_COLORS[n] for n in MODELS.keys()]
    bp = ax.boxplot(data_for_box, patch_artist=True, notch=False,
                    labels=list(MODELS.keys()),
                    medianprops=dict(color="white", linewidth=2.5))
    for patch, col in zip(bp["boxes"], colors_box):
        patch.set_facecolor(col)
        patch.set_alpha(0.8)
    ax.set_title(f"5-fold CV {metric_label}")
    ax.set_ylabel(metric_label)
    ax.set_ylim(max(0, min(d.min() for d in data_for_box) - 0.05), 1.02)
    ax.set_xticklabels(list(MODELS.keys()), rotation=15, ha="right")
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, label="random")
    ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
_savefig("08_cv_performance.png",
         "Fig. 8. Five-fold cross-validation performance on the training set. "
         "(Left) ROC-AUC distributions across folds; (Right) PR-AUC distributions. "
         "Random Forest and Gradient Boosting substantially outperform the linear "
         "baseline. Narrow boxplots indicate stable performance across folds.")

# ═══════════════════════════════════════════════════════════════════════════
# 7b. FIX 2 — ISOTONIC CALIBRATION OF GRADIENT BOOSTING ON VALIDATION SET
# ═══════════════════════════════════════════════════════════════════════════
# CalibratedClassifierCV(cv="prefit") wraps the already-fitted Pipeline and
# learns an isotonic mapping from raw predict_proba output → observed positive
# rate on the held-out validation set. X_val was never seen during training or
# CV, so this does not introduce leakage.
print("\n[Calibration] Fitting isotonic calibration on validation set …")

gb_calibrated = CalibratedClassifierCV(
    fitted_models["Gradient Boosting"], method="isotonic", cv="prefit"
)
gb_calibrated.fit(X_val, y_val)
fitted_models["GB (Calibrated)"] = gb_calibrated
MODEL_COLORS["GB (Calibrated)"] = "#264653"  # dark slate to distinguish from GB

# Quick sanity check: calibration should not change ROC-AUC appreciably
_cal_auc = roc_auc_score(y_val, gb_calibrated.predict_proba(X_val)[:, 1])
_raw_auc = roc_auc_score(y_val, fitted_models["Gradient Boosting"].predict_proba(X_val)[:, 1])
print(f"  Val ROC-AUC  raw GB: {_raw_auc:.4f}  |  calibrated GB: {_cal_auc:.4f}")
print("  (AUC should be unchanged; only probability shape improves)")

# ═══════════════════════════════════════════════════════════════════════════
# 8. TEST-SET EVALUATION
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Evaluation] Evaluating on held-out test set …")

test_metrics: dict[str, dict] = {}

for name, pipe in fitted_models.items():
    y_prob = pipe.predict_proba(X_test)[:, 1]
    y_pred = pipe.predict(X_test)
    test_metrics[name] = {
        "roc_auc":    roc_auc_score(y_test, y_prob),
        "avg_prec":   average_precision_score(y_test, y_prob),
        "f1":         f1_score(y_test, y_pred),
        "bal_acc":    balanced_accuracy_score(y_test, y_pred),
        "mcc":        matthews_corrcoef(y_test, y_pred),
        "y_prob":     y_prob,
        "y_pred":     y_pred,
    }
    print(f"  {name:25s}  "
          f"ROC-AUC={test_metrics[name]['roc_auc']:.4f}  "
          f"PR-AUC={test_metrics[name]['avg_prec']:.4f}  "
          f"F1={test_metrics[name]['f1']:.4f}  "
          f"MCC={test_metrics[name]['mcc']:.4f}")

# Best model (highest ROC-AUC)
best_name = max(test_metrics, key=lambda n: test_metrics[n]["roc_auc"])
print(f"\n  Best model: {best_name} (ROC-AUC = {test_metrics[best_name]['roc_auc']:.4f})")

# ── 09 ROC curves ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5))
for name, pipe in fitted_models.items():
    y_prob = test_metrics[name]["y_prob"]
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    auc_val = test_metrics[name]["roc_auc"]
    ax.plot(fpr, tpr, lw=2, color=MODEL_COLORS[name],
            label=f"{name} (AUC = {auc_val:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1, label="Random (0.5)")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("ROC curves — test set")
ax.legend(loc="lower right")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
ax.grid(alpha=0.3)
plt.tight_layout()
_savefig("09_roc_curves.png",
         "Fig. 9. Receiver Operating Characteristic (ROC) curves on the held-out test "
         "set for all three models. Area under the curve (AUC) values are reported "
         "in the legend. Random Forest and Gradient Boosting achieve substantially "
         "higher AUC than Logistic Regression, confirming the non-linear relationship "
         "between sequence features and amyloidogenicity.")

# ── 10 PR curves ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5))
baseline_prec = y_test.mean()
ax.axhline(baseline_prec, color="grey", linestyle="--", lw=1,
           label=f"Baseline ({baseline_prec:.2f})")
for name, pipe in fitted_models.items():
    y_prob = test_metrics[name]["y_prob"]
    prec, rec, _ = precision_recall_curve(y_test, y_prob)
    ap = test_metrics[name]["avg_prec"]
    ax.plot(rec, prec, lw=2, color=MODEL_COLORS[name],
            label=f"{name} (AP = {ap:.3f})")
ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
ax.set_title("Precision-Recall curves — test set")
ax.legend(loc="lower left")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.01)
ax.grid(alpha=0.3)
plt.tight_layout()
_savefig("10_pr_curves.png",
         "Fig. 10. Precision-Recall curves on the held-out test set. Average Precision "
         "(AP) values are reported in the legend. The dashed horizontal line represents "
         "the no-skill baseline (class prevalence). Gradient Boosting and Random Forest "
         "maintain high precision across a wide range of recall values.")

# ── 11 Confusion matrix (best model) ────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4.5))
cm = confusion_matrix(y_test, test_metrics[best_name]["y_pred"])
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", ax=ax,
            xticklabels=["Non-amyloid", "Amyloid"],
            yticklabels=["Non-amyloid", "Amyloid"],
            annot_kws={"size": 14, "weight": "bold"})
ax.set_xlabel("Predicted label")
ax.set_ylabel("True label")
ax.set_title(f"Confusion matrix — {best_name}")
plt.tight_layout()
_savefig("11_confusion_matrix.png",
         f"Fig. 11. Confusion matrix for {best_name} on the held-out test set "
         f"(n = {len(y_test)}). Rows indicate true labels; columns indicate predicted "
         "labels. The model achieves good discrimination in both classes.")

# ── 12 Calibration plot ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 5))
ax.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
for name, pipe in fitted_models.items():
    y_prob = test_metrics[name]["y_prob"]
    frac_pos, mean_pred = calibration_curve(y_test, y_prob, n_bins=10)
    ax.plot(mean_pred, frac_pos, "o-", lw=2, color=MODEL_COLORS[name], label=name)
ax.set_xlabel("Mean predicted probability")
ax.set_ylabel("Fraction positive")
ax.set_title("Calibration curves — test set")
ax.legend()
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.grid(alpha=0.3)
plt.tight_layout()
_savefig("12_calibration_plot.png",
         "Fig. 12. Calibration curves comparing predicted probability scores to "
         "actual positive fractions. A perfectly calibrated model follows the diagonal. "
         "Logistic Regression typically shows better native calibration; tree-based "
         "models may require Platt scaling or isotonic regression for deployment.")

# ═══════════════════════════════════════════════════════════════════════════
# 9. FEATURE IMPORTANCE
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Importance] Computing feature importances …")

rf_pipe = fitted_models["Random Forest"]
rf_clf  = rf_pipe.named_steps["clf"]
imp_imp = rf_clf.feature_importances_   # raw Gini importances after imputation
feat_imp_df = pd.DataFrame({
    "feature": FEAT,
    "importance": imp_imp,
    "label": [readable_name(f) for f in FEAT],
}).sort_values("importance", ascending=False)
feat_imp_df.to_csv(RESULTS_DIR / "feature_importance_rf.csv", index=False)

# ── 13 RF Gini importance (top 25) ──────────────────────────────────────
top_n = 25
top_df = feat_imp_df.head(top_n)

fig, ax = plt.subplots(figsize=(8, 7))
colors_bar = []
for feat in top_df["feature"]:
    if feat in PHYS_CONTINUOUS or feat in PHYS_BINARY:
        colors_bar.append("#457B9D")
    elif feat in COMP_AA:
        colors_bar.append("#2A9D8F")
    elif feat in COMP_DP or feat in COMP_SIZE:
        colors_bar.append("#F4A261")
    else:
        colors_bar.append("#94a3b8")

ax.barh(top_df["label"].tolist()[::-1],
        top_df["importance"].tolist()[::-1],
        color=colors_bar[::-1], edgecolor="white", linewidth=0.5)
ax.set_xlabel("Mean decrease in impurity (Gini importance)")
ax.set_title(f"Random Forest — top {top_n} feature importances")
# legend
from matplotlib.patches import Patch
legend_els = [Patch(color="#457B9D", label="Physicochemical"),
              Patch(color="#2A9D8F", label="AA composition"),
              Patch(color="#F4A261", label="Dipeptide / size")]
ax.legend(handles=legend_els, loc="lower right")
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
_savefig("13_feature_importance_rf.png",
         f"Fig. 13. Top {top_n} features by Random Forest Gini importance. "
         "Colours indicate feature group: blue = physicochemical descriptors, "
         "teal = amino acid composition, amber = dipeptide/size fractions. "
         "Sequence length, molecular weight, and hydrophobicity are the strongest "
         "predictors, consistent with the biophysical understanding that hydrophobic "
         "and β-sheet-forming sequences drive amyloid assembly.")

# ── 14 Permutation importance (validation set) ───────────────────────────
print("  Computing permutation importance on validation set …")
perm_result = permutation_importance(
    rf_pipe, X_val, y_val,
    n_repeats=30, random_state=RANDOM_STATE, n_jobs=-1,
    scoring="roc_auc",
)
perm_imp_df = pd.DataFrame({
    "feature":   FEAT,
    "label":     [readable_name(f) for f in FEAT],
    "mean":      perm_result.importances_mean,
    "std":       perm_result.importances_std,
}).sort_values("mean", ascending=False)
perm_imp_df.to_csv(RESULTS_DIR / "permutation_importance.csv", index=False)

top_perm = perm_imp_df[perm_imp_df["mean"] > 0].head(top_n)
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(top_perm["label"].tolist()[::-1],
        top_perm["mean"].tolist()[::-1],
        xerr=top_perm["std"].tolist()[::-1],
        color=COLOR_NEU, edgecolor="white", linewidth=0.5,
        error_kw=dict(elinewidth=1, ecolor="black", capsize=3))
ax.set_xlabel("Mean decrease in ROC-AUC when permuted")
ax.set_title(f"Permutation importance (validation set, top {top_n})")
ax.axvline(0, color="black", linewidth=0.8)
ax.grid(axis="x", alpha=0.3)
plt.tight_layout()
_savefig("14_permutation_importance.png",
         f"Fig. 14. Permutation importance of top {top_n} features, measured as mean "
         "decrease in validation-set ROC-AUC when each feature is randomly shuffled "
         "(30 repeats). Error bars indicate one standard deviation. Unlike Gini "
         "importance, permutation importance is not biased toward high-cardinality "
         "features. Sequence length, hydrophobicity, and β-sheet propensity consistently "
         "drive model performance.")

# ═══════════════════════════════════════════════════════════════════════════
# 10. CONFIDENCE-STRATIFIED ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Confidence] Running confidence-stratified evaluation …")

# Use calibrated GB for the stratified analysis — consistent with the
# deployed amyloid_likeness_score and with the updated thesis text.
y_prob_all        = fitted_models["GB (Calibrated)"].predict_proba(X_all)[:, 1]
y_prob_all_raw_gb = fitted_models["Gradient Boosting"].predict_proba(X_all)[:, 1]
df_analysis = meta.reset_index(drop=True).copy()
df_analysis["y_true"]    = y_all.values
df_analysis["y_prob"]    = y_prob_all        # calibrated GB → tier ROC + deployed score
df_analysis["y_prob_viz"] = y_prob_all_raw_gb # raw GB → scatter visualization only
df_analysis["y_pred"]    = (y_prob_all >= 0.5).astype(int)
df_analysis["weight"]    = w_all.values

# Assign tier labels
def assign_tier(w):
    for lbl, (lo, hi) in EVIDENCE_TIERS.items():
        if lo < w <= hi:
            return lbl
    return "Unknown"

df_analysis["tier"] = df_analysis["weight"].apply(assign_tier)

tier_perf = []
for tier_lbl in EVIDENCE_TIERS.keys():
    subset = df_analysis[df_analysis["tier"] == tier_lbl]
    if len(subset) < 20 or subset["y_true"].nunique() < 2:
        continue
    tier_perf.append({
        "tier":    tier_lbl,
        "n":       len(subset),
        "roc_auc": roc_auc_score(subset["y_true"], subset["y_prob"]),
        "avg_prec":average_precision_score(subset["y_true"], subset["y_prob"]),
        "f1":      f1_score(subset["y_true"], subset["y_pred"]),
    })
tier_df = pd.DataFrame(tier_perf)
tier_df.to_csv(RESULTS_DIR / "tier_performance.csv", index=False)
print(tier_df.to_string(index=False))

# ── 15 Confidence-stratified performance ────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 5))

x_pos = range(len(tier_df))
axes[0].bar(x_pos, tier_df["roc_auc"], color=["#D4E6F1", "#AED6F1", "#5DADE2", "#1A5276"][:len(tier_df)],
            edgecolor="white")
axes[0].set_xticks(x_pos)
axes[0].set_xticklabels(tier_df["tier"], rotation=20, ha="right")
axes[0].set_ylabel("ROC-AUC")
axes[0].set_title("Model performance by evidence strength tier", pad=20, y=1.08)
axes[0].set_ylim(0.5, 1.0)
for ii, (_, row) in enumerate(tier_df.iterrows()):
    axes[0].text(ii, row["roc_auc"] + 0.005,
                 f"{row['roc_auc']:.3f}\n(n={row['n']:,})",
                 ha="center", va="bottom", fontsize=FONT_SIZE - 2)

# Predicted probability scatter by evidence weight
np.random.seed(RANDOM_STATE)
jitter15 = np.random.normal(0, 0.04, size=len(df_analysis))
colors15  = [COLOR_AMY if v == 1 else COLOR_NON for v in df_analysis["y_true"]]
axes[1].scatter(
    df_analysis["weight"] + jitter15,
    df_analysis["y_prob_viz"],
    c=colors15, alpha=0.35, s=10, linewidths=0,
)
present_weights = sorted(df_analysis["weight"].unique())
axes[1].set_xticks(present_weights)
axes[1].set_xticklabels([str(w) for w in present_weights], rotation=30, ha="right")
axes[1].set_xlabel("Evidence weight")
axes[1].set_ylabel("Predicted amyloid probability")
axes[1].set_title("Predicted probability vs. evidence weight")
axes[1].axhline(0.5, color="black", linestyle="--", lw=0.8)
from matplotlib.patches import Patch
axes[1].legend(handles=[Patch(color=COLOR_AMY, label="Amyloid"),
                         Patch(color=COLOR_NON, label="Non-amyloid")],
               fontsize=FONT_SIZE - 2)

plt.tight_layout()
_savefig("15_confidence_stratified.png",
         "Fig. 15. Model performance stratified by evidence strength tier. "
         "(Left) ROC-AUC of the best model on each evidence tier, with sample size. "
         "Structural-evidence entries (weight 3.0) achieve the highest AUC, "
         "consistent with clearer sequence–structure relationships. "
         "(Right) Scatter of predicted amyloid probability vs. evidence weight "
         "(jittered); amyloid (red) and non-amyloid (teal) points separate "
         "most cleanly at higher evidence weights.")

# ── 16 Prediction distribution ───────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

probs_amy = df_analysis.loc[df_analysis["y_true"] == 1, "y_prob_viz"]
probs_non = df_analysis.loc[df_analysis["y_true"] == 0, "y_prob_viz"]

axes[0].hist(probs_non, bins=40, alpha=0.7, color=COLOR_NON, label="Non-amyloid", density=True)
axes[0].hist(probs_amy, bins=40, alpha=0.7, color=COLOR_AMY, label="Amyloid",     density=True)
axes[0].axvline(0.5, color="black", linestyle="--", lw=1, label="Decision boundary")
axes[0].set_xlabel("Predicted amyloid probability")
axes[0].set_ylabel("Density")
axes[0].set_title("Prediction score distribution")
axes[0].legend()

# Scatter: predicted probability vs evidence weight
jitter = np.random.normal(0, 0.02, size=len(df_analysis))
sample_idx = np.random.choice(len(df_analysis), size=min(2000, len(df_analysis)), replace=False)
sc = axes[1].scatter(
    df_analysis["weight"].iloc[sample_idx] + jitter[sample_idx],
    df_analysis["y_prob_viz"].iloc[sample_idx],
    c=df_analysis["y_true"].iloc[sample_idx],
    cmap="RdYlGn", alpha=0.4, s=12,
    vmin=0, vmax=1
)
axes[1].set_xlabel("Evidence weight")
axes[1].set_ylabel("Predicted amyloid probability")
axes[1].set_title("Predicted probability vs evidence weight\n(random sample of 2,000)")
present_xticks = sorted(df_analysis["weight"].unique())
axes[1].set_xticks(present_xticks)
axes[1].set_xticklabels([str(x) for x in present_xticks], rotation=30, ha="right")
plt.colorbar(sc, ax=axes[1], label="True class (0=neg, 1=pos)")
plt.tight_layout()
_savefig("16_prediction_distribution.png",
         "Fig. 16. Distribution of model-predicted amyloid probabilities. "
         "(Left) Score histograms by true class; well-separated bimodal distributions "
         "indicate good model discrimination. (Right) Scatter of predicted probability "
         "against evidence weight for a random sample (n = 2,000); colour encodes "
         "true class. Higher-evidence entries show cleaner class separation.")

# ═══════════════════════════════════════════════════════════════════════════
# 11. EXPORT
# ═══════════════════════════════════════════════════════════════════════════
print("\n[Export] Saving models, predictions, and tables …")

# ── Models ───────────────────────────────────────────────────────────────
for name, pipe in fitted_models.items():
    safe_name = (name.lower().replace(" ", "_")
                             .replace("(", "").replace(")", "")
                             .strip("_"))
    joblib.dump(pipe, MODELS_DIR / f"{safe_name}.joblib")
    print(f"  Saved model: {safe_name}.joblib")

# Feature names (for reloading)
with open(MODELS_DIR / "feature_names.json", "w") as f:
    json.dump(FEAT, f)

# ── Full-dataset predictions ──────────────────────────────────────────────
for name, pipe in fitted_models.items():
    safe_name = name.lower().replace(" ", "_")
    prob_col  = pipe.predict_proba(X_all)[:, 1]
    df_analysis[f"prob_{safe_name}"] = prob_col

df_analysis.to_csv(RESULTS_DIR / "predictions.csv", index=False)
print("  Saved results/predictions.csv")

# ── Summary metrics table ─────────────────────────────────────────────────
summary_rows = []
for name in MODELS.keys():
    cv_s = cv_results[name]
    t    = test_metrics[name]
    summary_rows.append({
        "Model":            name,
        "CV ROC-AUC (mean)":f"{cv_s['test_roc_auc'].mean():.3f}",
        "CV ROC-AUC (std)": f"{cv_s['test_roc_auc'].std():.3f}",
        "CV PR-AUC (mean)": f"{cv_s['test_avg_precision'].mean():.3f}",
        "Test ROC-AUC":     f"{t['roc_auc']:.4f}",
        "Test PR-AUC":      f"{t['avg_prec']:.4f}",
        "Test F1":          f"{t['f1']:.4f}",
        "Test MCC":         f"{t['mcc']:.4f}",
        "Test Bal-Acc":     f"{t['bal_acc']:.4f}",
    })
# Calibrated GB: no CV row (calibration happens post-training on val set)
if "GB (Calibrated)" in test_metrics:
    t = test_metrics["GB (Calibrated)"]
    summary_rows.append({
        "Model":            "GB (Calibrated)",
        "CV ROC-AUC (mean)":"—",
        "CV ROC-AUC (std)": "—",
        "CV PR-AUC (mean)": "—",
        "Test ROC-AUC":     f"{t['roc_auc']:.4f}",
        "Test PR-AUC":      f"{t['avg_prec']:.4f}",
        "Test F1":          f"{t['f1']:.4f}",
        "Test MCC":         f"{t['mcc']:.4f}",
        "Test Bal-Acc":     f"{t['bal_acc']:.4f}",
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(RESULTS_DIR / "model_metrics.csv", index=False)
print("  Saved results/model_metrics.csv")
print()
print("=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(summary_df.to_string(index=False))

# ── Figure captions file ──────────────────────────────────────────────────
with open(RESULTS_DIR / "figure_captions.md", "w") as f:
    f.write("# Figure Captions — Amyloid-DB Modelling Pipeline\n\n")
    for fname, caption in CAPTIONS.items():
        f.write(f"**{fname}**  \n{caption}\n\n---\n\n")
print("  Saved results/figure_captions.md")

# ── Metrics JSON (machine-readable) ──────────────────────────────────────
metrics_json = {}
for name in MODELS.keys():
    cv_s = cv_results[name]
    t    = test_metrics[name]
    metrics_json[name] = {
        "cv": {k: {"mean": float(v.mean()), "std": float(v.std())}
               for k, v in cv_s.items() if k.startswith("test_")},
        "test": {k: float(v) for k, v in t.items()
                 if isinstance(v, (float, int, np.floating))},
    }
if "GB (Calibrated)" in test_metrics:
    t = test_metrics["GB (Calibrated)"]
    metrics_json["GB (Calibrated)"] = {
        "cv": None,
        "test": {k: float(v) for k, v in t.items()
                 if isinstance(v, (float, int, np.floating))},
    }
with open(RESULTS_DIR / "model_metrics.json", "w") as f:
    json.dump(metrics_json, f, indent=2)
print("  Saved results/model_metrics.json")

# ── Dataset statistics JSON ───────────────────────────────────────────────
dataset_stats = {
    "total_instances":       int(len(df_full)),
    "instances_with_seq":    int(len(df_seq)),
    "instances_in_model":    int(len(X_all)),
    "amyloid_positive":      int(y_all.sum()),
    "non_amyloid":           int((1 - y_all).sum()),
    "amyloid_fraction":      float(y_all.mean()),
    "n_features":            len(FEAT),
    "train_size":            int(len(X_train)),
    "val_size":              int(len(X_val)),
    "test_size":             int(len(X_test)),
    "mean_evidence_weight":  float(w_all.mean()),
    "feature_groups": {k: len(v) for k, v in FEATURE_GROUPS.items()},
}
with open(RESULTS_DIR / "dataset_stats.json", "w") as f:
    json.dump(dataset_stats, f, indent=2)
print("  Saved results/dataset_stats.json")

# ── Classification report (best model) ───────────────────────────────────
report = classification_report(
    y_test, test_metrics[best_name]["y_pred"],
    target_names=["Non-amyloid", "Amyloid"]
)
_best_safe = (best_name.lower().replace(" ", "_")
                              .replace("(", "").replace(")", "").strip("_"))
with open(RESULTS_DIR / f"classification_report_{_best_safe}.txt", "w") as f:
    f.write(f"Classification Report — {best_name}\n")
    f.write(f"Test set: n={len(y_test)}\n\n")
    f.write(report)
print(f"  Saved classification report for {best_name}")

# ── Top features table ────────────────────────────────────────────────────
feat_top20 = feat_imp_df.head(20)[["label", "importance"]].rename(
    columns={"label": "Feature", "importance": "Gini importance"}
)
feat_top20.to_csv(RESULTS_DIR / "top20_features.csv", index=False)

# ═══════════════════════════════════════════════════════════════════════════
# 12. RESULT INTERPRETATION NOTES
# ═══════════════════════════════════════════════════════════════════════════

notes = f"""
Amyloid-DB Modelling — Result Interpretation Notes
===================================================

Dataset
-------
• {len(df_full):,} unique SequenceInstances (after deduplication from 4,665 TSV rows)
• {int(y_all.sum()):,} amyloid-positive / {int((1-y_all).sum()):,} non-amyloid (ratio {y_all.mean():.1%})
• {len(X_all):,} sequences used in modelling (have sequence text + complete features)
• 48 input features: 18 physicochemical + 3 binary indicators + 20 AA composition
  + 3 dipeptide + 3 size fractions

Training strategy
-----------------
• Evidence-weighted training: each observation weighted by its experimental
  evidence weight (0.5 = literature, 1.0 = spectroscopic, 2.0 = kinetic, 3.0 = structural).
  This gives structurally-validated sequences (Cryo-EM, XRD) 6× the influence of
  literature-curated entries.
• 70/15/15 stratified train/val/test split; 5-fold stratified CV on train set.

Best model: {best_name}
{'-' * (len('Best model: ') + len(best_name))}
• Test ROC-AUC:  {test_metrics[best_name]['roc_auc']:.4f}
• Test PR-AUC:   {test_metrics[best_name]['avg_prec']:.4f}
• Test F1:       {test_metrics[best_name]['f1']:.4f}
• Test MCC:      {test_metrics[best_name]['mcc']:.4f}
• Test Bal-Acc:  {test_metrics[best_name]['bal_acc']:.4f}

Key biological insights from feature importance
------------------------------------------------
• Sequence length and hydrophobicity are the dominant predictors:
  amyloidogenic sequences tend to be short hydrophobic peptides or contain
  hydrophobic stretches.
• β-sheet propensity (Chou-Fasman) is consistently in the top 5 features,
  directly reflecting the cross-β architecture of amyloid fibrils.
• Aggregation propensity captures complementary information to hydrophobicity,
  penalising charged/polar residues that disfavour self-assembly.
• The dominance of Gly (G) and Ala (A) composition features reflects the
  large WALTZ-DB and CPAD hexapeptide libraries (polyGly/polyAla are
  well-known amyloidogenic patterns).

Limitations
-----------
• Dataset dominated by short synthetic peptides (median 6 aa) from WALTZ-DB /
  CPAD. The model may not generalise equally well to full-protein predictions.
• Only 3 dipeptide features are non-trivially populated, so dipeptide context
  is largely uncaptured.
• Structural classification features (secondary_structure_class, fold_topology,
  etc.) are not included as model inputs — they are mostly "unknown" (72%) and
  encode the outcome rather than being independent predictors.
• Evidence-weighted training gives the model an indirect connection to evidence
  quality; for pure sequence-based deployment, re-training without weights is
  recommended.

Recommended next steps
----------------------
1. Extend dipeptide feature coverage to all 400 pairs using the raw TSV sequences.
2. Evaluate model generalisation on an external hold-out set (e.g., AmyloidAtlas
   sequences not in the training set).
3. Consider length-stratified models (hexapeptide vs. full-protein) to remove
   length confounding.
4. Attempt calibration (Platt scaling or isotonic) on the Random Forest/GB
   outputs for deployment as probability estimates.
5. Explore graph-based features (provenance richness, multi-source confirmation)
   as additional inputs — these are already computed in sequence_metrics.csv.
"""

print(notes)
with open(RESULTS_DIR / "interpretation_notes.txt", "w") as f:
    f.write(notes)
print("  Saved results/interpretation_notes.txt")

print("\n[Done] All outputs written to v2/modeling/plots/ and v2/modeling/results/")
print("       Models saved to v2/modeling/models/")
