# Amyloid-DB — Machine Learning Pipeline

Binary classification of amyloid/non-amyloid sequences using sequence-derived
physicochemical and compositional features. Includes evidence-weighted training,
isotonic calibration, and stratified evaluation by experimental evidence tier.

This folder contains the updated modelling pipeline described in
*Mapping the Sequence–Structure–Function Space of Amyloids: A Data-Centric Study
of Aggregation Criteria* — with two key improvements over the initial version:
(1) exclusion of literature-tier entries to eliminate circular labelling artefacts,
and (2) isotonic calibration of the best model to produce a deployable
`amyloid_likeness_score`.

---

## Key Results

| Model | CV ROC-AUC | Test ROC-AUC | Test PR-AUC | Test F1 | Test MCC |
|---|---|---|---|---|---|
| Logistic Regression | 0.914 ± 0.008 | 0.903 | 0.912 | 0.818 | 0.653 |
| Random Forest | 0.945 ± 0.011 | 0.941 | 0.949 | 0.866 | 0.741 |
| **Gradient Boosting** | **0.951 ± 0.012** | **0.942** | **0.948** | **0.887** | **0.782** |
| GB (Calibrated) | — | 0.940 | 0.937 | 0.886 | 0.774 |

Test set: n = 540 (267 amyloid-positive, 273 non-amyloid).  
**GB (Calibrated)** is the deployed model — Gradient Boosting with isotonic
calibration applied on the held-out validation set. Its output is the
`amyloid_likeness_score` stored in the database.

**Performance by evidence tier (GB Calibrated):**

| Tier | n | ROC-AUC | PR-AUC | F1 |
|---|---|---|---|---|
| Spectroscopic (weight 1.0) | 2,032 | 0.970 | 0.959 | 0.920 |
| Structural (weight 3.0) | 1,302 | 0.995 | 0.995 | 0.974 |

---

## Dataset

- **Source:** `../database/consensus_unified.tsv` (4,348 unique SequenceInstances)
- **Modelling subset:** 3,594 sequences — experimentally supported entries with
  complete sequence features (`evidence_weight ≥ 1.0`)
- **Excluded:** 397 literature-tier entries (`evidence_weight < 1.0`), comprising
  AmyPro literature-curated records and WALTZ-DB entries whose labels were assigned
  by computational tools (WALTZ/TANGO) using the same physicochemical features as
  the model inputs. Excluding these removes the circular-labelling artefact.
- **Class balance:** 1,776 amyloid-positive (49.4%) / 1,818 non-amyloid

### Data split

| Subset | n | Purpose |
|---|---|---|
| Train | 2,515 (70%) | Model fitting and 5-fold cross-validation |
| Validation | 539 (15%) | Permutation importance; isotonic calibration of GB |
| Test | 540 (15%) | Final held-out evaluation (never seen during training) |

Stratified by class label, `random_state = 42`.

---

## Features (47 total)

| Group | n | Description |
|---|---|---|
| Physicochemical (continuous) | 18 | Hydrophobicity, β-sheet propensity, aggregation propensity, charge, aromatic / aliphatic / polar / charged fractions, molecular weight, length |
| Physicochemical (binary) | 3 | `has_polyq`, `has_glycine_rich`, `has_proline` |
| AA composition | 20 | Single-residue frequencies for all 20 standard amino acids |
| Dipeptide | 3 | AA, AQ, QA dipeptide frequencies (only non-trivially populated pairs) |
| Size fractions | 3 | Tiny (A, G, S), small, large residue fractions |

**Top predictors** (Random Forest Gini): sequence length · molecular weight ·
β-sheet propensity · aggregation propensity · Val frequency · hydrophobicity (mean/max).

---

## Training Strategy

- **Literature-tier filtering:** entries with `evidence_weight < 1.0` are excluded
  before training. See `feature_engineering.py`.
- **Evidence-weighted training:** each sequence is weighted by its `evidence_weight`
  during model fitting — structural entries (Cryo-EM/XRD, weight 3.0) receive 6×
  the influence of spectroscopic assays (ThT, weight 1.0).
- **sklearn Pipeline:** `SimpleImputer → StandardScaler → Classifier` — preprocessing
  is fit inside each CV fold to prevent data leakage.
- **Isotonic calibration:** `CalibratedClassifierCV(method='isotonic', cv='prefit')`
  applied to the fitted Gradient Boosting model on the held-out validation set.
  Calibrated probabilities reflect observed positive rates in experimental data.

---

## Files

### Scripts

| File | Description |
|---|---|
| `modeling_pipeline.py` | Complete pipeline: EDA → training → evaluation → export (run this) |
| `feature_engineering.py` | Data loading, literature-tier filtering, and feature matrix preparation |
| `pipeline_config.py` | All configuration: paths, feature lists, split sizes, evidence tier boundaries, colour palette |
| `requirements.txt` | Python package dependencies |

### Models (`models/`)

| File | Description |
|---|---|
| `gb_calibrated.joblib` | **Deployed model.** Isotonically calibrated GB — outputs `amyloid_likeness_score` |
| `gradient_boosting.joblib` | Fitted sklearn Pipeline: imputer + scaler + GradientBoostingClassifier |
| `random_forest.joblib` | Fitted sklearn Pipeline: imputer + scaler + RandomForestClassifier |
| `logistic_regression.joblib` | Fitted sklearn Pipeline: imputer + scaler + LogisticRegression |
| `feature_names.json` | Ordered list of 47 feature column names used during training |

Load a model:
```python
import joblib, json, pandas as pd

model    = joblib.load("models/gb_calibrated.joblib")
features = json.load(open("models/feature_names.json"))

# X must be a DataFrame with columns matching features
score = model.predict_proba(X[features])[:, 1]  # amyloid_likeness_score
```

### Results (`results/`)

| File | Description |
|---|---|
| `model_metrics.csv` | CV and test metrics for all 4 models (human-readable) |
| `model_metrics.json` | Same metrics in JSON with full precision |
| `classification_report_gradient_boosting.txt` | Per-class precision, recall, F1 for GB on the test set |
| `tier_performance.csv` | ROC-AUC, PR-AUC, F1 by evidence tier (GB Calibrated) |
| `dataset_stats.json` | Dataset counts: total, modelling subset, class balance, feature counts, split sizes |
| `feature_importance_rf.csv` | Gini importance for all 47 features (Random Forest) |
| `permutation_importance.csv` | Permutation importance ± std for top features (validation set, 30 repeats) |
| `top20_features.csv` | Top 20 features by Gini importance with human-readable names |
| `predictions.csv` | Full dataset scores: `y_true`, `y_prob` (calibrated GB), `y_pred`, `weight`, `tier`, and raw probabilities for all 4 models |

### Plots (`plots/`)

| File | Description |
|---|---|
| `01_class_distribution.png` | Amyloid / non-amyloid counts and evidence weight by class |
| `02_sequence_length_distribution.png` | Sequence length distribution by class (linear + log) |
| `03_source_database_composition.png` | Sequences per source database; evidence type breakdown |
| `04_evidence_method_composition.png` | Method counts; evidence tier distribution |
| `05_confidence_distribution.png` | Confidence score and evidence weight distributions |
| `06_physicochemical_by_class.png` | Key physicochemical descriptors by amyloid class (box plots) |
| `07_feature_correlation.png` | Pearson correlation matrix of physicochemical features |
| `08_cv_performance.png` | 5-fold CV ROC-AUC and PR-AUC distributions, all models |
| `09_roc_curves.png` | ROC curves on test set, all models |
| `10_pr_curves.png` | Precision-Recall curves on test set, all models |
| `11_confusion_matrix.png` | Confusion matrix for Gradient Boosting on test set |
| `12_calibration_plot.png` | Calibration curves: predicted probability vs. observed positive rate |
| `13_feature_importance_rf.png` | Top 25 features by Random Forest Gini importance |
| `14_permutation_importance.png` | Top 25 features by permutation importance (validation set) |
| `15_confidence_stratified.png` | ROC-AUC by evidence tier; predicted probability vs. evidence weight |
| `16_prediction_distribution.png` | Score histograms by class; probability vs. evidence weight scatter |

---

## Reproducing the Results

```bash
# Install dependencies
pip install -r requirements.txt

# Run from the modeling_2/ directory
python modeling_pipeline.py
```

Outputs are written to `plots/`, `results/`, and `models/`.  
All results are reproducible with `random_state = 42`.

---

## Limitations

- **Short peptide dominance:** median sequence length is 6 aa (WALTZ-DB / CPAD
  hexapeptide libraries). Generalisation to full-protein predictions should be
  verified on an independent hold-out set.
- **Three dipeptide features only:** only AA, AQ, QA pairs are non-trivially
  populated; richer dipeptide context is uncaptured.
- **Structural classification excluded:** `secondary_structure_class` and related
  fields are not model inputs — they are 72% unknown and encode the outcome rather
  than being independent predictors.
- **Evidence-weighted training for deployment:** for scoring sequences without known
  evidence weights, re-training with uniform weights is recommended.

---

*Data: Amyloid-DB consensus dataset (8 integrated sources: CPAD, WALTZ-DB, AmyPro,
AmyloGraph, AmyloidAtlas, AmyloidExplorer, CPAD-Structures, Cross-Beta DB).
Pipeline last run: 2026-05-31.*

