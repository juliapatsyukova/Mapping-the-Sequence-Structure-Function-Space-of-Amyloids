# Mapping the Sequence–Structure–Function Space of Amyloids: A Data-Centric Study of Aggregation Criteria

This repository contains the analytical components of a data-centric study on amyloidogenicity: a machine learning classification pipeline and a knowledge graph representation of cross-database evidence. The project addresses data heterogeneity in amyloid research by building a unified, evidence-weighted resource from 8 public databases.

---

## Repository Structure

```
/
├── Database/           # Source data and initial pipeline (see note below)
├── Methods Mapping/    # Evidence standardization scripts (see note below)
├── Modelling/          # ML classification pipeline
└── Knowledge Graph/    # Property graph representation
```

---

## Database & Methods Mapping

The integrated Amyloid-DB database was developed in collaboration. The full pipeline — including parsers for all 8 source databases, evidence weighting, deduplication, conflict resolution, and SQLite export — is available in a separate repository:

> **[github.com/juliapatsyukova/Amyloid-DB/tree/merged_version/v2](https://github.com/juliapatsyukova/Amyloid-DB/tree/merged_version/v2)**

The `Database/` and `Methods Mapping/` folders in this repository contain earlier-stage scripts used during the initial data collection phase.

---

## Modelling

Binary classification of amyloid/non-amyloid sequences using 47 sequence-derived features, with evidence-weighted training and isotonic calibration to produce a deployable `amyloid_likeness_score`.

### Structure

```
Modelling/
├── modeling_pipeline.py    # Complete pipeline: EDA → training → evaluation → export
├── feature_engineering.py  # Data loading, literature-tier filtering, feature matrix
├── pipeline_config.py      # Feature lists, split sizes, evidence tier boundaries
├── requirements.txt
│
├── Models/
│   ├── gb_calibrated.joblib       # Deployed model (isotonically calibrated GB)
│   ├── gradient_boosting.joblib
│   ├── random_forest.joblib
│   ├── logistic_regression.joblib
│   └── feature_names.json
│
├── Plots/                  # 16 EDA and evaluation figures (300 DPI PNG)
│
└── Results/
    ├── model_metrics.csv
    ├── model_metrics.json
    ├── tier_performance.csv
    ├── dataset_stats.json
    ├── feature_importance_rf.csv
    ├── permutation_importance.csv
    ├── top20_features.csv
    └── predictions.csv
```

### Dataset

- **3,594 experimentally supported sequences** from the integrated database
- 397 literature-tier entries excluded (`evidence_weight < 1.0`) to avoid circular labelling — these entries carry labels assigned by computational prediction tools using the same physicochemical features as model inputs
- **Class balance:** 1,776 amyloid-positive (49.4%) / 1,818 non-amyloid
- **Split:** 70% train / 15% validation / 15% test, stratified, `random_state = 42`

### Features (47 total)

| Group | n | Description |
|---|---|---|
| Physicochemical | 21 | Hydrophobicity, β-sheet propensity, aggregation propensity, charge, aromatic / aliphatic / polar fractions, molecular weight, length, binary indicators |
| AA composition | 20 | Single-residue frequencies for all 20 standard amino acids |
| Dipeptide | 3 | AA, AQ, QA dipeptide frequencies |
| Size fractions | 3 | Tiny, small, large residue fractions |

### Results

| Model | CV ROC-AUC | Test ROC-AUC | Test PR-AUC | Test F1 | Test MCC |
|---|---|---|---|---|---|
| Logistic Regression | 0.914 ± 0.008 | 0.903 | 0.912 | 0.818 | 0.653 |
| Random Forest | 0.945 ± 0.011 | 0.941 | 0.949 | 0.866 | 0.741 |
| **Gradient Boosting** | **0.951 ± 0.012** | **0.942** | **0.948** | **0.887** | **0.782** |
| GB (Calibrated) | — | 0.940 | 0.937 | 0.886 | 0.774 |

**GB (Calibrated)** — Gradient Boosting with isotonic calibration on the held-out validation set — is the deployed model. Its output probability is the `amyloid_likeness_score`.

**Performance by evidence tier (GB Calibrated):**

| Tier | n | ROC-AUC | PR-AUC | F1 |
|---|---|---|---|---|
| Spectroscopic (weight 1.0) | 2,032 | 0.970 | 0.959 | 0.920 |
| Structural (weight 3.0) | 1,302 | 0.995 | 0.995 | 0.974 |

### Usage

```bash
cd Modelling/
pip install -r requirements.txt
python modeling_pipeline.py
```

### Load the deployed model

```python
import joblib, json

model    = joblib.load("Modelling/Models/gb_calibrated.joblib")
features = json.load(open("Modelling/Models/feature_names.json"))

score = model.predict_proba(X[features])[:, 1]  # amyloid_likeness_score (0–1)
```

---

## Knowledge Graph

A property graph encoding cross-database evidence and provenance for all 4,348 sequence instances, with two coexisting layers: a deduplicated consensus layer and a per-source observation layer.

### Structure

```
Knowledge Graph/
├── graph_nodes.csv               # 10,846 nodes with attributes
├── graph_edges.csv               # 30,431 edges with attributes
├── amyloid_graph_turingdb.jsonl  # TuringDB-ready import file
├── sequence_metrics.csv          # Per-sequence provenance and uncertainty scores
├── graph_build.py                # Builds nodes/edges from consensus TSV
├── compute_metrics.py            # Computes per-sequence provenance metrics
├── convert_to_turingdb.py        # Converts to TuringDB format
├── schema.sql                    # Normalized relational schema
├── schema.dbml                   # DBML schema for dbdiagram.io
├── GRAPH_NOTES.md                # Graph architecture and statistics
├── SCHEMA_NOTES.md               # Schema design decisions
└── QUERY_MANUAL.md               # 28 annotated Cypher queries
```

### Graph Statistics

| | |
|---|---|
| Total nodes | 10,846 |
| Total edges | 30,431 |
| SequenceInstance nodes | 4,348 |
| Observation nodes | 4,671 |
| Source databases | 8 |
| Proteins | 657 |
| Publications | 627 |
| Solved structures | 520 |

### Graph Architecture

- **Consensus layer** (`graph_layer = 'consensus'`): one `SequenceInstance` per unique biological identity `(sequence, uniprot_id, region_start, region_end)`, carrying the evidence-weighted classification result
- **Observation layer** (`graph_layer = 'observation'`): one `Observation` node per source record, enabling cross-database agreement and provenance queries

### Load into TuringDB

```bash
cp "Knowledge Graph/amyloid_graph_turingdb.jsonl" ~/.turing/data/
turingdb
> LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db
```

### Example queries

```cypher
-- All high-confidence amyloid regions of α-synuclein
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE p.protein_name CONTAINS 'SYUA'
  AND s.is_amyloid = 'True'
  AND toFloat(s.confidence) >= 70
RETURN s.name, s.sequence, s.confidence, s.secondary_structure_class
ORDER BY toFloat(s.confidence) DESC

-- Sequences confirmed by two or more databases
MATCH (s:SequenceInstance)-[:OBSERVED_IN]->(db:SourceDatabase)
WHERE s.graph_layer = 'consensus'
WITH s, count(db) AS n_db
WHERE n_db >= 2
RETURN s.name, s.is_amyloid, n_db
ORDER BY n_db DESC
```

For 28 annotated queries covering proteins, evidence tracing, cross-database comparison, structures, and publications — see `Knowledge Graph/QUERY_MANUAL.md`.

---

*Data sources: WALTZ-DB 2.0, Cross-Beta DB, AmyLoad, AmyloidExplorer, AmyloidAtlas, Amylobase, AmyPro, AmyloGraph, CPAD 2.0.*
