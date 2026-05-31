# Mapping the Sequence–Structure–Function Space of Amyloids: A Data-Centric Study of Aggregation Criteria

A complete data infrastructure for amyloidogenicity research: an integration pipeline for 8 public databases, a machine learning classifier with calibrated probability scores, and a property graph representation for evidence provenance queries.

---

## Project Overview

The study of amyloid proteins is hampered by data heterogeneity — dozens of databases use different experimental methods, vocabularies, and definitions of "amyloidogenicity." This project builds a unified, evidence-weighted resource by:

1. **Integrating** 8 public amyloid databases into a single deduplicated dataset
2. **Standardizing** experimental evidence into a quantitative confidence tier (0.5–3.0)
3. **Resolving conflicts** using a weighted-consensus algorithm
4. **Predicting** amyloidogenicity with an evidence-weighted ML classifier and calibrated probability score (`amyloid_likeness_score`)
5. **Representing** all cross-database provenance as a queryable knowledge graph

---

## Repository Structure

```
v2/
├── amyloid_pipeline/          # Integration pipeline (Python package)
│   ├── parsers/               # 8 database-specific parsers
│   ├── features/              # Physicochemical + composition feature computation
│   ├── export/                # TSV, FASTA, SQLite exporters
│   ├── config.py              # Evidence weights, method mapping rules
│   ├── models.py              # AmyloidEntry dataclass
│   ├── filters.py             # Chainable filtering system
│   ├── unifier.py             # Deduplication and consensus resolution
│   └── run_pipeline.py        # Pipeline orchestration
│
├── database/                  # Pipeline outputs
│   ├── consensus_unified.tsv  # All 4,348 deduplicated entries with features
│   ├── amyloid_positive.tsv   # Amyloid-forming entries only
│   ├── non_amyloid.tsv        # Non-amyloid entries only
│   ├── amyloid_sequences.fasta
│   ├── non_amyloid_sequences.fasta
│   ├── amyloid_database.sqlite
│   └── conflicts.json         # Unresolved label conflicts
│
├── modeling_2/                # ML classification pipeline
│   ├── modeling_pipeline.py   # Complete pipeline: EDA → training → evaluation
│   ├── feature_engineering.py # Data loading and feature matrix preparation
│   ├── pipeline_config.py     # Feature lists, split sizes, evidence tiers
│   ├── models/                # Serialized sklearn models (.joblib)
│   ├── plots/                 # 16 EDA and evaluation figures
│   └── results/               # Metrics, predictions, feature importance
│
└── output/                    # Knowledge graph
    ├── graph_nodes.csv        # 10,846 nodes
    ├── graph_edges.csv        # 30,431 edges
    ├── amyloid_graph_turingdb.jsonl  # TuringDB-ready import file
    ├── sequence_metrics.csv   # Per-sequence uncertainty and provenance scores
    ├── schema.sql             # Normalized relational schema
    └── QUERY_MANUAL.md        # 28 annotated Cypher queries
```

---

## Component 1 — Integration Pipeline (`amyloid_pipeline/`)

### Supported Databases

| Database | Content |
|---|---|
| WALTZ-DB 2.0 | Experimentally characterized hexapeptides |
| Cross-Beta DB | High-quality cross-β structures |
| AmyLoad | Peptides with verified amyloid propensity |
| AmyloidExplorer | Structural data with disease context |
| AmyloidAtlas | Cryo-EM / NMR fibril structures |
| Amylobase | Aggregation kinetics data |
| AmyPro | Functional and pathogenic amyloids |
| AmyloGraph | Cross-seeding interactions |
| CPAD 2.0 (peptides + structures) | Aggregating peptides and structures |

### Evidence Weight Tiers

| Evidence type | Weight | Examples |
|---|---|---|
| Structural | 3.0 | Cryo-EM, XRD, ssNMR, microED |
| Kinetic | 2.0 | Aggregation kinetics, seeding assays |
| Staining / binding | 1.0 | ThT fluorescence, Congo Red |
| Literature-curated | 0.5 | Expert curation without direct assay |
| Computational | 0.0 | TANGO, WALTZ, PASTA predictions |

### Pipeline Workflow

```
Parse (8 formats) → Standardize → Evidence-weight → Deduplicate
    → Consensus resolution → Feature computation → Export
```

Deduplication key: `(sequence, uniprot_id, region_start, region_end)`.  
Conflict resolution: label with >50% higher cumulative weight wins; otherwise flagged in `conflicts.json`.

### Usage

```bash
pip install pandas openpyxl beautifulsoup4

python -m amyloid_pipeline \
    --waltzdb waltzdb.csv \
    --crossbeta crossbetadb.json \
    --amyload AmyLoad.csv \
    --amyloid-atlas AmyloidAtlas.tsv \
    --cpad-peptides aggregatingpeptides.xlsx \
    --cpad-structures amyloidstructure.xlsx \
    -o database/
```

### Python API

```python
from amyloid_pipeline import run_pipeline, AmyloidFilter

results = run_pipeline(
    input_files={'waltzdb': 'waltzdb.csv', 'crossbeta': 'data.json'},
    output_dir='database/',
    compute_features=True,
    export_sqlite=True
)

f = AmyloidFilter()
f.exclude_evidence_type('computational')
f.min_confidence(50)
filtered = f.apply(results['consensus'])
```

---

## Component 2 — ML Classification Pipeline (`modeling_2/`)

Binary classification of amyloid/non-amyloid sequences using 47 sequence-derived features, with evidence-weighted training and isotonic calibration.

### Dataset

- **3,594 experimentally supported sequences** (literature-tier entries excluded to avoid circular labelling)
- **Class balance:** 49.4% amyloid-positive / 50.6% non-amyloid
- **Split:** 70% train / 15% validation / 15% test, stratified, `random_state = 42`

### Results

| Model | CV ROC-AUC | Test ROC-AUC | Test F1 | Test MCC |
|---|---|---|---|---|
| Logistic Regression | 0.914 ± 0.008 | 0.903 | 0.818 | 0.653 |
| Random Forest | 0.945 ± 0.011 | 0.941 | 0.866 | 0.741 |
| **Gradient Boosting** | **0.951 ± 0.012** | **0.942** | **0.887** | **0.782** |
| GB (Calibrated) | — | 0.940 | 0.886 | 0.774 |

**GB (Calibrated)** — Gradient Boosting with isotonic calibration on the held-out validation set — is the deployed model. Its output probability is the `amyloid_likeness_score`.

### Usage

```bash
cd modeling_2/
pip install -r requirements.txt
python modeling_pipeline.py
```

### Load the deployed model

```python
import joblib, json

model    = joblib.load("modeling_2/models/gb_calibrated.joblib")
features = json.load(open("modeling_2/models/feature_names.json"))

score = model.predict_proba(X[features])[:, 1]  # amyloid_likeness_score
```

---

## Component 3 — Knowledge Graph (`output/`)

A property graph encoding cross-database evidence and provenance for all 4,348 sequence instances.

### Graph Statistics

| Metric | Value |
|---|---|
| Nodes | 10,846 |
| Edges | 30,431 |
| SequenceInstance nodes | 4,348 |
| Source databases | 8 |
| Proteins | 657 |
| Publications | 627 |
| Solved structures | 520 |

### Two-layer architecture

- **Consensus layer** — one `SequenceInstance` node per unique biological identity, carrying the evidence-weighted classification result
- **Observation layer** — one `Observation` node per source record, enabling cross-database agreement queries

### Load into TuringDB

```bash
cp output/amyloid_graph_turingdb.jsonl ~/.turing/data/
turingdb
> LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db
```

### Example query

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE p.protein_name CONTAINS 'TAU'
  AND s.is_amyloid = 'True'
  AND toFloat(s.confidence) >= 70
RETURN s.name, s.sequence, s.confidence, s.secondary_structure_class
ORDER BY toFloat(s.confidence) DESC
```

For 28 annotated queries see `output/QUERY_MANUAL.md`.

---

## Database Outputs (`database/`)

| File | Description |
|---|---|
| `consensus_unified.tsv` | 4,348 deduplicated entries with all 47 features |
| `amyloid_positive.tsv` | Amyloid-forming entries only |
| `non_amyloid.tsv` | Non-amyloid entries only |
| `amyloid_sequences.fasta` | Amyloid sequences in FASTA format |
| `non_amyloid_sequences.fasta` | Non-amyloid sequences in FASTA format |
| `amyloid_database.sqlite` | Indexed SQLite database with pre-built views |
| `conflicts.json` | 4 unresolved label conflict entries |

### SQLite quick start

```python
import sqlite3, pandas as pd

conn = sqlite3.connect('database/amyloid_database.sqlite')

df = pd.read_sql_query('''
    SELECT e.protein_name, e.sequence, p.hydrophobicity_mean, p.beta_propensity_mean
    FROM entries e
    JOIN physicochemical_features p ON e.id = p.entry_id
    WHERE e.organism LIKE '%Homo sapiens%'
      AND e.pathogenicity = 'pathogenic'
      AND e.structure_type = 'fibril'
''', conn)
```

---

## Installation

```bash
# Pipeline dependencies
pip install pandas openpyxl beautifulsoup4

# ML pipeline dependencies
pip install -r modeling_2/requirements.txt
```

---

*Data sources: WALTZ-DB 2.0, Cross-Beta DB, AmyLoad, AmyloidExplorer, AmyloidAtlas,
Amylobase, AmyPro, AmyloGraph, CPAD 2.0.*
