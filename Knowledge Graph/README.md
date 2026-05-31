# Amyloid-DB — Knowledge Graph

A property graph representation of the integrated Amyloid-DB dataset, encoding
experimental evidence, provenance, and cross-database relationships for 4,348
unique amyloid/non-amyloid sequence instances.

---

## Contents

| File | Description |
|---|---|
| `graph_nodes.csv` | All 10,846 nodes with type and JSON attributes |
| `graph_edges.csv` | All 30,431 edges with type and JSON attributes |
| `graph.jsonl` | Edge-only JSONL (string node IDs; legacy format) |
| `amyloid_graph_turingdb.jsonl` | Complete graph in TuringDB JSONL format (nodes + edges; use this for import) |
| `graph_build.py` | Builds `graph_nodes.csv` and `graph_edges.csv` from `database/consensus_unified.tsv` |
| `convert_to_turingdb.py` | Converts CSV graph to `amyloid_graph_turingdb.jsonl` |
| `compute_metrics.py` | Computes per-sequence uncertainty and provenance metrics; writes to `sequence_metrics.csv` and merges into `graph_nodes.csv` |
| `sequence_metrics.csv` | Per-sequence uncertainty and provenance scores (one row per SequenceInstance) |
| `schema.sql` | Normalized 10-table relational schema (SQLite / PostgreSQL) |
| `schema.dbml` | Same schema in DBML format for [dbdiagram.io](https://dbdiagram.io) |
| `db_scheme.md` | Schema table reference with column descriptions and notes |
| `GRAPH_NOTES.md` | Detailed graph architecture, statistics, layer semantics, and known limitations |
| `SCHEMA_NOTES.md` | Schema design decisions, ambiguities, and open recommendations |
| `QUERY_MANUAL.md` | 28 annotated Cypher queries for TuringDB (with biological interpretation for each) |
| `turingdb_layers.md` | Layer-by-layer query guide for TuringDB visual exploration |

---

## Graph Statistics

| Metric | Value |
|---|---|
| Total nodes | 10,846 |
| Total edges | 30,431 |
| Unique sequence instances | 4,348 |
| Source databases | 8 |
| Proteins | 657 |
| Publications | 627 |
| Solved structures | 520 |
| Experiment methods | 15 |

---

## Graph Architecture

The graph has **two coexisting layers**, both present in all output files and
differentiated by the `graph_layer` attribute on every node and edge.

### Layer 1 — Consensus Layer (`graph_layer = 'consensus'`)

Deduplicated, evidence-weighted view. One `SequenceInstance` node per unique
biological identity `(sequence, uniprot_id, region_start, region_end)`.

```
SequenceInstance ──PART_OF──────► Protein
SequenceInstance ──OBSERVED_IN──► SourceDatabase
SequenceInstance ──SUPPORTED_BY─► ExperimentMethod
SequenceInstance ──REPORTED_IN──► Publication
SequenceInstance ──HAS_STRUCTURE► Structure
```

**Use for:** clean analysis, classification queries, feature export.

### Layer 2 — Observation Layer (`graph_layer = 'observation'`)

Per-source provenance record, before consensus merging. One `Observation` node
per row in `consensus_unified.tsv`.

```
Observation ──OBSERVED_AS──► SequenceInstance
Observation ──FROM_SOURCE──► SourceDatabase
Observation ──SUPPORTED_BY─► ExperimentMethod
Observation ──REPORTED_IN──► Publication
```

**Use for:** tracing which database reported what, cross-database agreement checks.

A third sub-layer (`graph_layer = 'conflict'`) contains 6 stub `Observation` nodes
for the 3 sequence instances where two databases reported conflicting labels. These
are marked `is_consensus_winner = 'false'` and have incomplete attributes.

---

## Node Types

| Type | Count | Key properties |
|---|---|---|
| `SequenceInstance` | 4,348 | `sequence`, `is_amyloid`, `confidence`, `evidence_weight`, `evidence_type`, `secondary_structure_class`, `disease`, `n_observations`, `provenance_richness` |
| `Observation` | 4,671 | `source_db`, `experimental_label`, `raw_method`, `evidence_weight`, `record_id`, `is_consensus_winner` |
| `Protein` | 657 | `uniprot_id`, `protein_name`, `organism`, `protein_family` |
| `Publication` | 627 | `doi`, `pmid`, `reference` |
| `Structure` | 520 | `pdb_id`, `emdb_id`, `resolution`, `method_universal` |
| `ExperimentMethod` | 15 | `method_universal`, `evidence_type`, `evidence_weight` |
| `SourceDatabase` | 8 | `source_db` |

---

## Evidence Weight Tiers

Evidence weights assigned per experimental method:

| Method type | `evidence_weight` | Examples |
|---|---|---|
| Structural | 3.0 | Cryo-EM, XRD, microED, ssNMR |
| Kinetic | 2.0 | Aggregation kinetics |
| Staining / binding | 1.0 | ThT binding, Congo Red |
| Literature-curated | 0.5 | Literature review |
| Computational | 0.0 | Prediction tools |

---

## Sequence-Level Metrics (`sequence_metrics.csv`)

Computed by `compute_metrics.py` and stored on every `SequenceInstance` node:

| Metric | Description |
|---|---|
| `n_observations` | Total observations linked to this sequence |
| `n_sources` | Number of distinct source databases |
| `n_methods` | Number of distinct experimental methods |
| `n_labels` | Distinct labels reported (1 = agreement, 2+ = disagreement) |
| `n_publications` | Distinct publications |
| `conflict_score` | Fraction of observations opposing the consensus label (0–1) |
| `consensus_support_weight` | Total evidence weight supporting the consensus label |
| `method_diversity` | Simpson's diversity index across methods (0–1) |
| `provenance_richness` | Composite log-scaled score: `log2(1+n_sources) + log2(1+n_methods) + log2(1+n_publications)` |

**Informative ranges in this dataset:** `conflict_score > 0` for 3 sequences;
`provenance_richness` ranges 2.0–4.76; `n_sources ≥ 2` for 271 sequences.

---

## Loading into TuringDB

```bash
cp amyloid_graph_turingdb.jsonl ~/.turing/data/
turingdb
```
```
> LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db
```

> **Note:** Use `amyloid_graph_turingdb.jsonl` — not `graph.jsonl`.
> The latter is a legacy edge-only file with string IDs and is not
> compatible with direct TuringDB import.

All properties are stored as **strings** in TuringDB. Use explicit casts for comparisons:
- `toFloat(n.confidence) >= 70`
- `toInteger(n.region_start)`
- `n.is_amyloid = 'True'` (not `true`)

---

## Example Queries

### High-confidence amyloid regions of α-synuclein
```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE p.protein_name CONTAINS 'SYUA'
  AND s.is_amyloid = 'True'
  AND toFloat(s.confidence) >= 70
RETURN s.name, s.sequence, s.confidence, s.secondary_structure_class
ORDER BY toFloat(s.confidence) DESC
```

### Sequences confirmed by two or more databases
```cypher
MATCH (s:SequenceInstance)-[:OBSERVED_IN]->(db:SourceDatabase)
WHERE s.graph_layer = 'consensus'
WITH s, count(db) AS n_db
WHERE n_db >= 2
RETURN s.name, s.is_amyloid, n_db
ORDER BY n_db DESC
```

### Full evidence trail for one sequence
```cypher
MATCH (s:SequenceInstance)
WHERE s.name = 'TAU_HUMAN [306–311]'
OPTIONAL MATCH (s)-[:SUPPORTED_BY]->(m:ExperimentMethod)
OPTIONAL MATCH (s)-[:OBSERVED_IN]->(db:SourceDatabase)
OPTIONAL MATCH (s)-[:REPORTED_IN]->(pub:Publication)
OPTIONAL MATCH (s)-[:HAS_STRUCTURE]->(st:Structure)
RETURN s.is_amyloid, s.confidence, m.method_universal, db.source_db,
       pub.reference, st.pdb_id, st.resolution
```

For 28 annotated queries covering proteins, evidence tracing, cross-database
comparison, structures, and publications — see **`QUERY_MANUAL.md`**.

---

## Regenerating the Graph

```bash
# Step 1 — build nodes and edges from TSV
python graph_build.py

# Step 2 — add sequence-level metrics
python compute_metrics.py

# Step 3 — convert to TuringDB format
python convert_to_turingdb.py
```

Input: `../database/consensus_unified.tsv` and `../database/conflicts.json`
Output: `graph_nodes.csv`, `graph_edges.csv`, `sequence_metrics.csv`, `amyloid_graph_turingdb.jsonl`

---

## Key Limitations

- **Only consensus winners are captured as full Observations.** Within each deduplication
  group, non-winning records are discarded before TSV export and have no Observation nodes.
- **~26% of entries have no UniProt ID** — deduplication falls back to protein name or
  sequence + source, which can merge biologically distinct entries.
- **72% of SequenceInstances have `secondary_structure_class = 'unknown'`** — structural
  classification is inferred from text fields, not curated annotation.
- **ThT binding dominates** — 2,030 CPAD entries all use ThT binding; the `ExperimentMethod`
  node for ThT has ~4,000 edges and will dominate degree-based graph analytics.

For the full list of known issues and design decisions, see **`SCHEMA_NOTES.md`**
and **`GRAPH_NOTES.md`**.

---

*Data sources: CPAD 2.0, WALTZ-DB, AmyPro, AmyloGraph, AmyloidAtlas, AmyloidExplorer,
CPAD-Structures, Cross-Beta DB. Graph generated: 2026-04-09.*

