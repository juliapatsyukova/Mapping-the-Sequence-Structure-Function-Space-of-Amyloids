# Graph Notes — Amyloid Evidence / Provenance Graph

Generated: 2026-04-09
Data sources: `database/consensus_unified.tsv` (4,665 rows),
              `database/conflicts.json` (4 records)
Metrics added: 2026-04-09 — see Section 9

---

## 1. Graph Architecture — Two Layers

The graph has two coexisting layers, differentiated by the `graph_layer` attribute on
every node and edge. Both layers are present in the same output files.

```
LAYER 1 — Consensus Layer
───────────────────────────────────────────────────────────────────────────────
SequenceInstance ──PART_OF──► Protein
SequenceInstance ──OBSERVED_IN──► SourceDatabase
SequenceInstance ──SUPPORTED_BY──► ExperimentMethod
SequenceInstance ──REPORTED_IN──► Publication
SequenceInstance ──HAS_STRUCTURE──► Structure

LAYER 2 — Observation Layer
───────────────────────────────────────────────────────────────────────────────
Observation ──OBSERVED_AS──► SequenceInstance
Observation ──FROM_SOURCE──► SourceDatabase
Observation ──SUPPORTED_BY──► ExperimentMethod
Observation ──REPORTED_IN──► Publication
```

These layers use different edge types where semantics differ (`OBSERVED_IN` vs `OBSERVED_AS`,
`FROM_SOURCE`) and reuse the same edge type where the semantics are identical
(`SUPPORTED_BY`, `REPORTED_IN`). The `graph_layer` attribute on every edge disambiguates
which layer an edge belongs to in all cases.

---

## 2. Graph Statistics

### Combined (both layers)

| Metric | Value |
|---|---|
| Total nodes | 10,846 |
| Total edges | 30,431 |

### Node counts by type

| Node type | Count | Layer |
|---|---|---|
| SequenceInstance | 4,348 | consensus + conflict stubs |
| Observation | 4,671 | observation (4,665 TSV rows + 6 conflict stubs) |
| Protein | 657 | consensus |
| Publication | 627 | consensus + observation (merged) |
| Structure | 520 | consensus |
| ExperimentMethod | 15 | consensus + observation (merged) |
| SourceDatabase | 8 | consensus + observation (merged) |

### Edge counts by type and layer

| Edge type | Layer | Count | Meaning |
|---|---|---|---|
| OBSERVED_AS | observation | 4,665 | TSV-row Observation → SequenceInstance |
| OBSERVED_AS | conflict | 6 | Conflict-stub Observation → SequenceInstance |
| FROM_SOURCE | observation | 4,665 | TSV-row Observation → SourceDatabase |
| FROM_SOURCE | conflict | 6 | Conflict-stub Observation → SourceDatabase |
| OBSERVED_IN | consensus | 4,619 | SequenceInstance → SourceDatabase |
| SUPPORTED_BY | consensus | 4,627 | SequenceInstance → ExperimentMethod |
| SUPPORTED_BY | observation | 4,665 | Observation → ExperimentMethod |
| REPORTED_IN | consensus | 1,185 | SequenceInstance → Publication |
| REPORTED_IN | observation | 1,200 | Observation → Publication |
| PART_OF | consensus | 4,031 | SequenceInstance → Protein |
| HAS_STRUCTURE | consensus | 762 | SequenceInstance → Structure |

---

## 3. Layer 1 — Consensus Layer

### What it represents

The outcome of the pipeline's deduplication and evidence-weighted consensus resolution.
One `SequenceInstance` node per unique biological identity
`(sequence, uniprot_id, region_start, region_end)`.

### SequenceInstance nodes are source-agnostic

SequenceInstance nodes do NOT carry a `source_db` attribute. Which database(s) contributed
to a SequenceInstance is encoded only in:
- `OBSERVED_IN` edges (consensus layer): SequenceInstance → SourceDatabase
- `OBSERVED_AS` edges (observation layer): Observation → SequenceInstance

### Consensus attributes on SequenceInstance

The `is_amyloid`, `evidence_weight`, `confidence`, and `experimental_label` attributes on
a SequenceInstance reflect the consensus winner:
- For single-source entries: the one available observation.
- For multi-source entries with the same label: the highest-weight observation.
- For multi-source entries with conflicting labels (weight ratio > 1.5): the winning label's
  highest-weight observation.

The full 11-level structural classification (secondary_structure_class, interface_type, etc.)
is also carried on the SequenceInstance node and reflects the winning observation's parser inference.

---

## 4. Layer 2 — Observation Layer

### What it represents

The per-row provenance record, before consensus collapse. Each `Observation` node
corresponds to one row in `consensus_unified.tsv` — i.e., one consensus-winner record
per dedup-key group.

### Option A: Observation nodes (chosen over edge-attribute approach)

**Why not Option B (multiple attribute-rich parallel edges):**
Each TSV row is a coherent bundle — one source, one method, one publication, one structural
ID, one experimental label — all tied together. Scattering these across multiple parallel
edges between the same node pair makes it impossible to reconstruct which method, publication,
or PDB ID belongs to which source observation for a given SequenceInstance. Observation nodes
solve this by making each bundle a first-class entity that can be queried independently.

### Multi-source SequenceInstances

**300 SequenceInstance nodes have 2 or more Observation nodes** pointing at them. These
represent biologically important cases where the same region was independently characterized
in multiple databases. Example:

```
seq:f999f2133b62be4c  (SequenceInstance)
  ← OBSERVED_AS ← obs:2d75235f1b52e589  (source=CPAD,     label=non-amyloid, method=ThT binding)
  ← OBSERVED_AS ← obs:38eaaea74bcb6893  (source=WALTZ-DB, label=non-amyloid, method=Experimental)
```

This structure allows queries like: "Show me all sequences that both CPAD and WALTZ-DB
independently characterized with the same label" — a cross-database agreement query that
is not possible from the consensus layer alone.

### Observation node attributes

| Attribute | Description |
|---|---|
| `record_id` | Original identifier from the source DB |
| `source_db` | Source database name |
| `experimental_label` | Label as reported by this source (raw, not normalized) |
| `method_universal` | Normalized method name |
| `raw_method` | Method string as found in source |
| `evidence_type` | structural / kinetic / staining_binding / literature_curated / computational |
| `evidence_weight` | 0.0–3.0 |
| `confidence` | 0–100 |
| `pdb_id` / `emdb_id` | Structural depositions from this observation |
| `resolution` | Resolution from this observation |
| `doi` / `pmid` | Publication from this observation |
| `category` | Source-provided category (e.g., "hexapeptide") |
| `is_consensus_winner` | `"true"` for all TSV-derived observations; `"false"` for conflict stubs |
| `data_completeness` | `"full"` for TSV rows; `"incomplete"` for conflict stubs |
| `graph_layer` | `"observation"` or `"conflict"` |

### What is NOT captured

**Data limitation:** This is a partial pre-consensus reconstruction, not a full audit log.

The pipeline discards within-group non-winners before writing to TSV. Within each dedup-key
group, only the highest-evidence-weight observation is exported. All other observations in
the same group are lost (not recoverable from any available output file without re-running
the pipeline with a custom export hook).

Concretely, if WALTZ-DB contributed 3 records for the same hexapeptide and 2 of them were
discarded in favor of the third (as consensus winner), those 2 discarded records have no
Observation nodes. The graph captures the inter-group provenance (different sources seeing
the same sequence) but not the intra-group provenance (multiple records within one source).

---

## 5. Conflict Stubs

### What conflicts.json records

`conflicts.json` is generated by `unifier.detect_conflicts()` which runs BEFORE consensus
resolution. It records all dedup groups where two or more source records carry different
`experimental_label` values — including cases that are later resolved by weighting AND
cases that are unresolvable.

**All 3 biological conflicts (P05453 hexapeptides) were resolved by weighting and ARE
present in the consensus TSV** (non-amyloid won with a sufficient margin). They have full
SequenceInstance nodes already in the consensus layer.

The 4th entry in conflicts.json (`('', 'CPAD-Structures')`) is an artifact — it uses the
fallback key `(sequence, source_db)` with an empty sequence. It is skipped.

### Conflict Observation stubs

6 Observation stubs (2 per real conflict) are added to the graph:

- These attach to SequenceInstance nodes that already exist in the consensus layer.
- `is_consensus_winner = "false"` — these represent the rejected label-variant observations.
- `data_completeness = "incomplete"` — only source_db and experimental_label are available
  from conflicts.json. Method, evidence_weight, doi, etc. are absent.
- They carry `graph_layer = "conflict"`.

**Why include them:** They document that an opposing label existed. A SequenceInstance
classified as `non-amyloid` (is_amyloid=False) in the consensus layer, but with a
`conflict`-layer Observation claiming `amyloid`, should be treated with extra skepticism.

---

## 6. Distinguishing the Two Layers in Queries

### Filter by layer attribute

All edges carry a `graph_layer` attribute:

```
graph_layer = "consensus"    → Layer 1 (SequenceInstance-centric)
graph_layer = "observation"  → Layer 2 (Observation-centric, TSV-derived)
graph_layer = "conflict"     → Layer 2 extension (conflict stubs)
```

### Typical queries by layer

**Consensus layer only** — "What is the consensus result for sequence X?"
→ Filter SequenceInstance nodes. Read `is_amyloid`, `confidence`, `evidence_type` directly.

**Observation layer only** — "How many independent sources saw this sequence?"
→ Count OBSERVED_AS edges pointing at a SequenceInstance.

**Cross-layer** — "Show sequences where consensus says amyloid but an observation says non-amyloid"
→ SequenceInstance.is_amyloid = True AND any connected Observation.experimental_label = "non-amyloid"

**Conflict detection** — "Which sequences have conflict-layer observations?"
→ Observation nodes WHERE graph_layer = "conflict" AND OBSERVED_AS → SequenceInstance

---

## 7. Mapping to TuringDB

TuringDB ingests `amyloid_graph_turingdb.jsonl` — a self-contained file containing both
nodes and edges in TuringDB JSONL format. Each line is one JSON object.

Node format:
```json
{"type": "node", "id": "<int>", "labels": ["<NodeType>"], "properties": {...}}
```

Relationship format:
```json
{"type": "relationship", "id": "<int>", "label": "<EDGE_TYPE>",
 "start": {"id": "<int>"}, "end": {"id": "<int>"}, "properties": {...}}
```

Integer node IDs are stable within a single file (assigned sequentially by
`convert_to_turingdb.py`). The original string node IDs are preserved as
`properties._original_id` on each node for traceability back to `graph_nodes.csv`.

The `properties.graph_layer` field on every edge enables layer-scoped queries in TuringDB
without needing separate graph namespaces.

Import command:
```
cp amyloid_graph_turingdb.jsonl ~/.turing/data/
turingdb
> LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db
```

To regenerate `amyloid_graph_turingdb.jsonl` from the CSV sources:
```
python convert_to_turingdb.py
```

Node type labels used in the graph (as-is — no renaming):

| Node type | Count | Primary display field |
|---|---|---|
| `SequenceInstance` | 4,348 | `name` = protein + region or sequence |
| `Observation` | 4,671 | `name` = source_db \| label \| method |
| `Protein` | 657 | `name` = protein_name |
| `SourceDatabase` | 8 | `name` = source_db |
| `ExperimentMethod` | 15 | `name` = method_universal |
| `Publication` | 627 | `name` = reference / DOI / PMID |
| `Structure` | 520 | `name` = pdb_id (resolution Å) |

Note: `graph.jsonl` is an older edge-only file generated by `graph_build.py` and uses string
node IDs. It is not suitable for direct TuringDB import. Use `amyloid_graph_turingdb.jsonl`
for all TuringDB workflows.

---

## 8. Limitations

### 8.1 Only consensus winners are available as full Observations

Within each dedup-key group, only the highest-evidence-weight record is exported to the TSV.
All other records in the same group are not recoverable without re-running the pipeline.

### 8.2 300 multi-source SequenceInstances, not all sources

For the 297 biological identities that appear in multiple TSV rows (from different dedup
key paths), each row becomes a distinct Observation node. However, this is still an
undercount of true multi-source coverage: sources whose records were merged under the
same dedup key (and the non-winner dropped) are invisible.

### 8.3 SequenceInstance count vs Observation count

4,348 SequenceInstances vs 4,671 Observations (including 6 conflict stubs):
- ~4,048 SequenceInstances have exactly 1 Observation (single-source, no dedup overlap)
- ~300 SequenceInstances have 2+ Observations (multi-source overlap captured)

### 8.4 Conflict stub completeness

6 Observation stubs from conflict entries have `data_completeness = "incomplete"`.
Method, evidence_weight, confidence, doi, and pmid are all empty strings for these nodes.
They are present to document label disagreement, not to add quantitative evidence.

### 8.5 ThT binding hub remains dominant

2,030 entries from CPAD all use `method_universal = "ThT binding"`. The ThT binding
ExperimentMethod node has approximately 2,000 SUPPORTED_BY edges in each layer
(~4,000 total). Degree-based graph analytics will be dominated by this node.

### 8.6 Structural classification fields are mostly "unknown"

72% of SequenceInstances have `secondary_structure_class = "unknown"`. All classification
was inferred from text fields by the pipeline parser — not from authoritative annotation.

---

## 9. Sequence-Level Uncertainty and Disagreement Metrics

Computed by `compute_metrics.py`. Written to:
- `sequence_metrics.csv` — one row per SequenceInstance, all metrics plus identity fields
- `graph_nodes.csv` — metrics merged into `attributes_json` of every `SequenceInstance` node

### 9.1 Notation

For a SequenceInstance `S`, let `O(S)` = all Observation nodes with `OBSERVED_AS → S`
(both `"observation"` and `"conflict"` graph layers).

Per observation `o ∈ O(S)`:

| Symbol | Source field | Notes |
|---|---|---|
| `label(o)` | `experimental_label` | Always non-empty in this dataset |
| `source(o)` | `source_db` | |
| `method(o)` | `method_universal` | Empty string for the 6 conflict stubs |
| `weight(o)` | `evidence_weight` | Float; 0.0 for conflict stubs (field is empty) |
| `doi(o)`, `pmid(o)` | `doi`, `pmid` | Either may be empty |

`L(S)` = `SequenceInstance.experimental_label` — the consensus label.

`W_win`   = Σ `weight(o)` for `o` where `label(o) = L(S)`  
`W_total` = Σ `weight(o)` for all `o ∈ O(S)`  
`W_runner_up` = second-largest cumulative label weight (0 if only one label exists)  

---

### 9.2 Metric Definitions

#### `n_observations`
```
n_observations = |O(S)|
```
Total Observation nodes linked to this SequenceInstance, including conflict stubs
(`data_completeness = "incomplete"`).

---

#### `n_sources`
```
n_sources = |{ source(o) : o ∈ O(S) }|
```
Distinct `source_db` values across all observations.

---

#### `n_methods`
```
n_methods = |{ method(o) : o ∈ O(S),  method(o) ≠ "" }|
```
Distinct non-empty `method_universal` values. Conflict stubs contribute 0 (empty method).

---

#### `n_labels`
```
n_labels = |{ label(o) : o ∈ O(S) }|
```
Distinct `experimental_label` values. `1` = all observations agree on label category.
`2+` = label disagreement exists.

---

#### `n_publications`
```
pub_id(o) = "pub:doi:<sha8(doi)>"   if doi(o) ≠ ""
          = "pub:pmid:<sha8(pmid)>" if pmid(o) ≠ "" and doi(o) = ""
          = ""                       otherwise

n_publications = |{ pub_id(o) : o ∈ O(S),  pub_id(o) ≠ "" }|
```
Distinct publications linked to this SequenceInstance through its observations.

---

#### `conflict_score`
```
n_win          = |{ o ∈ O(S) : label(o) = L(S) }|
conflict_score = 1 - (n_win / n_observations)
```
Count-based fraction of observations that disagree with the consensus label.  
Range `[0, 1]`. `0` = all observations support consensus. `>0` = label disagreement exists.

**Why count-based, not weight-based:** The only opposing observations in this dataset are
conflict stubs with `weight = 0`. A weight-based formula would report `0.0` for all entries
including genuine conflicts. Count-based correctly detects label disagreement regardless of
evidence strength, and is the appropriate primary signal for disagreement.

**Empirical distribution:** `0.0` for 4,345 SequenceInstances; `0.333` for the 3 P05453
conflict sequences (1 opposing stub out of 3 total observations each).

---

#### `consensus_support_weight`
```
consensus_support_weight = W_win = Σ weight(o)  for o where label(o) = L(S)
```
Total evidence weight supporting the consensus label. Reflects experimental strength:
structural methods contribute 3.0 per observation, kinetic 2.0, staining/binding 1.0,
literature-curated 0.5, computational 0.0.

---

#### `opposing_support_weight`
```
opposing_support_weight = W_total - W_win
                        = Σ weight(o)  for o where label(o) ≠ L(S)
```
Total evidence weight supporting all non-consensus labels.

**Dataset note:** `0.0` for all entries in the current database. The only opposing
observations (conflict stubs) have `weight = 0` because their evidence details were not
retained in `conflicts.json`. This metric is meaningful in principle and will be non-zero
if the pipeline is extended to export full pre-consensus observations.

---

#### `confidence_margin`
```
label_weights = { lbl : Σ weight(o) for o where label(o) = lbl }
sorted descending → [w_1st, w_2nd, ...]

if n_labels = 1:
    confidence_margin = w_1st        ← UNOPPOSED SUPPORT, not a competitive margin
else:
    confidence_margin = w_1st - w_2nd  ← COMPETITIVE MARGIN (winner minus runner-up)
```
**Interpretation caveat:** when `n_labels = 1` (4,345 of 4,348 entries), this value is
the full `consensus_support_weight` — it reflects how much evidence backs the consensus,
not how much more evidence the winner has over a competitor. No runner-up exists. Check
`n_labels` before interpreting this as a competitive quantity.

For the 3 conflict entries: `confidence_margin = 3.0` (competitive: `W_non-amyloid = 3.0`
minus `W_amyloid = 0.0`). The margin equals the winner's weight because the runner-up
has no weighted evidence.

---

#### `normalized_margin`
```
W_runner_up = label_weights sorted descending[1], or 0.0 if n_labels = 1

if W_total > 0:
    normalized_margin = (W_win - W_runner_up) / W_total
else:
    normalized_margin = 0.0
```
The competitive weight advantage of the consensus label as a fraction of all evidence weight.  
Range `[-1, 1]`. `1.0` = all weight is on the consensus label. `0.0` = tied. Negative = consensus
label is the loser by weight (possible if a conflict stub overrides weighting, but not observed here).

**Dataset note:** trivially `1.0` for all 4,348 entries in the current database. Because
`W_runner_up = 0` for every entry (no opposing observations have positive weight),
`normalized_margin = W_win / W_total = 1.0` always. This metric will differentiate entries
when full pre-consensus observations with real evidence weights are included.

---

#### `method_diversity`
```
O_m(S) = { o ∈ O(S) : method(o) ≠ "" }
n_m    = |O_m(S)|

if n_m ≤ 1:
    method_diversity = 0.0
else:
    count_per_method_m = |{ o ∈ O_m : method(o) = m }|
    D                  = Σ_m (count_per_method_m / n_m)²    ← Simpson's concentration
    method_diversity   = 1 - D                               ← Simpson's diversity
```
Range `[0, 1)`. `0` = all observations share one method (or only one observation with a method).
`0.5` = two observations using two different methods.
`→1` = many observations each using a different method (maximum diversity).

**Empirical distribution:** `0.0` for 4,073 entries (single method); `0.5` for 275 entries
(two observations, two distinct methods).

---

#### `provenance_richness`
```
provenance_richness = log2(1 + n_sources)
                    + log2(1 + n_methods)
                    + log2(1 + n_publications)
```
Log-scaled composite of three independent provenance dimensions. Log scaling gives diminishing
returns for each additional source/method/publication (the second source is worth less than the
first; the third less than the second).

Range `[0, ∞)`. Typical values in this dataset:

| Profile | Value |
|---|---|
| 1 source, 1 method, 0 publications | `log2(2) + log2(2) + 0 = 2.0` |
| 1 source, 1 method, 1 publication | `1 + 1 + 1 = 3.0` |
| 2 sources, 2 methods, 0 publications | `log2(3) + log2(3) + 0 ≈ 3.17` |
| 2 sources, 2 methods, 2 publications | `≈ 3.17 + log2(3) ≈ 4.75` |

**Observed range:** `2.0` (min) to `4.755` (max). Mean `2.34`.

---

### 9.3 Dataset-Specific Observations

The following metrics are **trivially constant** in the current database due to data
availability constraints. They are correctly computed and correctly defined — the
uniformity is a property of the data, not a bug.

| Metric | Observed value | Reason |
|---|---|---|
| `normalized_margin` | `1.0` for all entries | All opposing observations are conflict stubs with `weight = 0`; `W_runner_up = 0` always |
| `opposing_support_weight` | `0.0` for all entries | Same reason |
| `confidence_margin` (competitive case) | Equals `consensus_support_weight` | Same reason — no weighted runner-up |

The **informative metrics** for this dataset are:

| Metric | What it differentiates |
|---|---|
| `conflict_score` | 3 entries with label disagreement vs. 4,345 clean |
| `n_sources` | 271 entries confirmed by 2 databases vs. 4,077 single-source |
| `n_methods` | 275 entries with method diversity vs. 4,073 single-method |
| `method_diversity` | Degree of method spread (0 or 0.5 in this dataset) |
| `consensus_support_weight` | Absolute evidence strength (range 0.5–18.0) |
| `provenance_richness` | Composite quality score (range 2.0–4.755) |
| `n_publications` | 0–5 distinct publications per entry |

---

### 9.4 Output Files

**`sequence_metrics.csv`** columns:

| Column | Type | Description |
|---|---|---|
| `seq_id` | str | SequenceInstance node ID |
| `sequence` | str | Amino acid sequence |
| `uniprot_id` | str | UniProt accession (may be empty) |
| `region_start` | str | 1-based start position |
| `region_end` | str | 1-based end position |
| `is_amyloid` | str | Consensus classification |
| `experimental_label` | str | Consensus label string |
| `n_observations` | int | |
| `n_sources` | int | |
| `n_methods` | int | |
| `n_labels` | int | |
| `n_publications` | int | |
| `conflict_score` | float | [0, 1] |
| `consensus_support_weight` | float | [0, ∞) |
| `opposing_support_weight` | float | [0, ∞) |
| `confidence_margin` | float | [0, ∞); see n_labels caveat |
| `normalized_margin` | float | [-1, 1]; trivially 1.0 in this dataset |
| `method_diversity` | float | [0, 1) |
| `provenance_richness` | float | [0, ∞) |

**`graph_nodes.csv`** — all columns above (except `seq_id`, `sequence`, `uniprot_id`,
`region_start`, `region_end`) are merged into the `attributes_json` of every
`SequenceInstance` node. Other node types are unchanged.
