# Amyloid-DB TuringDB — Practical Query Manual

**Graph:** `amyloid_db`  
**Nodes:** 10,846 · **Edges:** 30,431  
**Last export:** `amyloid_graph_turingdb.jsonl`

---

## 1. Graph Structure

### Node types

| Label | Count | What it represents |
|---|---|---|
| `SequenceInstance` | 4,348 | A unique biological sequence (or region) with its consensus amyloid classification and all structural metadata. The backbone of the graph. |
| `Observation` | 4,671 | One record from one source database — a raw experimental report before consensus merging. |
| `Protein` | 657 | A protein entity (UniProt accession + name + organism). One protein can have many sequence regions. |
| `Publication` | 627 | A publication (DOI / PMID / author-year reference string). |
| `Structure` | 520 | A solved 3-D structure (PDB or EMDB entry), with resolution and method. |
| `ExperimentMethod` | 15 | A normalised experimental method (e.g. Cryo-EM, ThT binding). |
| `SourceDatabase` | 8 | One of eight contributing databases (CPAD, WALTZ-DB, AmyPro, AmyloGraph, AmyloidAtlas, AmyloidExplorer, CPAD-Structures, Cross-Beta DB). |

### Relationships

| Edge | Direction | Meaning |
|---|---|---|
| `PART_OF` | `SequenceInstance → Protein` | This region belongs to this protein. Carries `region_start`, `region_end`. |
| `OBSERVED_IN` | `SequenceInstance → SourceDatabase` | The consensus sequence was observed in this database. Carries `experimental_label`. |
| `SUPPORTED_BY` | `SequenceInstance → ExperimentMethod` | The consensus classification is supported by this method. Carries `evidence_weight`, `confidence`. |
| `REPORTED_IN` | `SequenceInstance → Publication` | The sequence is reported in this publication. |
| `HAS_STRUCTURE` | `SequenceInstance → Structure` | A solved 3-D structure exists for this region. Carries `pdb_id`, `resolution`. |
| `OBSERVED_AS` | `Observation → SequenceInstance` | This raw observation maps to this consensus sequence. |
| `FROM_SOURCE` | `Observation → SourceDatabase` | This observation came from this database. |

### Two graph layers

- **Consensus layer** (`graph_layer = 'consensus'`): deduplicated, evidence-weighted nodes.  
  Always filter to this layer for clean analysis.
- **Observation layer** (`graph_layer = 'observation'`): one node per source record.  
  Use this to trace where each piece of evidence came from.
- **Conflict layer** (`graph_layer = 'conflict'`): 6 stub observations for unresolved contradictions.

### Key property notes

All properties are **strings** in TuringDB — including numbers and booleans.
- Numeric comparisons: `toFloat(n.confidence) > 70`
- Boolean: `n.is_amyloid = 'True'`  (not `true`)
- Empty values: `n.protein_name <> ''`
- `is_amyloid`: `'True'` / `'False'`
- `experimental_label`: `'amyloid'` / `'non-amyloid'`
- `evidence_type`: `'structural'` / `'staining_binding'` / `'kinetic'` / `'literature_curated'`
- `aggregate_type`: `'synthetic'` / `'in_vitro'` / `'recombinant'` / `'ex_vivo'` / `'unknown'`
- `structure_type`: `'fibril'` / `'crystal'` / `'aggregate'` / `'unknown'`

---

## 2. Quick Orientation Queries

### Q1 — Count nodes by type

```cypher
MATCH (n)
RETURN labels(n)[0] AS node_type, count(n) AS total
ORDER BY total DESC
```

**What it does:** Returns a table of node types and counts.  
**Expect:** 7 rows: SequenceInstance 4348, Observation 4671, Protein 657, etc.  
**Use this first** to confirm the graph loaded correctly.

---

### Q2 — Count edges by type

```cypher
MATCH ()-[r]->()
RETURN type(r) AS edge_type, count(r) AS total
ORDER BY total DESC
```

**Expect:** 7 edge types. SUPPORTED_BY (9,292) is the largest because both `SequenceInstance` and `Observation` nodes connect to methods.

---

### Q3 — List all source databases

```cypher
MATCH (db:SourceDatabase)
RETURN db.source_db AS database
ORDER BY database
```

**Expect:** 8 rows — AmyloGraph, AmyloidAtlas, AmyloidExplorer, AmyPro, CPAD, CPAD-Structures, Cross-Beta DB, WALTZ-DB.

---

### Q4 — List all experiment methods with their evidence weights

```cypher
MATCH (m:ExperimentMethod)
RETURN m.method_universal AS method,
       m.evidence_type AS evidence_type,
       m.evidence_weight AS weight
ORDER BY toFloat(m.evidence_weight) DESC
```

**What it does:** Shows the 15 normalised methods ranked by evidence weight.  
**Biologically useful:** Cryo-EM (weight 3.0) is the strongest; ThT binding (1.0) is weaker — understanding this lets you weight results accordingly.

---

## 3. Protein and Sequence Exploration

### Q5 — All sequence regions of a protein

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE p.protein_name CONTAINS 'TAU'
RETURN p.protein_name AS protein,
       s.region_start AS start,
       s.region_end   AS end,
       s.sequence     AS sequence,
       s.is_amyloid   AS amyloid,
       s.confidence   AS confidence,
       s.structure_type AS structure
ORDER BY toInteger(s.region_start)
```

**What it does:** Lists all classified sequence windows across the Tau protein, ordered by position.  
**Biologically useful:** Reveals which parts of Tau are amyloidogenic vs non-amyloid — directly relevant to Alzheimer's and frontotemporal dementia research.  
**Expect:** ~328 rows covering known Tau aggregation regions.  
**Adapt for other proteins:** Replace `'TAU'` with `'SYUA'` (α-synuclein), `'A4_HUMAN'` (APP/Aβ), `'PRIO_HUMAN'` (prion), `'FUS'`, `'TADBP'` (TDP-43).

---

### Q6 — Top 10 proteins by number of sequence regions

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE s.graph_layer = 'consensus'
RETURN p.protein_name AS protein,
       p.organism AS organism,
       count(s) AS n_regions,
       sum(CASE WHEN s.is_amyloid = 'True' THEN 1 ELSE 0 END) AS n_amyloid,
       sum(CASE WHEN s.is_amyloid = 'False' THEN 1 ELSE 0 END) AS n_non_amyloid
ORDER BY n_regions DESC
LIMIT 10
```

**What it does:** Ranks proteins by total experimental coverage, splitting amyloid vs non-amyloid counts.  
**Expect:** ERF3_YEAST (460 regions), TAU_HUMAN (328), APOA1_HUMAN (245), CSGA_ECOLI (211), TADBP_HUMAN (206) lead the list.  
**Biologically useful:** Shows where experimental effort has concentrated. ERF3_YEAST's dominance reflects its use as a yeast prion model.

---

### Q7 — Find high-confidence amyloid regions in a protein

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE p.protein_name CONTAINS 'SYUA'
  AND s.is_amyloid = 'True'
  AND toFloat(s.confidence) >= 70
RETURN s.name AS region,
       s.sequence AS sequence,
       s.confidence AS confidence,
       s.evidence_weight AS weight,
       s.secondary_structure_class AS structure_class,
       s.n_observations AS n_obs
ORDER BY toFloat(s.confidence) DESC
```

**What it does:** Finds the most reliably amyloid-positive regions of α-synuclein.  
**Biologically useful:** Pinpoints the NAC region (residues ~61–95) that drives Parkinson's-associated aggregation. High confidence + cross_beta class = strong structural evidence.

---

### Q8 — Sequences with conflicting evidence (amyloid vs non-amyloid disagreement)

```cypher
MATCH (s:SequenceInstance)
WHERE toFloat(s.conflict_score) > 0
RETURN s.name AS sequence_region,
       s.conflict_score AS conflict_score,
       s.confidence AS confidence,
       s.experimental_label AS consensus_label,
       s.n_observations AS n_obs,
       s.n_sources AS n_sources
ORDER BY toFloat(s.conflict_score) DESC
```

**What it does:** Finds sequences where different sources reported opposite labels.  
**Expect:** 3 sequences, all ERF3_YEAST regions with conflict_score 0.33.  
**Biologically useful:** These are scientifically contested regions worth examining manually. The `opposing_support_weight` property tells you how strong the minority evidence is.

---

## 4. Evidence Tracing

### Q9 — Full evidence trail for one sequence

```cypher
MATCH (s:SequenceInstance)
WHERE s.name = 'TAU_HUMAN [306–311]'
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
OPTIONAL MATCH (s)-[:SUPPORTED_BY]->(m:ExperimentMethod)
OPTIONAL MATCH (s)-[:OBSERVED_IN]->(db:SourceDatabase)
OPTIONAL MATCH (s)-[:REPORTED_IN]->(pub:Publication)
OPTIONAL MATCH (s)-[:HAS_STRUCTURE]->(st:Structure)
RETURN s.is_amyloid AS amyloid,
       s.confidence AS confidence,
       s.evidence_weight AS weight,
       p.protein_name AS protein,
       p.organism AS organism,
       m.method_universal AS method,
       db.source_db AS database,
       pub.reference AS publication,
       st.pdb_id AS pdb,
       st.resolution AS resolution_A
```

**What it does:** One-stop summary of all evidence associated with a single sequence region.  
**Biologically useful:** Gives the full provenance picture: which databases reported it, by what method, in which publication, with what structure.  
**Adapt:** Replace the `name` value with any label from Q5 output.

---

### Q10 — All sequences supported by Cryo-EM

```cypher
MATCH (s:SequenceInstance)-[:SUPPORTED_BY]->(m:ExperimentMethod)
WHERE m.method_universal = 'Cryo-EM'
  AND s.graph_layer = 'consensus'
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
RETURN p.protein_name AS protein,
       s.name AS region,
       s.is_amyloid AS amyloid,
       s.structure_type AS structure_type,
       s.secondary_structure_class AS ss_class,
       toFloat(s.confidence) AS confidence
ORDER BY confidence DESC
LIMIT 30
```

**What it does:** Returns the highest-confidence Cryo-EM supported regions.  
**Expect:** ~902 matches. Filtering to `confidence >= 70` reduces to the most reliable subset.  
**Biologically useful:** Cryo-EM has the highest evidence weight (3.0) — these sequences have the strongest direct structural evidence for fibril architecture.

---

### Q11 — Sequences with multiple independent methods (method diversity)

```cypher
MATCH (s:SequenceInstance)
WHERE s.graph_layer = 'consensus'
  AND toFloat(s.method_diversity) > 0
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
RETURN s.name AS region,
       p.protein_name AS protein,
       s.n_methods AS n_methods,
       s.method_diversity AS method_diversity,
       s.confidence AS confidence,
       s.is_amyloid AS amyloid
ORDER BY toFloat(s.method_diversity) DESC
LIMIT 20
```

**What it does:** Finds sequences characterised by more than one type of experiment.  
**Biologically useful:** Multi-method sequences are the most reliable entries in the database. `method_diversity > 0` means at least two different method types were used.

---

### Q12 — Trace Observations back to raw source records

```cypher
MATCH (obs:Observation)-[:OBSERVED_AS]->(s:SequenceInstance)
WHERE s.name CONTAINS 'TAU_HUMAN'
OPTIONAL MATCH (obs)-[:FROM_SOURCE]->(db:SourceDatabase)
OPTIONAL MATCH (obs)-[:SUPPORTED_BY]->(m:ExperimentMethod)
RETURN obs.record_id AS source_record_id,
       db.source_db AS database,
       obs.experimental_label AS raw_label,
       m.method_universal AS method,
       obs.raw_method AS raw_method_string,
       obs.evidence_weight AS weight,
       s.name AS maps_to_sequence
ORDER BY database
```

**What it does:** For each Tau region, shows the raw per-database observation that contributed to the consensus, including the original record ID and verbatim method description.  
**Biologically useful:** Lets you trace any consensus result back to its original source record for manual verification.

---

## 5. Cross-Database Comparison

### Q13 — How many sequences does each database contribute?

```cypher
MATCH (s:SequenceInstance)-[:OBSERVED_IN]->(db:SourceDatabase)
WHERE s.graph_layer = 'consensus'
RETURN db.source_db AS database,
       count(s) AS n_sequences,
       sum(CASE WHEN s.is_amyloid = 'True' THEN 1 ELSE 0 END) AS amyloid_positive,
       sum(CASE WHEN s.is_amyloid = 'False' THEN 1 ELSE 0 END) AS non_amyloid,
       avg(toFloat(s.confidence)) AS avg_confidence
ORDER BY n_sequences DESC
```

**What it does:** Compares the contribution of each database — size, label balance, and average confidence.  
**Biologically useful:** WALTZ-DB and CPAD contribute large numbers of synthetic peptides; CPAD-Structures and AmyloidAtlas contribute structural data. Knowing this helps interpret results.

---

### Q14 — Sequences reported by two or more databases

```cypher
MATCH (s:SequenceInstance)-[:OBSERVED_IN]->(db:SourceDatabase)
WHERE s.graph_layer = 'consensus'
WITH s, collect(db.source_db) AS databases, count(db) AS n_db
WHERE n_db >= 2
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
RETURN s.name AS region,
       p.protein_name AS protein,
       databases AS source_databases,
       s.is_amyloid AS amyloid,
       s.confidence AS confidence
ORDER BY n_db DESC
LIMIT 20
```

**What it does:** Finds sequences with multi-database agreement.  
**Biologically useful:** Sequences confirmed in two or more independent databases have higher evidential credibility — important for selecting benchmark cases.

---

### Q15 — Label agreement across databases for one protein

```cypher
MATCH (obs:Observation)-[:OBSERVED_AS]->(s:SequenceInstance)-[:PART_OF]->(p:Protein)
MATCH (obs)-[:FROM_SOURCE]->(db:SourceDatabase)
WHERE p.protein_name CONTAINS 'B2MG'
RETURN s.name AS region,
       db.source_db AS database,
       obs.experimental_label AS label,
       obs.method_universal AS method
ORDER BY s.region_start, database
```

**What it does:** Shows how different databases labelled each region of β-2-microglobulin.  
**Biologically useful:** β-2-microglobulin (B2MG_HUMAN) causes dialysis-related amyloidosis — cross-database label comparison reveals where databases agree and where they diverge on this clinically important protein.

---

## 6. Structure-Focused Queries

### Q16 — All structures with high-resolution data

```cypher
MATCH (s:SequenceInstance)-[r:HAS_STRUCTURE]->(st:Structure)
WHERE st.pdb_id <> ''
  AND toFloat(st.resolution) <= 4
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
RETURN p.protein_name AS protein,
       s.name AS region,
       st.pdb_id AS pdb_id,
       st.resolution AS resolution_A,
       st.method_universal AS structure_method,
       s.secondary_structure_class AS ss_class,
       s.structure_type AS structure_type
ORDER BY toFloat(st.resolution)
LIMIT 25
```

**What it does:** Finds sequence regions with atomically resolved structures (≤4 Å).  
**Biologically useful:** High-resolution structures (Cryo-EM, microED) give direct evidence for cross-beta architecture. Resolution ≤ 4 Å typically allows side-chain modelling.

---

### Q17 — Protein with the most solved structures

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
MATCH (s)-[:HAS_STRUCTURE]->(st:Structure)
RETURN p.protein_name AS protein,
       count(DISTINCT st) AS n_structures,
       collect(DISTINCT st.pdb_id)[..5] AS sample_pdb_ids,
       avg(toFloat(st.resolution)) AS avg_resolution_A
ORDER BY n_structures DESC
LIMIT 10
```

**What it does:** Ranks proteins by structural coverage in the PDB.  
**Biologically useful:** Highly structured proteins (Tau, α-synuclein, Aβ) have many polymorphic fibril structures — this query reveals the structural diversity landscape per protein.

---

### Q18 — Explore structural polymorphs (multiple structures for one protein region)

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
MATCH (s)-[:HAS_STRUCTURE]->(st:Structure)
WHERE p.protein_name CONTAINS 'TAU'
RETURN s.name AS region,
       collect(st.pdb_id + ' (' + st.resolution + 'Å)') AS structures,
       count(st) AS n_structures,
       s.secondary_structure_class AS ss_class
ORDER BY n_structures DESC
LIMIT 15
```

**What it does:** For Tau, finds regions that have multiple distinct PDB entries — evidence for polymorphism.  
**Biologically useful:** Tau fibril polymorphs are a key feature of different tauopathies (Alzheimer's vs CBD vs PSP). Multiple structures on the same region = polymorph candidates.  
**Adapt:** Replace `'TAU'` with `'SYUA'` for α-synuclein polymorphs relevant to Parkinson's vs MSA.

---

## 7. Publication-Centric Exploration

### Q19 — Publications with the most associated sequences

```cypher
MATCH (s:SequenceInstance)-[:REPORTED_IN]->(pub:Publication)
WHERE pub.reference <> ''
RETURN pub.reference AS publication,
       pub.doi AS doi,
       count(s) AS n_sequences,
       sum(CASE WHEN s.is_amyloid = 'True' THEN 1 ELSE 0 END) AS amyloid,
       sum(CASE WHEN s.is_amyloid = 'False' THEN 1 ELSE 0 END) AS non_amyloid
ORDER BY n_sequences DESC
LIMIT 10
```

**What it does:** Finds the most data-rich publications in the graph.  
**Expect:** Lövestam et al. 2022 (46 sequences), Lövestam et al. 2024 (40), Wilkinson et al. 2023 (30).  
**Biologically useful:** Large systematic structure papers (often Cryo-EM Tau/α-synuclein studies) dominate. Useful for literature-guided graph navigation.

---

### Q20 — All sequences from one publication

```cypher
MATCH (s:SequenceInstance)-[:REPORTED_IN]->(pub:Publication)
WHERE pub.reference CONTAINS 'Lövestam'
  AND pub.reference CONTAINS '2022'
OPTIONAL MATCH (s)-[:PART_OF]->(p:Protein)
OPTIONAL MATCH (s)-[:HAS_STRUCTURE]->(st:Structure)
RETURN p.protein_name AS protein,
       s.name AS region,
       s.is_amyloid AS amyloid,
       s.confidence AS confidence,
       st.pdb_id AS pdb_id,
       st.resolution AS resolution_A
ORDER BY p.protein_name, toInteger(s.region_start)
```

**What it does:** Reconstructs the complete sequence/structure dataset from a specific paper.  
**Biologically useful:** Lövestam et al. 2022 is a major Tau polymorph Cryo-EM paper. Querying by publication lets you audit what the graph captured from a paper you know well.

---

## 8. Advanced Graph Analysis

### Q21 — Proteins connected via shared publications (protein co-publication network)

```cypher
MATCH (p1:Protein)<-[:PART_OF]-(s1:SequenceInstance)-[:REPORTED_IN]->(pub:Publication)
    <-[:REPORTED_IN]-(s2:SequenceInstance)-[:PART_OF]->(p2:Protein)
WHERE p1.protein_name < p2.protein_name
RETURN p1.protein_name AS protein_A,
       p2.protein_name AS protein_B,
       count(DISTINCT pub) AS shared_publications,
       collect(DISTINCT pub.reference)[..3] AS sample_publications
ORDER BY shared_publications DESC
LIMIT 15
```

**What it does:** Finds pairs of proteins that appear together in the same publications.  
**Biologically useful:** Co-published proteins are often functionally related or compared experimentally. Tau + α-synuclein co-appearing in publications signals cross-disease comparisons.

---

### Q22 — Sequences with the strongest evidence (composite provenance score)

```cypher
MATCH (s:SequenceInstance)
WHERE s.graph_layer = 'consensus'
  AND s.is_amyloid = 'True'
  AND s.protein_name <> ''
RETURN s.name AS region,
       s.protein_name AS protein,
       toFloat(s.provenance_richness) AS provenance_score,
       toFloat(s.confidence) AS confidence,
       s.n_observations AS observations,
       s.n_sources AS databases,
       s.n_methods AS methods,
       s.n_publications AS publications,
       s.secondary_structure_class AS ss_class
ORDER BY provenance_score DESC, confidence DESC
LIMIT 20
```

**What it does:** Ranks amyloid-positive sequences by `provenance_richness` — a composite score combining evidence diversity.  
**Biologically useful:** The top sequences represent the best-characterised amyloid regions across the entire database. Use this as a shortlist for benchmark studies or structural modelling.

---

### Q23 — Graph neighbourhood of α-synuclein (2-hop expansion)

```cypher
MATCH path = (p:Protein)-[:PART_OF*..1]-(s:SequenceInstance)-[*1..2]-(n)
WHERE p.protein_name = 'SYUA_HUMAN'
  AND s.is_amyloid = 'True'
RETURN p, s, n
LIMIT 50
```

**What it does:** Returns the local 2-hop neighbourhood around amyloid-positive α-synuclein regions — sequences + their methods, databases, structures, and publications.  
**TuringDB tip:** Run this in the visual explorer to get an interactive subgraph of the α-synuclein amyloid landscape. Keep `LIMIT 50` to avoid overloading the canvas.

---

### Q24 — Find the most central sequences (highest degree / most connections)

```cypher
MATCH (s:SequenceInstance)
WHERE s.graph_layer = 'consensus'
OPTIONAL MATCH (s)-[r]->()
WITH s, count(r) AS out_degree
OPTIONAL MATCH ()-[r2]->(s)
WITH s, out_degree, count(r2) AS in_degree
RETURN s.name AS region,
       s.protein_name AS protein,
       out_degree + in_degree AS total_degree,
       out_degree,
       in_degree,
       s.is_amyloid AS amyloid
ORDER BY total_degree DESC
LIMIT 15
```

**What it does:** Finds the most connected SequenceInstance nodes — those linked to the most methods, databases, publications, and structures simultaneously.  
**Biologically useful:** High-degree sequences are the most experimentally characterised; they are good candidates for detailed biological review or use as positive controls.

---

### Q25 — Shortest path between two proteins via shared evidence

```cypher
MATCH path = shortestPath(
  (p1:Protein)-[*..6]-(p2:Protein)
)
WHERE p1.protein_name = 'TAU_HUMAN'
  AND p2.protein_name = 'SYUA_HUMAN'
RETURN path
```

**What it does:** Finds the shortest connection between Tau and α-synuclein through the graph.  
**Biologically useful:** Tau and α-synuclein co-aggregate in some mixed proteinopathies. A short path via shared publications, methods, or databases signals experimental cross-links in the literature.  
**TuringDB tip:** This renders as a path in the visual explorer. Click each intermediate node to inspect its properties.

---

### Q26 — Disease-annotated amyloid sequences

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE s.disease <> ''
  AND s.is_amyloid = 'True'
RETURN s.disease AS disease,
       p.protein_name AS protein,
       count(s) AS n_regions,
       avg(toFloat(s.confidence)) AS avg_confidence,
       collect(DISTINCT s.structure_type) AS structure_types
ORDER BY n_regions DESC
LIMIT 20
```

**What it does:** Aggregates disease annotations across amyloid-positive sequences.  
**Biologically useful:** Directly maps proteins to diseases in the graph. Useful for identifying which proteins have the most disease-associated amyloid regions and what structural types they form.

---

### Q27 — Sequences with mutation annotations

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE s.mutation <> ''
  AND s.graph_layer = 'consensus'
RETURN p.protein_name AS protein,
       s.name AS region,
       s.mutation AS mutation,
       s.is_amyloid AS amyloid,
       s.pathogenicity AS pathogenicity,
       s.confidence AS confidence
ORDER BY p.protein_name
```

**What it does:** Lists all sequence regions with recorded mutations.  
**Biologically useful:** Pathogenic mutations (e.g. in SYUA, PRIO) have well-established amyloidogenicity. `pathogenicity = 'pathogenic'` identifies clinically relevant variants.

---

### Q28 — Ex-vivo sequences (patient-derived material)

```cypher
MATCH (s:SequenceInstance)-[:PART_OF]->(p:Protein)
WHERE s.aggregate_type = 'ex_vivo'
  AND s.graph_layer = 'consensus'
RETURN p.protein_name AS protein,
       p.organism AS organism,
       s.name AS region,
       s.is_amyloid AS amyloid,
       s.disease AS disease,
       s.confidence AS confidence,
       s.n_observations AS n_obs
ORDER BY p.protein_name
```

**What it does:** Filters to sequences derived from biological specimens (patient tissue, biopsy, etc.), as opposed to synthetic or recombinant peptides.  
**Biologically useful:** Ex-vivo sequences have the highest translational relevance — they represent real disease material, not laboratory constructs.  
**Expect:** ~136 sequences. These are the most clinically grounded entries.

---

## 9. TuringDB Visual Exploration Tips

### Starting points for visual exploration

1. **Start with a protein, not the whole graph.**  
   Run Q5 (`MATCH (s)-[:PART_OF]->(p) WHERE p.protein_name CONTAINS '...'`) to get a manageable protein-centric subgraph. Then click sequence nodes to expand their evidence connections.

2. **Use `LIMIT` aggressively.**  
   The full graph has 10,846 nodes. Never run `MATCH (n) RETURN n` without a limit.  
   Start with `LIMIT 25`, expand to 50–100 only after confirming what you're looking at.

3. **Filter to the consensus layer first.**  
   Add `WHERE n.graph_layer = 'consensus'` to any SequenceInstance query to exclude Observation nodes and halve the result size before you visualise.

4. **Pin anchor nodes.**  
   In TuringDB's visual canvas, pin a central node (e.g. a Protein or Publication) and expand its immediate neighbours. This prevents the graph from re-layouting when you add more nodes.

---

### Inspecting node properties

When you click a node in TuringDB, the **property panel** shows all stored fields. Key fields to look for:

| Node type | Most informative properties to check |
|---|---|
| SequenceInstance | `sequence`, `protein_name`, `is_amyloid`, `confidence`, `evidence_weight`, `n_observations`, `provenance_richness`, `secondary_structure_class`, `disease`, `mutation` |
| Protein | `protein_name`, `organism`, `protein_family`, `uniprot_id` |
| Structure | `pdb_id`, `resolution`, `method_universal` |
| Publication | `reference`, `doi`, `pmid` |
| Observation | `source_db`, `experimental_label`, `raw_method`, `record_id` |
| ExperimentMethod | `method_universal`, `evidence_weight`, `evidence_type` |

The `_original_id` field on each node (e.g. `seq:00042fb17fd6a9f6`, `protein:TAU_HUMAN`) is your traceability link back to the source CSV if you ever need to cross-check data.

---

### Avoiding overwhelming expansions

- **Expand one hop at a time.** Click a SequenceInstance → inspect its Protein first. Then expand to methods, then to publications. Never expand all neighbours at once on a hub node like WALTZ-DB (which connects to thousands of sequences).
- **Filter before visualising.** Use `WHERE` clauses in Cypher before rendering. Visual expansion ("expand all") from a Protein node like ERF3_YEAST (460 regions) will flood the canvas.
- **Collapse observation nodes.** Observation nodes outnumber SequenceInstance nodes. If you don't need per-source details, use `WHERE ... graph_layer = 'consensus'` to suppress them.
- **Use `DISTINCT` on large path queries.** Paths through shared publications can explode combinatorially. Always add `DISTINCT` or `LIMIT` when collecting across multiple `MATCH` hops.

---

### Useful filtering combinations

```cypher
-- Only well-evidenced amyloid-positive fibril sequences
WHERE s.is_amyloid = 'True'
  AND s.structure_type = 'fibril'
  AND toFloat(s.confidence) >= 70
  AND s.secondary_structure_class = 'cross_beta'

-- Only patient-relevant sequences
WHERE s.aggregate_type = 'ex_vivo'
  AND s.pathogenicity = 'pathogenic'

-- Only high-resolution structural entries
WHERE st.pdb_id <> ''
  AND toFloat(st.resolution) <= 3

-- Only observations with full provenance
WHERE obs.data_completeness = 'full'
  AND obs.is_consensus_winner = 'true'
```

---

## 10. Quick Reference Card

```
Node types:        SequenceInstance · Protein · Observation ·
                   Publication · Structure · ExperimentMethod · SourceDatabase

Relationships:     (Seq)-[:PART_OF]->(Protein)
                   (Seq)-[:SUPPORTED_BY]->(ExperimentMethod)
                   (Seq)-[:OBSERVED_IN]->(SourceDatabase)
                   (Seq)-[:REPORTED_IN]->(Publication)
                   (Seq)-[:HAS_STRUCTURE]->(Structure)
                   (Obs)-[:OBSERVED_AS]->(Seq)
                   (Obs)-[:FROM_SOURCE]->(SourceDatabase)

is_amyloid:        'True' / 'False'       (string)
experimental_label:'amyloid' / 'non-amyloid'
evidence_type:     'structural' · 'staining_binding' · 'kinetic' · 'literature_curated'
aggregate_type:    'synthetic' · 'in_vitro' · 'recombinant' · 'ex_vivo' · 'unknown'
structure_type:    'fibril' · 'crystal' · 'aggregate' · 'unknown'
graph_layer:       'consensus' · 'observation' · 'conflict'

Numeric casts:     toFloat(n.confidence) · toFloat(n.evidence_weight)
                   toInteger(n.region_start) · toFloat(n.provenance_richness)

Key proteins:      ERF3_YEAST (460 regions) · TAU_HUMAN (328) · APOA1_HUMAN (245)
                   TADBP_HUMAN (206) · A4_HUMAN (147) · B2MG_HUMAN (131)
                   FUS_HUMAN (102) · SYUA_HUMAN (93) · PRIO_HUMAN (87)

Key methods:       ThT binding (4,060) · Electron microscopy (2,390)
                   Cryo-EM (902) · Aggregation kinetics (518) · NMR solid-state (105)

Databases:         CPAD · WALTZ-DB · AmyPro · AmyloGraph
                   AmyloidAtlas · AmyloidExplorer · CPAD-Structures · Cross-Beta DB
```
