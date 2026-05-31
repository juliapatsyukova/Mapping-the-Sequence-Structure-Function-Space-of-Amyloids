# Schema Notes — Amyloid Database

Generated: 2026-04-09

---

## 1. Ambiguities

### 1.1 The biological identity key is ambiguous for ~26% of entries

`AmyloidEntry.get_dedup_key()` uses a three-tier fallback:

1. `(sequence, uniprot_id, region_start, region_end)` — preferred
2. `(sequence, protein_name, region_start, region_end)` — when uniprot_id is missing
3. `(sequence, source_db)` — last resort when no protein identifier is available

**Problem:** Tier 3 entries from the same source DB with the same sequence but from different
proteins are silently merged into one entry. This is biologically wrong — the same hexapeptide
can be amyloidogenic in one protein and non-amyloidogenic in another.

~1,216 of 4,585 consensus entries (26.5%) have no `uniprot_id`. Among these, many are
WALTZ-DB hexapeptides where the protein identifier was simply not populated from the source.

**Recommendation:** Backfill `uniprot_id` from the UniProt API using `protein_name` before
deduplication. The pipeline has `fetch_sequences_and_info()` for this purpose — it should be
run before the deduplication phase, not independently.

### 1.2 `experimental_label` is not an enum

Values observed in the current database include: `"amyloid"`, `"non-amyloid"`,
`"amyloid-forming"`, `"fibril"`, `"unclear"`, `"positive"`, `"negative"`, plus empty strings.
There is no validation or normalization of this field — it is stored as-is from each parser.

This means the consensus voting in `DatasetUnifier.build_consensus()` can compare `"amyloid"`
vs `"amyloid-forming"` and treat them as different labels, potentially creating false conflicts.

**Recommendation:** Normalize `experimental_label` to a controlled vocabulary
(`amyloid | non_amyloid | unclear`) in each parser before the labels are used for voting.

### 1.3 `resolution` is stored as free text, not a numeric

Values in the current database include `"2.1"`, `"3.40"`, `"3.4 Å"`, `"NA"`, empty string.
This makes range queries and sorting impossible without string parsing.

**Recommendation:** Parse to `REAL` at ingest time; store as `resolution_angstrom REAL`.

### 1.4 `protein_family` is assigned by keyword matching, not UniProt hierarchy

`config.PROTEIN_FAMILIES` uses simple substring matching against `protein_name + disease`.
This can produce incorrect or missing assignments:
- `"Belongs to the CsgA/CsgB family"` (raw UniProt description) does not match any family keyword.
- `"IAPP_HUMAN"` matches `"amyloid-beta"` family via the `"iapp"` keyword — correct, but the
  family name `"amyloid-beta"` is misleading for IAPP (islet amyloid polypeptide).

**Recommendation:** Replace with UniProt-based family classification using the UniProt REST API
or a curated mapping table.

---

## 2. Design Tradeoffs

### 2.1 Flat `entries` table vs. normalized schema

**Current pipeline:** A single flat `entries` table with all 40+ columns, plus two satellite
tables for computed features. This is simple to export and query but:

- Protein metadata (organism, family) is repeated for every region of the same protein.
- Method/evidence metadata is embedded per-row rather than referenced.
- Publications and structures have no independent existence.

**Proposed schema (this document):** Normalized into 10 tables. This adds join complexity but:
- Protein edits propagate automatically.
- Method confidence/weight parameters live in one place (`evidence_methods`).
- Publications can be augmented (e.g., add journal, year) without touching entry rows.

The tradeoff is verbosity in queries — the `ml_features` view restores the flat layout for
machine learning use.

### 2.2 Consensus winner selection discards provenance

`build_consensus()` picks the **single best entry** for each dedup group and discards all
others. In the normalized schema, `raw_observations` preserves all source records
(with `seq_id = NULL` for dropped losers), but in the current SQLite export, losers are
completely gone — only the winner row is written.

This means the current `amyloid_database.sqlite` cannot answer: "Was this sequence ever
reported as non-amyloid by any source?" without re-running the pipeline.

**Recommendation:** Export `raw_observations` alongside the consensus table so provenance
is auditable without re-running.

### 2.3 Structural classification is entirely inferred, not authoritative

All 11 structural classification levels are assigned by regex matching against `raw_method`
and `notes` text fields in `BaseParser._infer_all_structural_fields()`. This means:
- They reflect the description in the source, not a curated expert annotation.
- A source that says "cryo-EM fibril study" will get `secondary_structure_class = cross_beta`
  even if the structure turns out to be something unusual.
- The vast majority of entries are classified as `"unknown"` (3,308/4,585 = 72%).

**Recommendation:** Add a `classification_source` column to record whether the classification
came from inference, manual curation, or an authoritative external source (e.g., AmyloidAtlas
which is already treated as ground-truth `cross_beta`).

### 2.4 `structure_type` (legacy) and levels 1–11 coexist

The pipeline carries both the old single-level `structure_type` (`fibril/oligomer/crystal/...`)
and the new 11-level hierarchy. The two are not synchronized — `structure_type` is inferred
from `raw_method` text independently of level 1 (`secondary_structure_class`). They can
be inconsistent: an entry could have `structure_type = fibril` and
`secondary_structure_class = unknown`.

**Recommendation:** Once the 11-level classification is trusted, deprecate `structure_type`
and derive it from levels 1–11 via a view or computed column.

### 2.5 Conflict resolution is aggressive

The 1.5× weight threshold in `build_consensus()` means that if the "amyloid" label has
cumulative weight 1.4 and "non-amyloid" has weight 1.0, the entry is **dropped** (neither
label wins). This currently removes 4 entries from the database. This is a conservative
but scientifically defensible choice — unresolved conflicts should not silently become
one answer.

However, the dropped entries are not prominently exposed in the output (only `conflicts.json`).
Downstream users may not realize these entries exist.

---

## 3. Things That Are Messy in the Current Pipeline

### 3.1 `sequence_length` is redundant but computed separately

`sequence_length` in the `entries` table is computed as `len(entry.sequence) if entry.sequence else 0`
during the SQLite insert, separately from `PhysicochemicalFeatures.length` which is also
`len(seq)` after `clean_sequence()`. If `clean_sequence()` strips characters, these two
values will differ silently.

### 3.2 Dipeptide frequencies in the TSV are sparsely populated

The TSV output (`consensus_unified.tsv`) contains `comp_dp_AA`, `comp_dp_AQ`, `comp_dp_QA`
columns, which are the only top-25 dipeptides that had non-zero values in the actual dataset.
These are dataset-specific artifacts — any re-run with a different input would produce
different column names. The SQLite schema does not store dipeptide frequencies at all.

### 3.3 The `UNIQUE(record_id, source_db, sequence)` constraint in the current SQLite

The current constraint `UNIQUE(record_id, source_db, sequence)` in the `entries` table is
**not the biological deduplication key** — it is a technical insert-level guard. Two entries
from different `source_db` values with the same sequence would both be inserted even if they
represent the same biological region. The actual deduplication happens in Python before
SQLite export; the constraint only prevents double-insertion of the same record within one
source.

The proposed schema replaces this with the correct biological natural key:
`UNIQUE(sequence, uniprot_id, region_start, region_end)`.

### 3.4 `INSERT OR REPLACE` in the SQLite exporter silently clobbers rows

`cursor.execute('INSERT OR REPLACE INTO entries ...')` in `sqlite.py` uses SQLite's
REPLACE semantics: if the UNIQUE constraint fires, the old row is **deleted** and a new one
is **inserted** (with a new autoincrement ID). This invalidates any foreign keys from
`physicochemical_features` and `sequence_composition` that referenced the old row's `id`.
The feature rows use `INSERT OR REPLACE` too, so they get re-inserted, but only if the
`entry_id` is correctly captured from `cursor.lastrowid` after the REPLACE. If two entries
collide, the feature data for the first is orphaned.

### 3.5 `method_universal` and `evidence_type` can be empty strings

When no method rule matches in `map_method_to_universal()`, the function returns empty
strings rather than `None`. These appear as `""` in the database. The `evidence_type`
column is used for filtering (e.g., `WHERE evidence_type != "computational"`), so empty
strings pass through as non-computational, even for truly unknown methods.

### 3.6 PDB enrichment is best-effort and can silently fail

`run_pipeline.py` Phase 3c wraps the entire `enrich_entries()` call in a broad `except Exception`
that logs a warning and continues. A failed enrichment produces a database with fewer
structure annotations than intended, with no record of which entries were affected.

---

## 4. Recommendations (Not Implemented)

1. **Normalize `experimental_label`** to a controlled 3-value vocabulary before consensus voting.

2. **Run sequence fetching before deduplication**, not as an optional post-parse phase.
   This ensures UniProt IDs are available for all entries before the dedup key is computed.

3. **Store `raw_observations`** (all pre-consensus records) in the SQLite export so provenance
   is queryable without re-running the pipeline.

4. **Parse `resolution` to REAL** at ingest. Store as `resolution_angstrom REAL`.

5. **Add `classification_source` flag** to structural classification fields:
   `'inferred'` | `'authoritative'` | `'manual'` so downstream users know which annotations
   to trust.

6. **Replace the `UNIQUE(record_id, source_db, sequence)` constraint** with the biological
   key `UNIQUE(sequence, uniprot_id, region_start, region_end)`.

7. **Use `INSERT OR IGNORE` instead of `INSERT OR REPLACE`** for the feature tables to
   avoid silent row deletion on constraint violation.

8. **Export dipeptide frequencies as a sparse separate table**, not as TSV columns.
   The current approach creates different column sets on each pipeline run.
