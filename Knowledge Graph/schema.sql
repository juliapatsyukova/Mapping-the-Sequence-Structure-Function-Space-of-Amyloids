-- =============================================================================
-- Amyloid Database Schema
-- Reverse-engineered from pipeline codebase and existing SQLite database.
-- Compatible with: SQLite 3.x, PostgreSQL 14+
-- Generated: 2026-04-09
-- =============================================================================
-- Notes on compatibility:
--   - AUTOINCREMENT is SQLite syntax; use SERIAL or GENERATED ALWAYS AS IDENTITY in PostgreSQL
--   - BOOLEAN is stored as INTEGER (0/1) in SQLite; native BOOLEAN in PostgreSQL
--   - JSON arrays in label_conflicts.labels/sources: TEXT in SQLite, JSON/JSONB in PostgreSQL
--   - Inline comments with -- are valid in both dialects
-- =============================================================================


-- =============================================================================
-- LAYER 1: REFERENCE / LOOKUP TABLES
-- =============================================================================

CREATE TABLE IF NOT EXISTS source_databases (
    source_db       TEXT        NOT NULL,
    entry_count     INTEGER,                     -- informational; not enforced
    description     TEXT,
    --
    PRIMARY KEY (source_db)
);

-- Populated from config.METHODS_RULES
CREATE TABLE IF NOT EXISTS evidence_methods (
    method_id       INTEGER     NOT NULL,
    method_universal TEXT       NOT NULL,        -- e.g. "Cryo-EM", "ThT binding"
    evidence_type   TEXT        NOT NULL,        -- structural | kinetic | staining_binding | literature_curated | computational
    evidence_weight REAL        NOT NULL,        -- 3.0 | 2.0 | 1.0 | 0.5 | 0.0
    tier            INTEGER     NOT NULL,        -- 0-5 (config.TIER_CONFIDENCE)
    confidence      INTEGER     NOT NULL,        -- 10|35|50|70|80|90
    pattern         TEXT,                        -- regex pattern from METHODS_RULES (for reference)
    --
    PRIMARY KEY (method_id),
    UNIQUE (method_universal),
    CHECK (evidence_type IN ('structural','kinetic','staining_binding','literature_curated','computational')),
    CHECK (evidence_weight >= 0.0),
    CHECK (confidence BETWEEN 0 AND 100),
    CHECK (tier BETWEEN 0 AND 5)
);

-- Pre-populated values from config.py
INSERT OR IGNORE INTO evidence_methods VALUES (1,  'XRD (cross-β diffraction)',              'structural',          3.0, 4, 90, '\bxrd\b|\bx[-\s]?ray diffraction\b|\bfiber diffraction\b');
INSERT OR IGNORE INTO evidence_methods VALUES (2,  'Cryo-EM',                                'structural',          3.0, 3, 70, '\bcryo[-\s]?em\b');
INSERT OR IGNORE INTO evidence_methods VALUES (3,  'Electron microscopy (TEM/EM)',            'structural',          3.0, 3, 70, '\btem\b|\btransmission electron microscopy\b|\belectron microscopy\b');
INSERT OR IGNORE INTO evidence_methods VALUES (4,  'AFM',                                    'structural',          3.0, 3, 70, '\bafm\b|\batomic force microscopy\b');
INSERT OR IGNORE INTO evidence_methods VALUES (5,  'NMR (solid-state)',                      'structural',          3.0, 3, 70, '\bsolid[-\s]?state nmr\b|\bssnmr\b');
INSERT OR IGNORE INTO evidence_methods VALUES (6,  'NMR (solution)',                         'kinetic',             2.0, 2, 50, '\bsolution nmr\b');
INSERT OR IGNORE INTO evidence_methods VALUES (7,  'NMR (unspecified)',                      'structural',          3.0, 3, 70, '\bnmr\b');
INSERT OR IGNORE INTO evidence_methods VALUES (8,  'X-ray (structural)',                     'structural',          3.0, 3, 70, '\bcrystallograph\b|\bx-ray\b|\bxray\b');
INSERT OR IGNORE INTO evidence_methods VALUES (9,  'Electron crystallography / microED',     'structural',          3.0, 3, 70, '\bmicroed\b|\belectron crystallography\b');
INSERT OR IGNORE INTO evidence_methods VALUES (10, 'FTIR',                                   'staining_binding',    1.0, 5, 80, '\bftir\b|\binfrared\b');
INSERT OR IGNORE INTO evidence_methods VALUES (11, 'CD',                                     'staining_binding',    1.0, 5, 80, '\bcircular dichroism\b|\bcd\b');
INSERT OR IGNORE INTO evidence_methods VALUES (12, 'Raman',                                  'staining_binding',    1.0, 5, 80, '\braman\b');
INSERT OR IGNORE INTO evidence_methods VALUES (13, 'ThT binding',                            'staining_binding',    1.0, 1, 35, '\bthioflavin[-\s]?t\b|\btht\b|\bth[-\s]?t\b');
INSERT OR IGNORE INTO evidence_methods VALUES (14, 'Congo Red binding',                      'staining_binding',    1.0, 1, 35, '\bcongo[-\s]?red\b');
INSERT OR IGNORE INTO evidence_methods VALUES (15, 'Proteostat binding',                     'staining_binding',    1.0, 1, 35, '\bproteostat\b');
INSERT OR IGNORE INTO evidence_methods VALUES (16, 'Aggregation kinetics',                   'kinetic',             2.0, 2, 50, '\bkinetic\b|\blag\b|\baggregation rate\b|\bseeding\b|\baggregation assay\b');
INSERT OR IGNORE INTO evidence_methods VALUES (17, 'PASTA 2.0 (prediction)',                 'computational',       0.0, 0, 10, '\bpasta\b');
INSERT OR IGNORE INTO evidence_methods VALUES (18, 'AGGRESCAN (prediction)',                 'computational',       0.0, 0, 10, '\baggrescan\b');
INSERT OR IGNORE INTO evidence_methods VALUES (19, 'TANGO (prediction)',                     'computational',       0.0, 0, 10, '\btango\b');
INSERT OR IGNORE INTO evidence_methods VALUES (20, 'WALTZ (prediction)',                     'computational',       0.0, 0, 10, '\bwaltz\b');
INSERT OR IGNORE INTO evidence_methods VALUES (21, 'Literature-curated',                     'literature_curated',  0.5, 0, 10, '\bliterature[-\s]?curated\b|\bcurated\b');


-- =============================================================================
-- LAYER 2: BIOLOGICAL IDENTITY TABLES
-- =============================================================================

CREATE TABLE IF NOT EXISTS proteins (
    uniprot_id      TEXT        NOT NULL,        -- e.g. "P05067", "Q13148"
    protein_name    TEXT,                        -- as reported by source (not normalized)
    organism        TEXT,                        -- e.g. "Homo sapiens"
    protein_family  TEXT,                        -- from config.PROTEIN_FAMILIES keyword match
    --
    PRIMARY KEY (uniprot_id)
);

CREATE INDEX IF NOT EXISTS idx_proteins_family ON proteins (protein_family);
CREATE INDEX IF NOT EXISTS idx_proteins_organism ON proteins (organism);


CREATE TABLE IF NOT EXISTS sequence_instances (
    seq_id                      INTEGER     NOT NULL,
    -- Biological identity (natural key)
    sequence                    TEXT        NOT NULL,   -- cleaned AA sequence
    uniprot_id                  TEXT,                   -- NULL for ~1216/4585 entries
    region_start                INTEGER,                -- 1-based; NULL when not available
    region_end                  INTEGER,                -- 1-based; NULL when not available
    sequence_length             INTEGER     NOT NULL,   -- len(sequence)
    -- Consensus classification
    is_amyloid                  BOOLEAN     NOT NULL,   -- consensus result
    experimental_label          TEXT,                   -- winning label string
    category                    TEXT,                   -- source category (e.g. "hexapeptide")
    -- Evidence summary (from winning observation)
    evidence_weight             REAL        NOT NULL DEFAULT 0.0,
    confidence                  INTEGER     NOT NULL DEFAULT 0,
    evidence_type               TEXT,
    method_universal            TEXT,
    raw_method                  TEXT,
    -- Biological context
    disease                     TEXT,
    tissue                      TEXT,
    mutation                    TEXT,
    notes                       TEXT,
    -- Structural IDs (from winning observation or enrichment)
    pdb_id                      TEXT,
    emdb_id                     TEXT,
    resolution                  TEXT,
    -- =====================================================================
    -- Multi-dimensional structural classification (11 levels)
    -- All values are inferred by BaseParser._infer_all_structural_fields()
    -- =====================================================================
    -- Level 1
    secondary_structure_class   TEXT        NOT NULL DEFAULT 'unknown',  -- cross_beta | cross_alpha | unknown
    -- Level 2
    interface_type              TEXT        NOT NULL DEFAULT 'unknown',  -- steric_zipper | LARKS | cross_alpha_sheet | unknown
    -- Level 3
    strand_arrangement          TEXT        NOT NULL DEFAULT 'unknown',  -- parallel_in_register | antiparallel | out_of_register | mixed | unknown
    -- Level 4
    fold_topology               TEXT        NOT NULL DEFAULT 'unknown',  -- beta_arcade | greek_key | beta_solenoid | beta_helix | beta_sandwich | superpleated | other | unknown
    -- Level 5
    fold_shape                  TEXT        NOT NULL DEFAULT 'unknown',  -- J_shaped | S_shaped | C_shaped | LS_shaped | kidney | elongated | compact | other | unknown
    -- Level 6
    zipper_class                TEXT,                                    -- "1"-"8" or "1"-"15"; empty string = unknown
    -- Level 7
    protofilament_symmetry      TEXT        NOT NULL DEFAULT 'unknown',  -- C1 | C2 | C3 | screw_21 | unknown
    -- Level 8
    protofilament_count         INTEGER     NOT NULL DEFAULT 0,          -- 0 = unknown
    -- Level 9
    twist_handedness            TEXT        NOT NULL DEFAULT 'unknown',  -- left | right | unknown
    -- Level 10
    aggregate_type              TEXT        NOT NULL DEFAULT 'unknown',  -- ex_vivo | in_vitro | recombinant | synthetic | unknown
    -- Level 11
    polymorph_name              TEXT,                                    -- e.g. "PHF", "rod", "twister"; empty = unknown
    -- Additional classification
    pathogenicity               TEXT        NOT NULL DEFAULT 'unknown',  -- pathogenic | functional | non_pathogenic | unknown
    structure_type              TEXT        NOT NULL DEFAULT 'unknown',  -- fibril | oligomer | crystal | aggregate | monomer | unknown (legacy)
    --
    PRIMARY KEY (seq_id),
    UNIQUE (sequence, uniprot_id, region_start, region_end),
    FOREIGN KEY (uniprot_id) REFERENCES proteins (uniprot_id),
    CHECK (is_amyloid IN (0, 1)),
    CHECK (confidence BETWEEN 0 AND 100),
    CHECK (evidence_weight >= 0.0),
    CHECK (secondary_structure_class IN ('cross_beta','cross_alpha','unknown')),
    CHECK (strand_arrangement IN ('parallel_in_register','antiparallel','out_of_register','mixed','unknown')),
    CHECK (twist_handedness IN ('left','right','unknown')),
    CHECK (pathogenicity IN ('pathogenic','functional','non_pathogenic','unknown')),
    CHECK (structure_type IN ('fibril','oligomer','crystal','aggregate','monomer','unknown'))
);

-- Autoincrement (SQLite syntax; adapt for PostgreSQL)
-- PostgreSQL: use GENERATED ALWAYS AS IDENTITY or SERIAL instead of INTEGER NOT NULL with AUTOINCREMENT

CREATE INDEX IF NOT EXISTS idx_si_uniprot ON sequence_instances (uniprot_id);
CREATE INDEX IF NOT EXISTS idx_si_is_amyloid ON sequence_instances (is_amyloid);
CREATE INDEX IF NOT EXISTS idx_si_evidence_type ON sequence_instances (evidence_type);
CREATE INDEX IF NOT EXISTS idx_si_confidence ON sequence_instances (confidence);
CREATE INDEX IF NOT EXISTS idx_si_secondary_struct ON sequence_instances (secondary_structure_class);
CREATE INDEX IF NOT EXISTS idx_si_pathogenicity ON sequence_instances (pathogenicity);
CREATE INDEX IF NOT EXISTS idx_si_aggregate_type ON sequence_instances (aggregate_type);
CREATE INDEX IF NOT EXISTS idx_si_fold_topology ON sequence_instances (fold_topology);
CREATE INDEX IF NOT EXISTS idx_si_structure_type ON sequence_instances (structure_type);
CREATE INDEX IF NOT EXISTS idx_si_pdb ON sequence_instances (pdb_id) WHERE pdb_id IS NOT NULL AND pdb_id != '';


-- =============================================================================
-- LAYER 3: EVIDENCE / PROVENANCE TABLES
-- =============================================================================

CREATE TABLE IF NOT EXISTS publications (
    pub_id          INTEGER     NOT NULL,
    doi             TEXT,                        -- without https://doi.org/ prefix
    pmid            TEXT,
    reference_text  TEXT,
    --
    PRIMARY KEY (pub_id),
    UNIQUE (doi, pmid)
);

CREATE INDEX IF NOT EXISTS idx_pub_doi ON publications (doi) WHERE doi IS NOT NULL AND doi != '';
CREATE INDEX IF NOT EXISTS idx_pub_pmid ON publications (pmid) WHERE pmid IS NOT NULL AND pmid != '';


CREATE TABLE IF NOT EXISTS structures (
    struct_id       INTEGER     NOT NULL,
    pdb_id          TEXT,                        -- 4-character PDB accession
    emdb_id         TEXT,
    resolution      TEXT,                        -- free text; not parsed to float
    method_universal TEXT,
    --
    PRIMARY KEY (struct_id),
    UNIQUE (pdb_id, emdb_id)
);

CREATE INDEX IF NOT EXISTS idx_struct_pdb ON structures (pdb_id) WHERE pdb_id IS NOT NULL AND pdb_id != '';


CREATE TABLE IF NOT EXISTS raw_observations (
    obs_id          INTEGER     NOT NULL,
    -- Link to consensus (NULL if this observation was a conflict loser)
    seq_id          INTEGER,
    -- Source identity
    record_id       TEXT,                        -- original identifier from source DB
    source_db       TEXT        NOT NULL,
    -- Labels as reported by THIS source (may differ from consensus)
    experimental_label TEXT,
    category        TEXT,
    -- Method evidence
    raw_method      TEXT,
    method_universal TEXT,
    evidence_type   TEXT,
    evidence_weight REAL,
    confidence      INTEGER,
    -- Structural refs
    pdb_id          TEXT,
    emdb_id         TEXT,
    resolution      TEXT,
    -- Publication refs
    doi             TEXT,
    pmid            TEXT,
    reference       TEXT,
    -- Context
    disease         TEXT,
    tissue          TEXT,
    mutation        TEXT,
    notes           TEXT,
    --
    PRIMARY KEY (obs_id),
    FOREIGN KEY (seq_id) REFERENCES sequence_instances (seq_id),
    FOREIGN KEY (source_db) REFERENCES source_databases (source_db),
    CHECK (confidence IS NULL OR confidence BETWEEN 0 AND 100),
    CHECK (evidence_weight IS NULL OR evidence_weight >= 0.0)
);

CREATE INDEX IF NOT EXISTS idx_ro_seq_id ON raw_observations (seq_id);
CREATE INDEX IF NOT EXISTS idx_ro_source_db ON raw_observations (source_db);
CREATE INDEX IF NOT EXISTS idx_ro_record_id ON raw_observations (record_id);


-- =============================================================================
-- LAYER 4: COMPUTED FEATURE TABLES (1:1 with sequence_instances)
-- =============================================================================

CREATE TABLE IF NOT EXISTS physicochemical_features (
    seq_id                      INTEGER     NOT NULL,
    -- Basic
    length                      INTEGER     NOT NULL,
    molecular_weight            REAL        NOT NULL,   -- Da; sum(residue MW) - (n-1)*18.015
    -- Hydrophobicity (Kyte-Doolittle scale)
    hydrophobicity_mean         REAL        NOT NULL,
    hydrophobicity_std          REAL        NOT NULL,
    hydrophobicity_max          REAL        NOT NULL,
    hydrophobicity_min          REAL        NOT NULL,
    -- Charge (pH 7: R,K=+1; D,E=-1; H=+0.1)
    net_charge                  REAL        NOT NULL,
    positive_residues           INTEGER     NOT NULL,   -- count of R, K, H
    negative_residues           INTEGER     NOT NULL,   -- count of D, E
    charge_density              REAL        NOT NULL,   -- net_charge / length
    -- Structure propensity
    beta_propensity_mean        REAL        NOT NULL,   -- Chou-Fasman
    beta_propensity_max         REAL        NOT NULL,
    aggregation_propensity_mean REAL        NOT NULL,   -- custom scale (see features/physicochemical.py)
    aggregation_propensity_max  REAL        NOT NULL,
    -- Composition fractions
    aromatic_fraction           REAL        NOT NULL,   -- F, W, Y
    aliphatic_fraction          REAL        NOT NULL,   -- A, V, I, L
    polar_fraction              REAL        NOT NULL,   -- S, T, N, Q
    charged_fraction            REAL        NOT NULL,   -- D, E, K, R, H
    -- Sequence motifs
    has_polyq                   BOOLEAN     NOT NULL,   -- QQQ or longer
    has_glycine_rich            BOOLEAN     NOT NULL,   -- G fraction > 20%
    has_proline                 BOOLEAN     NOT NULL,   -- any P
    --
    PRIMARY KEY (seq_id),
    FOREIGN KEY (seq_id) REFERENCES sequence_instances (seq_id) ON DELETE CASCADE
);


CREATE TABLE IF NOT EXISTS sequence_composition (
    seq_id          INTEGER     NOT NULL,
    -- Per-residue frequencies (0.0-1.0), 20 standard AAs
    aa_A            REAL        NOT NULL DEFAULT 0.0,
    aa_C            REAL        NOT NULL DEFAULT 0.0,
    aa_D            REAL        NOT NULL DEFAULT 0.0,
    aa_E            REAL        NOT NULL DEFAULT 0.0,
    aa_F            REAL        NOT NULL DEFAULT 0.0,
    aa_G            REAL        NOT NULL DEFAULT 0.0,
    aa_H            REAL        NOT NULL DEFAULT 0.0,
    aa_I            REAL        NOT NULL DEFAULT 0.0,
    aa_K            REAL        NOT NULL DEFAULT 0.0,
    aa_L            REAL        NOT NULL DEFAULT 0.0,
    aa_M            REAL        NOT NULL DEFAULT 0.0,
    aa_N            REAL        NOT NULL DEFAULT 0.0,
    aa_P            REAL        NOT NULL DEFAULT 0.0,
    aa_Q            REAL        NOT NULL DEFAULT 0.0,
    aa_R            REAL        NOT NULL DEFAULT 0.0,
    aa_S            REAL        NOT NULL DEFAULT 0.0,
    aa_T            REAL        NOT NULL DEFAULT 0.0,
    aa_V            REAL        NOT NULL DEFAULT 0.0,
    aa_W            REAL        NOT NULL DEFAULT 0.0,
    aa_Y            REAL        NOT NULL DEFAULT 0.0,
    -- Grouped fractions
    tiny_fraction   REAL        NOT NULL DEFAULT 0.0,   -- A, G, S
    small_fraction  REAL        NOT NULL DEFAULT 0.0,   -- A, C, D, G, N, P, S, T, V
    large_fraction  REAL        NOT NULL DEFAULT 0.0,   -- F, H, K, R, W, Y
    -- Note: dipeptide frequencies are computed in pipeline but NOT exported to DB
    --       (only 3 sparse dipeptides appear in current TSV output)
    --
    PRIMARY KEY (seq_id),
    FOREIGN KEY (seq_id) REFERENCES sequence_instances (seq_id) ON DELETE CASCADE
);


-- =============================================================================
-- LAYER 5: AUXILIARY — LABEL CONFLICTS
-- =============================================================================

CREATE TABLE IF NOT EXISTS label_conflicts (
    conflict_id     INTEGER     NOT NULL,
    -- Biological identity (raw key from unifier.detect_conflicts())
    dedup_key       TEXT        NOT NULL,        -- string repr of the Python tuple key
    sequence        TEXT,                        -- extracted from dedup_key
    uniprot_id      TEXT,                        -- extracted from dedup_key
    region_start    INTEGER,                     -- extracted from dedup_key
    region_end      INTEGER,                     -- extracted from dedup_key
    -- Conflict data
    labels          TEXT        NOT NULL,        -- JSON array: e.g. ["amyloid","non-amyloid"]
    sources         TEXT        NOT NULL,        -- JSON array: e.g. ["WALTZ-DB","WALTZ-DB"]
    -- Resolution status
    resolved        BOOLEAN     NOT NULL DEFAULT 0,
    resolution_notes TEXT,
    --
    PRIMARY KEY (conflict_id),
    UNIQUE (dedup_key)
    -- NOT linked to sequence_instances by design:
    -- conflicted entries are excluded from the main table until manually resolved
);


-- =============================================================================
-- VIEWS (equivalent to current SQLite views, adapted to new table name)
-- =============================================================================

CREATE VIEW IF NOT EXISTS amyloid_entries AS
    SELECT * FROM sequence_instances WHERE is_amyloid = 1;

CREATE VIEW IF NOT EXISTS experimental_only AS
    SELECT * FROM sequence_instances
    WHERE evidence_type NOT IN ('computational', 'literature_curated')
       OR evidence_type IS NULL;

CREATE VIEW IF NOT EXISTS high_confidence AS
    SELECT * FROM sequence_instances WHERE confidence >= 70;

CREATE VIEW IF NOT EXISTS cross_beta_fibrils AS
    SELECT * FROM sequence_instances WHERE secondary_structure_class = 'cross_beta';

CREATE VIEW IF NOT EXISTS cross_alpha_fibrils AS
    SELECT * FROM sequence_instances WHERE secondary_structure_class = 'cross_alpha';

CREATE VIEW IF NOT EXISTS ex_vivo_structures AS
    SELECT * FROM sequence_instances WHERE aggregate_type IN ('ex_vivo', 'patient_derived');

CREATE VIEW IF NOT EXISTS steric_zippers AS
    SELECT * FROM sequence_instances WHERE interface_type = 'steric_zipper';

CREATE VIEW IF NOT EXISTS with_fold_topology AS
    SELECT * FROM sequence_instances WHERE fold_topology != 'unknown';

CREATE VIEW IF NOT EXISTS ml_features AS
    SELECT
        si.*,
        pf.length, pf.molecular_weight,
        pf.hydrophobicity_mean, pf.hydrophobicity_std, pf.hydrophobicity_max, pf.hydrophobicity_min,
        pf.net_charge, pf.positive_residues, pf.negative_residues, pf.charge_density,
        pf.beta_propensity_mean, pf.beta_propensity_max,
        pf.aggregation_propensity_mean, pf.aggregation_propensity_max,
        pf.aromatic_fraction, pf.aliphatic_fraction, pf.polar_fraction, pf.charged_fraction,
        pf.has_polyq, pf.has_glycine_rich, pf.has_proline,
        sc.aa_A, sc.aa_C, sc.aa_D, sc.aa_E, sc.aa_F, sc.aa_G, sc.aa_H, sc.aa_I,
        sc.aa_K, sc.aa_L, sc.aa_M, sc.aa_N, sc.aa_P, sc.aa_Q, sc.aa_R, sc.aa_S,
        sc.aa_T, sc.aa_V, sc.aa_W, sc.aa_Y,
        sc.tiny_fraction, sc.small_fraction, sc.large_fraction
    FROM sequence_instances si
    LEFT JOIN physicochemical_features pf ON si.seq_id = pf.seq_id
    LEFT JOIN sequence_composition sc ON si.seq_id = sc.seq_id
    WHERE si.sequence IS NOT NULL AND length(si.sequence) >= 4;
