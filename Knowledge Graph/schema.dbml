// Amyloid Database Schema - DBML for dbdiagram.io
// Generated: 2026-04-09  |  Source: schema.sql
// Views and pre-populated inserts are in schema.sql only.
// Multi-column UNIQUE and CHECK constraints are captured in column notes.

// =============================================================================
// LAYER 1: REFERENCE / LOOKUP TABLES
// =============================================================================

Table source_databases {
    source_db       varchar     [pk, not null, note: 'e.g. waltzdb, amyloid_atlas']
    entry_count     int         [note: 'Informational count; not enforced by DB']
    description     varchar

    Note: 'One row per upstream data source. PK used as FK in raw_observations.'
}

Table evidence_methods {
    method_id       int         [pk, not null]
    method_universal varchar    [not null, unique, note: 'e.g. Cryo-EM, ThT binding']
    evidence_type   varchar     [not null, note: 'CHECK: structural | kinetic | staining_binding | literature_curated | computational']
    evidence_weight float       [not null, note: 'CHECK >= 0.0. Values: 3.0 | 2.0 | 1.0 | 0.5 | 0.0']
    tier            int         [not null, note: 'CHECK: 0-5 (config.TIER_CONFIDENCE)']
    confidence      int         [not null, note: 'CHECK: 0-100. Values: 10|35|50|70|80|90']
    pattern         varchar     [note: 'Regex pattern from config.METHODS_RULES (21 rules total)']

    Note: 'Pre-populated from config.py METHODS_RULES. 21 rows inserted at schema creation.'
}


// =============================================================================
// LAYER 2: BIOLOGICAL IDENTITY TABLES
// =============================================================================

Table proteins {
    uniprot_id      varchar     [pk, not null, note: 'e.g. P05067, Q13148']
    protein_name    varchar     [note: 'As reported by source; not normalized']
    organism        varchar     [note: 'e.g. Homo sapiens']
    protein_family  varchar     [note: 'Assigned by keyword match from config.PROTEIN_FAMILIES']

    indexes {
        protein_family [name: 'idx_proteins_family']
        organism [name: 'idx_proteins_organism']
    }

    Note: 'One row per UniProt accession. Referenced by sequence_instances.uniprot_id.'
}

Table sequence_instances {
    seq_id                      int         [pk, not null, note: 'AUTOINCREMENT (SQLite) / GENERATED ALWAYS AS IDENTITY (PostgreSQL)']
    sequence                    varchar     [not null, note: 'Cleaned AA sequence']
    uniprot_id                  varchar     [ref: > proteins.uniprot_id, note: 'NULL for ~26% of entries; no UniProt ID available from source']
    region_start                int         [note: '1-based; NULL when not available']
    region_end                  int         [note: '1-based; NULL when not available']
    sequence_length             int         [not null, note: 'len(sequence). UNIQUE(sequence, uniprot_id, region_start, region_end)']
    is_amyloid                  boolean     [not null, note: 'CHECK: 0 or 1. Consensus winner result.']
    experimental_label          varchar     [note: 'Winning label string; not normalized to enum']
    category                    varchar     [note: 'Source category e.g. hexapeptide']
    evidence_weight             float       [not null, default: 0.0, note: 'CHECK >= 0.0']
    confidence                  int         [not null, default: 0, note: 'CHECK: 0-100']
    evidence_type               varchar     [note: 'Inherited from winning observation method']
    method_universal            varchar     [note: 'Canonical method name']
    raw_method                  varchar     [note: 'Raw method text from source']
    disease                     varchar
    tissue                      varchar
    mutation                    varchar
    notes                       varchar
    pdb_id                      varchar     [note: '4-character PDB accession; from enrichment']
    emdb_id                     varchar
    resolution                  varchar     [note: 'Free text; not parsed to float. See SCHEMA_NOTES 1.3']
    secondary_structure_class   varchar     [not null, default: 'unknown', note: 'Level 1. CHECK: cross_beta | cross_alpha | unknown. 72% of entries = unknown.']
    interface_type              varchar     [not null, default: 'unknown', note: 'Level 2. steric_zipper | LARKS | cross_alpha_sheet | unknown']
    strand_arrangement          varchar     [not null, default: 'unknown', note: 'Level 3. CHECK: parallel_in_register | antiparallel | out_of_register | mixed | unknown']
    fold_topology               varchar     [not null, default: 'unknown', note: 'Level 4. beta_arcade | greek_key | beta_solenoid | beta_helix | beta_sandwich | superpleated | other | unknown']
    fold_shape                  varchar     [not null, default: 'unknown', note: 'Level 5. J_shaped | S_shaped | C_shaped | LS_shaped | kidney | elongated | compact | other | unknown']
    zipper_class                varchar     [note: 'Level 6. 1-8 or 1-15; empty = unknown']
    protofilament_symmetry      varchar     [not null, default: 'unknown', note: 'Level 7. C1 | C2 | C3 | screw_21 | unknown']
    protofilament_count         int         [not null, default: 0, note: 'Level 8. 0 = unknown']
    twist_handedness            varchar     [not null, default: 'unknown', note: 'Level 9. CHECK: left | right | unknown']
    aggregate_type              varchar     [not null, default: 'unknown', note: 'Level 10. ex_vivo | in_vitro | recombinant | synthetic | unknown']
    polymorph_name              varchar     [note: 'Level 11. e.g. PHF, rod, twister; empty = unknown']
    pathogenicity               varchar     [not null, default: 'unknown', note: 'CHECK: pathogenic | functional | non_pathogenic | unknown']
    structure_type              varchar     [not null, default: 'unknown', note: 'Legacy. CHECK: fibril | oligomer | crystal | aggregate | monomer | unknown. See SCHEMA_NOTES 2.4']

    indexes {
        uniprot_id [name: 'idx_si_uniprot']
        is_amyloid [name: 'idx_si_is_amyloid']
        evidence_type [name: 'idx_si_evidence_type']
        confidence [name: 'idx_si_confidence']
        secondary_structure_class [name: 'idx_si_secondary_struct']
        pathogenicity [name: 'idx_si_pathogenicity']
        aggregate_type [name: 'idx_si_aggregate_type']
        fold_topology [name: 'idx_si_fold_topology']
        structure_type [name: 'idx_si_structure_type']
        pdb_id [name: 'idx_si_pdb']
    }

    Note: 'One row per unique (sequence, uniprot_id, region_start, region_end). Consensus winner only; losers stored in raw_observations.'
}


// =============================================================================
// LAYER 3: EVIDENCE / PROVENANCE TABLES
// =============================================================================

Table publications {
    pub_id          int         [pk, not null]
    doi             varchar     [note: 'Without https://doi.org/ prefix']
    pmid            varchar
    reference_text  varchar

    indexes {
        doi [name: 'idx_pub_doi']
        pmid [name: 'idx_pub_pmid']
    }

    Note: 'UNIQUE(doi, pmid). One row per publication. Referenced by raw_observations via denormalized doi/pmid columns.'
}

Table structures {
    struct_id       int         [pk, not null]
    pdb_id          varchar     [note: '4-character PDB accession']
    emdb_id         varchar
    resolution      varchar     [note: 'Free text; not parsed. e.g. 2.1, 3.4 Angstrom, NA']
    method_universal varchar

    indexes {
        pdb_id [name: 'idx_struct_pdb']
    }

    Note: 'UNIQUE(pdb_id, emdb_id). Structural depositions referenced by sequence_instances.'
}

Table raw_observations {
    obs_id          int         [pk, not null]
    seq_id          int         [ref: > sequence_instances.seq_id, note: 'NULL if observation was a conflict loser dropped during consensus']
    record_id       varchar     [note: 'Original identifier from upstream source DB']
    source_db       varchar     [not null, ref: > source_databases.source_db]
    experimental_label varchar  [note: 'Label as reported by this source; may differ from consensus']
    category        varchar
    raw_method      varchar
    method_universal varchar
    evidence_type   varchar
    evidence_weight float       [note: 'CHECK: NULL or >= 0.0']
    confidence      int         [note: 'CHECK: NULL or 0-100']
    pdb_id          varchar
    emdb_id         varchar
    resolution      varchar
    doi             varchar
    pmid            varchar
    reference       varchar
    disease         varchar
    tissue          varchar
    mutation        varchar
    notes           varchar

    indexes {
        seq_id [name: 'idx_ro_seq_id']
        source_db [name: 'idx_ro_source_db']
        record_id [name: 'idx_ro_record_id']
    }

    Note: 'All pre-consensus observations, including consensus losers (seq_id = NULL). Preserves full provenance. See SCHEMA_NOTES 2.2.'
}


// =============================================================================
// LAYER 4: COMPUTED FEATURE TABLES (1:1 with sequence_instances)
// =============================================================================

Table physicochemical_features {
    seq_id                      int         [pk, not null, ref: - sequence_instances.seq_id, note: 'CASCADE DELETE']
    length                      int         [not null]
    molecular_weight            float       [not null, note: 'Da; sum(residue MW) minus (n-1)*18.015']
    hydrophobicity_mean         float       [not null]
    hydrophobicity_std          float       [not null]
    hydrophobicity_max          float       [not null]
    hydrophobicity_min          float       [not null]
    net_charge                  float       [not null]
    positive_residues           int         [not null, note: 'Count of R, K, H']
    negative_residues           int         [not null, note: 'Count of D, E']
    charge_density              float       [not null, note: 'net_charge / length']
    beta_propensity_mean        float       [not null, note: 'Chou-Fasman scale']
    beta_propensity_max         float       [not null]
    aggregation_propensity_mean float       [not null, note: 'Custom scale; see features/physicochemical.py']
    aggregation_propensity_max  float       [not null]
    aromatic_fraction           float       [not null, note: 'F, W, Y']
    aliphatic_fraction          float       [not null, note: 'A, V, I, L']
    polar_fraction              float       [not null, note: 'S, T, N, Q']
    charged_fraction            float       [not null, note: 'D, E, K, R, H']
    has_polyq                   boolean     [not null, note: 'True if QQQ or longer run present']
    has_glycine_rich            boolean     [not null, note: 'True if G fraction > 20%']
    has_proline                 boolean     [not null, note: 'True if any P present']

    Note: '1:1 with sequence_instances. ON DELETE CASCADE. Computed by features/physicochemical.py.'
}

Table sequence_composition {
    seq_id          int         [pk, not null, ref: - sequence_instances.seq_id, note: 'CASCADE DELETE']
    aa_A            float       [not null, default: 0.0]
    aa_C            float       [not null, default: 0.0]
    aa_D            float       [not null, default: 0.0]
    aa_E            float       [not null, default: 0.0]
    aa_F            float       [not null, default: 0.0]
    aa_G            float       [not null, default: 0.0]
    aa_H            float       [not null, default: 0.0]
    aa_I            float       [not null, default: 0.0]
    aa_K            float       [not null, default: 0.0]
    aa_L            float       [not null, default: 0.0]
    aa_M            float       [not null, default: 0.0]
    aa_N            float       [not null, default: 0.0]
    aa_P            float       [not null, default: 0.0]
    aa_Q            float       [not null, default: 0.0]
    aa_R            float       [not null, default: 0.0]
    aa_S            float       [not null, default: 0.0]
    aa_T            float       [not null, default: 0.0]
    aa_V            float       [not null, default: 0.0]
    aa_W            float       [not null, default: 0.0]
    aa_Y            float       [not null, default: 0.0]
    tiny_fraction   float       [not null, default: 0.0, note: 'A, G, S']
    small_fraction  float       [not null, default: 0.0, note: 'A, C, D, G, N, P, S, T, V']
    large_fraction  float       [not null, default: 0.0, note: 'F, H, K, R, W, Y']

    Note: '1:1 with sequence_instances. ON DELETE CASCADE. 20 AA frequencies + 3 grouped fractions. Dipeptide frequencies are computed in pipeline but not exported.'
}


// =============================================================================
// LAYER 5: AUXILIARY — LABEL CONFLICTS
// =============================================================================

Table label_conflicts {
    conflict_id     int         [pk, not null]
    dedup_key       varchar     [not null, unique, note: 'String repr of the Python tuple key e.g. (KNFNYN, P05453, 1, 6)']
    sequence        varchar     [note: 'Extracted from dedup_key']
    uniprot_id      varchar     [note: 'Extracted from dedup_key; NOT a FK to proteins by design']
    region_start    int         [note: 'Extracted from dedup_key']
    region_end      int         [note: 'Extracted from dedup_key']
    labels          varchar     [not null, note: 'JSON array of labels e.g. amyloid, non-amyloid']
    sources         varchar     [not null, note: 'JSON array of source DBs parallel to labels']
    resolved        boolean     [not null, default: false]
    resolution_notes varchar

    Note: 'Intentionally NOT FK-linked to sequence_instances. Conflicted entries are excluded from the main table until manually resolved. 4 entries in current dataset. See SCHEMA_NOTES 2.5.'
}


// =============================================================================
// TABLE GROUPS (logical layers)
// =============================================================================

TableGroup layer_1_reference {
    source_databases
    evidence_methods
}

TableGroup layer_2_biological_identity {
    proteins
    sequence_instances
}

TableGroup layer_3_evidence_provenance {
    publications
    structures
    raw_observations
}

TableGroup layer_4_computed_features {
    physicochemical_features
    sequence_composition
}

TableGroup layer_5_auxiliary {
    label_conflicts
}
