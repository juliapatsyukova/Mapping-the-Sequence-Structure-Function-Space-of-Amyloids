"""
graph_build.py
==============
Builds an Evidence / Provenance Graph from the Amyloid DB outputs.

The graph has TWO layers that coexist in the same output files:

  LAYER 1 — Consensus layer (original)
  -------------------------------------
  Nodes: SequenceInstance, Protein, SourceDatabase, ExperimentMethod,
         Publication, Structure
  Edges: PART_OF, OBSERVED_IN, SUPPORTED_BY, REPORTED_IN, HAS_STRUCTURE

  SequenceInstance nodes are SOURCE-AGNOSTIC. They carry the consensus
  classification (is_amyloid, evidence_weight, confidence) and all
  structural classification fields (11 levels).

  LAYER 2 — Pre-consensus Observation layer (extension)
  -------------------------------------------------------
  Nodes: Observation  (new)
  Edges: OBSERVED_AS  (Observation → SequenceInstance)
         FROM_SOURCE  (Observation → SourceDatabase)
         SUPPORTED_BY (Observation → ExperimentMethod)
         REPORTED_IN  (Observation → Publication)

  One Observation node per TSV row = one consensus-winner record per
  dedup-key group. When multiple TSV rows share the same biological
  identity (sequence + uniprot_id + region), they produce separate
  Observation nodes pointing at the same SequenceInstance node —
  representing genuinely multi-source evidence.

  Conflict stubs (from conflicts.json) produce minimal Observation nodes
  for the 3 biologically meaningful unresolved conflicts.

  DATA LIMITATION: This is a partial reconstruction. The pipeline only
  exports consensus winners. Within each dedup group, lower-evidence
  observations were discarded before export and are not recoverable
  from the available output files.

Outputs (written to --out dir):
    graph_nodes.csv   — node_id, node_type, attributes_json
    graph_edges.csv   — source, target, relation_type, attributes_json
    graph.jsonl       — TuringDB-ready: one JSON object per line

Usage:
    python graph_build.py [--tsv PATH] [--conflicts PATH] [--out DIR]

Defaults:
    --tsv       ../database/consensus_unified.tsv
    --conflicts ../database/conflicts.json
    --out       .  (same directory as this script)

Requirements:
    pandas >= 1.3
    networkx (optional; only used if --write-graphml is passed)
"""

import argparse
import ast
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha8(s: str) -> str:
    """Return first 16 hex chars of SHA-256 of s."""
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def safe(v) -> str:
    """Coerce a value to string; return empty string for NaN / None / 'nan'."""
    if v is None:
        return ""
    if isinstance(v, float) and math.isnan(v):
        return ""
    s = str(v).strip()
    return "" if s.lower() == "nan" else s


# ---------------------------------------------------------------------------
# Node ID constructors — canonical, deterministic
# ---------------------------------------------------------------------------

def seq_node_id(sequence: str, uniprot_id: str, region_start, region_end) -> str:
    """SequenceInstance node ID: hashed biological key."""
    key = f"{safe(sequence)}|{safe(uniprot_id)}|{safe(region_start)}|{safe(region_end)}"
    return f"seq:{sha8(key)}"


def obs_node_id(source_db: str, record_id: str, sequence: str,
                region_start, region_end) -> str:
    """Observation node ID: hashed per-row identity (one per TSV row)."""
    key = f"{safe(source_db)}|{safe(record_id)}|{safe(sequence)}|{safe(region_start)}|{safe(region_end)}"
    return f"obs:{sha8(key)}"


def protein_node_id(uniprot_id: str, protein_name: str) -> str:
    uid = safe(uniprot_id)
    if uid:
        return f"protein:{uid}"
    name = safe(protein_name)
    if name:
        return f"protein:name:{sha8(name)}"
    return ""


def method_node_id(method_universal: str) -> str:
    m = safe(method_universal)
    return f"method:{sha8(m)}" if m else ""


def pub_node_id(doi: str, pmid: str) -> str:
    d, p = safe(doi), safe(pmid)
    if d:
        return f"pub:doi:{sha8(d)}"
    if p:
        return f"pub:pmid:{sha8(p)}"
    return ""


def struct_node_id(pdb_id: str, emdb_id: str) -> str:
    pid, eid = safe(pdb_id), safe(emdb_id)
    if pid:
        return f"struct:pdb:{pid.upper()}"
    if eid:
        return f"struct:emdb:{eid}"
    return ""


# ---------------------------------------------------------------------------
# Graph accumulator
# ---------------------------------------------------------------------------

class GraphBuilder:
    def __init__(self):
        # node_id -> (node_type, attrs_dict)
        self.nodes: dict[str, tuple[str, dict]] = {}
        # list of (source_id, target_id, relation_type, attrs_dict)
        self.edges: list[tuple[str, str, str, dict]] = []

    def add_node(self, node_id: str, node_type: str, attrs: dict):
        if not node_id:
            return
        if node_id not in self.nodes:
            self.nodes[node_id] = (node_type, attrs)
        else:
            # Merge: fill in any missing values without overwriting existing
            existing = self.nodes[node_id][1]
            for k, v in attrs.items():
                if k not in existing or not existing[k]:
                    existing[k] = v

    def add_edge(self, source: str, target: str, relation: str, attrs: dict):
        if not source or not target:
            return
        self.edges.append((source, target, relation, attrs))

    def deduplicate_edges(self):
        seen: set[str] = set()
        deduped = []
        for src, tgt, rel, attrs in self.edges:
            key = (src, tgt, rel, json.dumps(attrs, sort_keys=True))
            if key not in seen:
                seen.add(key)
                deduped.append((src, tgt, rel, attrs))
        self.edges = deduped

    def node_type_counts(self) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for _, (nt, _) in self.nodes.items():
            counts[nt] += 1
        return dict(counts)

    def edge_type_counts(self) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for _, _, rel, _ in self.edges:
            counts[rel] += 1
        return dict(counts)


# ---------------------------------------------------------------------------
# Layer 1: Consensus graph (original, unchanged logic)
# ---------------------------------------------------------------------------

def build_consensus_layer(df: pd.DataFrame, g: GraphBuilder) -> None:
    """
    Populate g with the consensus layer.

    One SequenceInstance per unique biological identity
    (sequence, uniprot_id, region_start, region_end).
    All structural classification and consensus evidence fields live here.
    SequenceInstance nodes are source-agnostic.
    """
    for _, row in df.iterrows():
        sequence     = safe(row.get("sequence", ""))
        uniprot_id   = safe(row.get("uniprot_id", ""))
        region_start = safe(row.get("region_start", ""))
        region_end   = safe(row.get("region_end", ""))
        protein_name = safe(row.get("protein_name", ""))
        organism     = safe(row.get("organism", ""))
        protein_fam  = safe(row.get("protein_family", ""))
        source_db    = safe(row.get("source_db", ""))
        is_amyloid   = safe(row.get("is_amyloid", ""))
        exp_label    = safe(row.get("experimental_label", ""))
        evidence_wt  = safe(row.get("evidence_weight", ""))
        confidence   = safe(row.get("confidence", ""))
        evidence_tp  = safe(row.get("evidence_type", ""))
        method_univ  = safe(row.get("method_universal", ""))
        raw_method   = safe(row.get("raw_method", ""))
        pdb_id       = safe(row.get("pdb_id", ""))
        emdb_id      = safe(row.get("emdb_id", ""))
        resolution   = safe(row.get("resolution", ""))
        doi          = safe(row.get("doi", ""))
        pmid         = safe(row.get("pmid", ""))
        reference    = safe(row.get("reference", ""))
        disease      = safe(row.get("disease", ""))
        mutation     = safe(row.get("mutation", ""))
        sec_struct   = safe(row.get("secondary_structure_class", ""))
        iface_type   = safe(row.get("interface_type", ""))
        strand_arr   = safe(row.get("strand_arrangement", ""))
        fold_top     = safe(row.get("fold_topology", ""))
        fold_shp     = safe(row.get("fold_shape", ""))
        zipper_cls   = safe(row.get("zipper_class", ""))
        proto_sym    = safe(row.get("protofilament_symmetry", ""))
        proto_cnt    = safe(row.get("protofilament_count", ""))
        twist_hand   = safe(row.get("twist_handedness", ""))
        agg_type     = safe(row.get("aggregate_type", ""))
        polymorph    = safe(row.get("polymorph_name", ""))
        pathogenic   = safe(row.get("pathogenicity", ""))
        struct_type  = safe(row.get("structure_type", ""))

        # -- SequenceInstance (source-agnostic, carries consensus result) --
        sid = seq_node_id(sequence, uniprot_id, region_start, region_end)
        g.add_node(sid, "SequenceInstance", {
            "sequence": sequence,
            "uniprot_id": uniprot_id,
            "protein_name": protein_name,
            "organism": organism,
            "region_start": region_start,
            "region_end": region_end,
            "sequence_length": str(len(sequence)) if sequence else "",
            "is_amyloid": is_amyloid,
            "experimental_label": exp_label,
            "evidence_weight": evidence_wt,
            "confidence": confidence,
            "evidence_type": evidence_tp,
            "disease": disease,
            "mutation": mutation,
            "secondary_structure_class": sec_struct,
            "interface_type": iface_type,
            "strand_arrangement": strand_arr,
            "fold_topology": fold_top,
            "fold_shape": fold_shp,
            "zipper_class": zipper_cls,
            "protofilament_symmetry": proto_sym,
            "protofilament_count": proto_cnt,
            "twist_handedness": twist_hand,
            "aggregate_type": agg_type,
            "polymorph_name": polymorph,
            "pathogenicity": pathogenic,
            "structure_type": struct_type,
            "graph_layer": "consensus",
        })

        # -- Protein --
        pid = protein_node_id(uniprot_id, protein_name)
        if pid:
            g.add_node(pid, "Protein", {
                "uniprot_id": uniprot_id,
                "protein_name": protein_name,
                "organism": organism,
                "protein_family": protein_fam,
            })
            g.add_edge(sid, pid, "PART_OF", {
                "region_start": region_start,
                "region_end": region_end,
                "graph_layer": "consensus",
            })

        # -- SourceDatabase --
        if source_db:
            db_id = f"db:{source_db}"
            g.add_node(db_id, "SourceDatabase", {"source_db": source_db})
            g.add_edge(sid, db_id, "OBSERVED_IN", {
                "experimental_label": exp_label,
                "source_db": source_db,
                "graph_layer": "consensus",
            })

        # -- ExperimentMethod --
        mid = method_node_id(method_univ)
        if mid:
            g.add_node(mid, "ExperimentMethod", {
                "method_universal": method_univ,
                "evidence_type": evidence_tp,
                "evidence_weight": evidence_wt,
            })
            g.add_edge(sid, mid, "SUPPORTED_BY", {
                "raw_method": raw_method,
                "evidence_weight": evidence_wt,
                "confidence": confidence,
                "graph_layer": "consensus",
            })

        # -- Publication --
        pubid = pub_node_id(doi, pmid)
        if pubid:
            g.add_node(pubid, "Publication", {
                "doi": doi,
                "pmid": pmid,
                "reference": reference,
            })
            g.add_edge(sid, pubid, "REPORTED_IN", {
                "doi": doi,
                "pmid": pmid,
                "graph_layer": "consensus",
            })

        # -- Structure --
        stid = struct_node_id(pdb_id, emdb_id)
        if stid:
            g.add_node(stid, "Structure", {
                "pdb_id": pdb_id,
                "emdb_id": emdb_id,
                "resolution": resolution,
                "method_universal": method_univ,
            })
            g.add_edge(sid, stid, "HAS_STRUCTURE", {
                "pdb_id": pdb_id,
                "emdb_id": emdb_id,
                "resolution": resolution,
                "graph_layer": "consensus",
            })


# ---------------------------------------------------------------------------
# Layer 2: Pre-consensus Observation layer (extension)
# ---------------------------------------------------------------------------

def build_observation_layer(df: pd.DataFrame, g: GraphBuilder) -> None:
    """
    Extend g with Observation nodes — one per TSV row.

    Each TSV row is a consensus-winner record from one dedup-key group.
    When multiple rows share the same biological identity, they produce
    separate Observation nodes all pointing at the same SequenceInstance.

    Observation node edges:
        Observation -OBSERVED_AS->  SequenceInstance
        Observation -FROM_SOURCE->  SourceDatabase
        Observation -SUPPORTED_BY-> ExperimentMethod
        Observation -REPORTED_IN->  Publication

    Note: HAS_STRUCTURE is not repeated from Observation because structures
    are already anchored to the SequenceInstance in the consensus layer and
    Observation nodes do not add new structural identity information.
    """
    for idx, row in df.iterrows():
        sequence     = safe(row.get("sequence", ""))
        uniprot_id   = safe(row.get("uniprot_id", ""))
        region_start = safe(row.get("region_start", ""))
        region_end   = safe(row.get("region_end", ""))
        source_db    = safe(row.get("source_db", ""))
        record_id    = safe(row.get("record_id", ""))
        exp_label    = safe(row.get("experimental_label", ""))
        evidence_wt  = safe(row.get("evidence_weight", ""))
        confidence   = safe(row.get("confidence", ""))
        evidence_tp  = safe(row.get("evidence_type", ""))
        method_univ  = safe(row.get("method_universal", ""))
        raw_method   = safe(row.get("raw_method", ""))
        pdb_id       = safe(row.get("pdb_id", ""))
        emdb_id      = safe(row.get("emdb_id", ""))
        resolution   = safe(row.get("resolution", ""))
        doi          = safe(row.get("doi", ""))
        pmid         = safe(row.get("pmid", ""))
        reference    = safe(row.get("reference", ""))
        category     = safe(row.get("category", ""))

        # -- Observation node --
        oid = obs_node_id(source_db, record_id, sequence, region_start, region_end)
        g.add_node(oid, "Observation", {
            "record_id": record_id,
            "source_db": source_db,
            "experimental_label": exp_label,
            "method_universal": method_univ,
            "raw_method": raw_method,
            "evidence_type": evidence_tp,
            "evidence_weight": evidence_wt,
            "confidence": confidence,
            "pdb_id": pdb_id,
            "emdb_id": emdb_id,
            "resolution": resolution,
            "doi": doi,
            "pmid": pmid,
            "category": category,
            "is_consensus_winner": "true",   # all TSV rows are consensus winners
            "data_completeness": "full",
            "graph_layer": "observation",
        })

        # Edge: Observation → SequenceInstance (OBSERVED_AS)
        sid = seq_node_id(sequence, uniprot_id, region_start, region_end)
        g.add_edge(oid, sid, "OBSERVED_AS", {
            "experimental_label": exp_label,
            "source_db": source_db,
            "graph_layer": "observation",
        })

        # Edge: Observation → SourceDatabase (FROM_SOURCE)
        if source_db:
            db_id = f"db:{source_db}"
            # SourceDatabase node was already registered in consensus layer;
            # add_node will merge without overwriting.
            g.add_node(db_id, "SourceDatabase", {"source_db": source_db})
            g.add_edge(oid, db_id, "FROM_SOURCE", {
                "source_db": source_db,
                "graph_layer": "observation",
            })

        # Edge: Observation → ExperimentMethod (SUPPORTED_BY)
        mid = method_node_id(method_univ)
        if mid:
            g.add_node(mid, "ExperimentMethod", {
                "method_universal": method_univ,
                "evidence_type": evidence_tp,
                "evidence_weight": evidence_wt,
            })
            g.add_edge(oid, mid, "SUPPORTED_BY", {
                "raw_method": raw_method,
                "evidence_weight": evidence_wt,
                "confidence": confidence,
                "graph_layer": "observation",
            })

        # Edge: Observation → Publication (REPORTED_IN)
        pubid = pub_node_id(doi, pmid)
        if pubid:
            g.add_node(pubid, "Publication", {
                "doi": doi,
                "pmid": pmid,
                "reference": reference,
            })
            g.add_edge(oid, pubid, "REPORTED_IN", {
                "doi": doi,
                "pmid": pmid,
                "graph_layer": "observation",
            })


def build_conflict_stubs(conflicts: list[dict], g: GraphBuilder) -> None:
    """
    Add minimal Observation stubs for unresolved label conflicts.

    Conflicts were dropped from consensus. Only source, label, and sequence
    identity are recoverable from conflicts.json. Full evidence attributes
    are not available (data_completeness: incomplete).

    Skips the CPAD-Structures artifact conflict (empty sequence key).
    Creates a SequenceInstance node marked graph_layer=conflict for each
    real conflict, plus one Observation node per source/label pair.
    """
    for conflict in conflicts:
        raw_key = conflict.get("key", "")
        labels  = conflict.get("labels", [])
        sources = conflict.get("sources", [])

        # Parse the Python-tuple-repr key: "('SEQ', 'UNIPROT', start, end)"
        # or the 2-element fallback key "('', 'source_db')"
        try:
            parsed = ast.literal_eval(raw_key)
        except (ValueError, SyntaxError):
            continue

        if not isinstance(parsed, tuple) or len(parsed) < 2:
            continue

        # Skip the CPAD-Structures artifact: 2-element key where first elem is empty
        if len(parsed) == 2 and not parsed[0]:
            continue

        # Real biological conflict: 4-element key (sequence, uniprot_id, start, end)
        if len(parsed) != 4:
            continue

        sequence, uniprot_id, region_start, region_end = parsed
        sequence     = str(sequence) if sequence else ""
        uniprot_id   = str(uniprot_id) if uniprot_id else ""
        region_start = str(region_start) if region_start is not None else ""
        region_end   = str(region_end) if region_end is not None else ""

        # SequenceInstance node for this conflict (NOT in consensus)
        sid = seq_node_id(sequence, uniprot_id, region_start, region_end)
        g.add_node(sid, "SequenceInstance", {
            "sequence": sequence,
            "uniprot_id": uniprot_id,
            "region_start": region_start,
            "region_end": region_end,
            "sequence_length": str(len(sequence)) if sequence else "",
            "is_amyloid": "unresolved",
            "experimental_label": "conflict",
            "evidence_weight": "",
            "confidence": "",
            "graph_layer": "conflict",
        })

        # One Observation stub per (source, label) pair
        for i, (src, lbl) in enumerate(zip(sources, labels)):
            oid = obs_node_id(
                src,
                f"conflict_{sha8(raw_key)}_{i}",
                sequence,
                region_start,
                region_end,
            )
            g.add_node(oid, "Observation", {
                "record_id": "",
                "source_db": src,
                "experimental_label": lbl,
                "method_universal": "",
                "raw_method": "",
                "evidence_type": "",
                "evidence_weight": "",
                "confidence": "",
                "pdb_id": "",
                "emdb_id": "",
                "doi": "",
                "pmid": "",
                "is_consensus_winner": "false",
                "data_completeness": "incomplete",
                "graph_layer": "conflict",
            })

            # OBSERVED_AS → conflict SequenceInstance
            g.add_edge(oid, sid, "OBSERVED_AS", {
                "experimental_label": lbl,
                "source_db": src,
                "graph_layer": "conflict",
            })

            # FROM_SOURCE → SourceDatabase
            if src:
                db_id = f"db:{src}"
                g.add_node(db_id, "SourceDatabase", {"source_db": src})
                g.add_edge(oid, db_id, "FROM_SOURCE", {
                    "source_db": src,
                    "graph_layer": "conflict",
                })


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def write_nodes_csv(g: GraphBuilder, path: str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "node_type", "attributes_json"])
        for node_id, (node_type, attrs) in sorted(g.nodes.items()):
            writer.writerow([node_id, node_type, json.dumps(attrs, ensure_ascii=False)])
    print(f"  Nodes: {len(g.nodes):,} → {path}")


def write_edges_csv(g: GraphBuilder, path: str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["source", "target", "relation_type", "attributes_json"])
        for src, tgt, rel, attrs in g.edges:
            writer.writerow([src, tgt, rel, json.dumps(attrs, ensure_ascii=False)])
    print(f"  Edges: {len(g.edges):,} → {path}")


def write_jsonl(g: GraphBuilder, path: str):
    with open(path, "w", encoding="utf-8") as f:
        for src, tgt, rel, attrs in g.edges:
            obj = {"source": src, "target": tgt, "type": rel, "properties": attrs}
            f.write(json.dumps(obj, ensure_ascii=False) + "\n")
    print(f"  JSONL:  {len(g.edges):,} lines → {path}")


def print_summary(g: GraphBuilder):
    print("  Node type counts:")
    for nt, cnt in sorted(g.node_type_counts().items()):
        print(f"    {nt}: {cnt:,}")
    print("  Edge type counts:")
    for rel, cnt in sorted(g.edge_type_counts().items()):
        print(f"    {rel}: {cnt:,}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    script_dir = Path(__file__).parent

    parser = argparse.ArgumentParser(
        description="Build two-layer Evidence/Provenance Graph from Amyloid DB outputs."
    )
    parser.add_argument(
        "--tsv",
        default=str(script_dir / "../database/consensus_unified.tsv"),
        help="Path to consensus_unified.tsv",
    )
    parser.add_argument(
        "--conflicts",
        default=str(script_dir / "../database/conflicts.json"),
        help="Path to conflicts.json",
    )
    parser.add_argument(
        "--out",
        default=str(script_dir),
        help="Output directory (default: same as this script)",
    )
    parser.add_argument(
        "--write-graphml",
        action="store_true",
        help="Also write a GraphML file (requires networkx)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Load TSV
    # -----------------------------------------------------------------------
    tsv_path = Path(args.tsv)
    if not tsv_path.exists():
        print(f"ERROR: TSV not found: {tsv_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {tsv_path} ...")
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    df.replace("nan", "", inplace=True)
    print(f"  Loaded {len(df):,} rows × {len(df.columns)} columns")

    # -----------------------------------------------------------------------
    # Load conflicts
    # -----------------------------------------------------------------------
    conflicts_path = Path(args.conflicts)
    conflicts: list[dict] = []
    if conflicts_path.exists():
        with open(conflicts_path, encoding="utf-8") as f:
            conflicts = json.load(f)
        print(f"  Loaded {len(conflicts)} conflict records from {conflicts_path}")
    else:
        print(f"  WARNING: conflicts.json not found at {conflicts_path} — skipping conflict stubs")

    # -----------------------------------------------------------------------
    # Build graph (both layers in one pass)
    # -----------------------------------------------------------------------
    g = GraphBuilder()

    print("\nLayer 1: Building consensus layer ...")
    build_consensus_layer(df, g)
    print(f"  After consensus layer:")
    print_summary(g)

    print("\nLayer 2: Building observation layer ...")
    build_observation_layer(df, g)

    print("\nLayer 2 (ext): Adding conflict stubs ...")
    build_conflict_stubs(conflicts, g)

    g.deduplicate_edges()

    print(f"\nFinal combined graph:")
    print_summary(g)

    # -----------------------------------------------------------------------
    # Export
    # -----------------------------------------------------------------------
    print("\nWriting outputs ...")
    write_nodes_csv(g, str(out_dir / "graph_nodes.csv"))
    write_edges_csv(g, str(out_dir / "graph_edges.csv"))
    write_jsonl(g,      str(out_dir / "graph.jsonl"))

    # Optional: GraphML
    if args.write_graphml:
        try:
            import networkx as nx
            G = nx.DiGraph()
            for node_id, (node_type, attrs) in g.nodes.items():
                G.add_node(node_id, node_type=node_type,
                           **{k: str(v) for k, v in attrs.items()})
            for src, tgt, rel, attrs in g.edges:
                G.add_edge(src, tgt, relation_type=rel,
                           **{k: str(v) for k, v in attrs.items()})
            gml_path = str(out_dir / "graph.graphml")
            nx.write_graphml(G, gml_path)
            print(f"  GraphML: {gml_path}")
        except ImportError:
            print("  WARNING: networkx not installed; skipping GraphML export.")

    print("\nDone.")


if __name__ == "__main__":
    main()
