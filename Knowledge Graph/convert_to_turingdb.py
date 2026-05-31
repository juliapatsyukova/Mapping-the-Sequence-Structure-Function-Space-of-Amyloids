"""
convert_to_turingdb.py
======================
Converts Amyloid DB graph files into TuringDB-compatible JSONL format.

Reads:
    graph_nodes.csv   — node registry (node_id, node_type, attributes_json)
    graph_edges.csv   — edge registry (source, target, relation_type, attributes_json)

Writes:
    amyloid_graph_turingdb.jsonl — TuringDB LOAD JSONL format

TuringDB JSONL format requirements
-----------------------------------
Nodes:
    {"type": "node", "id": "<int>", "labels": ["<NodeType>"], "properties": {...}}
    - id must be a non-negative integer (as string or int), starting at 0, incrementing by 1
    - labels must be an array with at least one element

Relationships:
    {"type": "relationship", "id": "<int>", "label": "<EDGE_TYPE>",
     "start": {"id": "<int>"}, "end": {"id": "<int>"}, "properties": {...}}
    - id must be a non-negative integer (as string or int), starting at 0, incrementing by 1
    - label must be a single string (not an array)
    - start/end reference node integer ids

Display labels
--------------
Each node properties dict includes three synonymous human-readable label fields so
the TuringDB UI and any downstream tool can surface readable text instead of numeric IDs:
    name          — primary display value
    label         — same value (alias)
    display_name  — same value (alias)

Label rules by node type:
    SequenceInstance  → protein_name [region] (enriched from PART_OF→Protein in two-pass load)
                        fallback chain: uniprot_id [region] → short seq (≤30 AA) → truncated seq… → Region N–M → Seq:{hash8}
    Observation       → "<source_db> | <experimental_label> | <method_universal>"
                        (method omitted if empty)
    Protein           → protein_name (fallback: uniprot_id)
    SourceDatabase    → source_db
    ExperimentMethod  → method_universal
    Publication       → reference string (fallback: DOI stripped of prefix, PMID, "publication")
    Structure         → pdb_id (resolution Å) (fallback: emdb_id, "structure")

Usage:
    python convert_to_turingdb.py [--nodes graph_nodes.csv] [--edges graph_edges.csv] [--out <dir>]

Load into TuringDB after running:
    cp amyloid_graph_turingdb.jsonl ~/.turing/data/
    turingdb
    > LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db
"""

import argparse
import csv
import html as html_mod
import json
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Display label computation
# ---------------------------------------------------------------------------

def _s(v) -> str:
    """Return stripped string, empty string for None/whitespace."""
    return str(v).strip() if v is not None else ""


def compute_display_label(node_type: str, attrs: dict, original_id: str) -> str:
    """Return a human-readable display label for a node."""

    if node_type == "SequenceInstance":
        seq          = _s(attrs.get("sequence"))
        uid          = _s(attrs.get("uniprot_id"))
        protein_name = _s(attrs.get("protein_name"))  # pre-enriched from PART_OF→Protein
        region_start = _s(attrs.get("region_start"))
        region_end   = _s(attrs.get("region_end"))

        def region_suffix() -> str:
            if region_start and region_end:
                return f" [{region_start}–{region_end}]"
            if region_start:
                return f" [from {region_start}]"
            if region_end:
                return f" [to {region_end}]"
            return ""

        # Short peptide with no protein context — sequence is the identifier
        if seq and len(seq) <= 30 and not uid and not protein_name:
            return seq
        if protein_name:
            return f"{protein_name}{region_suffix()}"
        if uid:
            return f"{uid}{region_suffix()}"
        if seq and len(seq) <= 30:
            return seq
        if seq:
            return f"{seq[:20]}…{region_suffix()}"
        if region_start or region_end:
            return f"Region {region_start or '?'}–{region_end or '?'}"
        return f"Seq:{original_id[-8:]}"

    if node_type == "Observation":
        src = _s(attrs.get("source_db"))
        lbl = _s(attrs.get("experimental_label"))
        mth = _s(attrs.get("method_universal"))
        parts = [p for p in [src, lbl, mth] if p]
        return " | ".join(parts) if parts else f"Observation:{original_id}"

    if node_type == "Protein":
        name = _s(attrs.get("protein_name"))
        if name:
            return name
        uid = _s(attrs.get("uniprot_id"))
        if uid:
            return uid
        return f"Protein:{original_id}"

    if node_type == "SourceDatabase":
        return _s(attrs.get("source_db")) or f"SourceDatabase:{original_id}"

    if node_type == "ExperimentMethod":
        return _s(attrs.get("method_universal")) or f"ExperimentMethod:{original_id}"

    if node_type == "Publication":
        ref = _s(attrs.get("reference"))
        if ref:
            return html_mod.unescape(ref)
        doi = _s(attrs.get("doi"))
        if doi:
            doi = doi.replace("https://doi.org/", "").replace("http://doi.org/", "")
            return doi
        pmid = _s(attrs.get("pmid"))
        if pmid:
            return f"PMID:{pmid}"
        return "publication"

    if node_type == "Structure":
        pdb        = _s(attrs.get("pdb_id"))
        emdb       = _s(attrs.get("emdb_id"))
        resolution = _s(attrs.get("resolution"))
        res_str    = f" ({resolution}Å)" if resolution else ""
        if pdb:
            return f"{pdb}{res_str}"
        if emdb:
            return f"{emdb}{res_str}"
        return "structure"

    # Unknown node type — use original string ID as fallback
    return original_id


def main():
    script_dir = Path(__file__).parent

    parser = argparse.ArgumentParser(
        description="Convert Amyloid DB graph to TuringDB JSONL format."
    )
    parser.add_argument(
        "--nodes",
        default=str(script_dir / "graph_nodes.csv"),
        help="Path to graph_nodes.csv",
    )
    parser.add_argument(
        "--edges",
        default=str(script_dir / "graph_edges.csv"),
        help="Path to graph_edges.csv",
    )
    parser.add_argument(
        "--out",
        default=str(script_dir),
        help="Output directory",
    )
    args = parser.parse_args()

    out_path = Path(args.out) / "amyloid_graph_turingdb.jsonl"

    # ------------------------------------------------------------------
    # Pass 1: Load all nodes into a map (labels computed after enrichment)
    # ------------------------------------------------------------------
    print(f"Pass 1: Loading nodes from {args.nodes} ...")
    node_info_map: dict[str, tuple[str, dict]] = {}  # original_id -> (node_type, props)
    node_order: list[str] = []  # preserve CSV order for stable integer IDs

    with open(args.nodes, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            original_id = row["node_id"]
            node_info_map[original_id] = (row["node_type"], json.loads(row["attributes_json"]))
            node_order.append(original_id)

    print(f"  {len(node_info_map):,} nodes loaded")

    # ------------------------------------------------------------------
    # Pass 2: Load edges; enrich SequenceInstance nodes via PART_OF→Protein
    # ------------------------------------------------------------------
    print(f"Pass 2: Loading edges from {args.edges} ...")
    raw_edges: list[dict] = []

    with open(args.edges, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            raw_edges.append({
                "source":   row["source"],
                "target":   row["target"],
                "relation": row["relation_type"],
                "props":    json.loads(row["attributes_json"]),
            })

    enriched = 0
    for edge in raw_edges:
        if edge["relation"] != "PART_OF":
            continue
        src_id, tgt_id = edge["source"], edge["target"]
        if src_id not in node_info_map or tgt_id not in node_info_map:
            continue
        src_type, src_props = node_info_map[src_id]
        tgt_type, tgt_props = node_info_map[tgt_id]
        if src_type == "SequenceInstance" and tgt_type == "Protein":
            protein_name = _s(tgt_props.get("protein_name"))
            if protein_name and not src_props.get("protein_name"):
                src_props["protein_name"] = protein_name
                enriched += 1

    print(f"  {len(raw_edges):,} edges loaded; "
          f"{enriched} SequenceInstance nodes enriched with protein_name")

    # ------------------------------------------------------------------
    # 3. Assign integer IDs and compute display labels
    # ------------------------------------------------------------------
    string_id_to_int: dict[str, int] = {}
    node_records = []

    for original_id in node_order:
        int_id = len(string_id_to_int)
        string_id_to_int[original_id] = int_id
        node_type, props = node_info_map[original_id]

        display = compute_display_label(node_type, props, original_id)
        props["name"]         = display
        props["label"]        = display
        props["display_name"] = display
        props["_original_id"] = original_id

        node_records.append({
            "type":       "node",
            "id":         str(int_id),
            "labels":     [node_type],
            "properties": props,
        })

    print(f"  {len(node_records):,} nodes indexed")

    # ------------------------------------------------------------------
    # 4. Convert edge records
    # ------------------------------------------------------------------
    edge_records = []
    skipped = 0

    for edge in raw_edges:
        src_str = edge["source"]
        tgt_str = edge["target"]

        if src_str not in string_id_to_int:
            print(f"  WARNING: unknown source node '{src_str}' — skipping edge",
                  file=sys.stderr)
            skipped += 1
            continue
        if tgt_str not in string_id_to_int:
            print(f"  WARNING: unknown target node '{tgt_str}' — skipping edge",
                  file=sys.stderr)
            skipped += 1
            continue

        edge_records.append({
            "type":       "relationship",
            "id":         str(len(edge_records)),
            "label":      edge["relation"],
            "start":      {"id": str(string_id_to_int[src_str])},
            "end":        {"id": str(string_id_to_int[tgt_str])},
            "properties": edge["props"],
        })

    print(f"  {len(edge_records):,} edges converted"
          + (f", {skipped} skipped (dangling references)" if skipped else ""))

    # ------------------------------------------------------------------
    # 5. Write output — nodes first, then edges
    # ------------------------------------------------------------------
    print(f"Writing {out_path} ...")
    with open(out_path, "w", encoding="utf-8") as f:
        for record in node_records:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")
        for record in edge_records:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")

    total = len(node_records) + len(edge_records)
    print(f"  Done: {total:,} lines ({len(node_records):,} nodes + {len(edge_records):,} edges)")
    print(f"\nTo load into TuringDB:")
    print(f"  cp {out_path} ~/.turing/data/")
    print(f"  turingdb")
    print(f"  > LOAD JSONL 'amyloid_graph_turingdb.jsonl' AS amyloid_db")


if __name__ == "__main__":
    main()
