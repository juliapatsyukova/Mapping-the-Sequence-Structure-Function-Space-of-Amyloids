"""
compute_metrics.py
==================
Computes sequence-level uncertainty and disagreement metrics for each
SequenceInstance node in the Amyloid DB graph.

Reads:
    graph_nodes.csv     — node registry (both layers)
    graph_edges.csv     — edge registry (both layers)

Writes:
    sequence_metrics.csv       — one row per SequenceInstance, all metrics
    graph_nodes.csv            — UPDATED in-place: SequenceInstance nodes
                                 gain metric attributes in attributes_json

Metric definitions
------------------
For SequenceInstance S, let O(S) = set of Observation nodes reachable
via OBSERVED_AS edges (both "observation" and "conflict" graph layers).

Per observation o:
    label(o)   = experimental_label
    source(o)  = source_db
    method(o)  = method_universal  (empty string for conflict stubs)
    weight(o)  = float(evidence_weight) if non-empty, else 0.0
    doi(o)     = doi field
    pmid(o)    = pmid field

L(S) = SequenceInstance.experimental_label  (consensus label)
W_win   = Σ weight(o)  for o where label(o) = L(S)
W_total = Σ weight(o)  for all o ∈ O(S)

Metrics:
    n_observations          = |O(S)|
    n_sources               = |{ source(o) }|
    n_methods               = |{ method(o) : method(o) ≠ "" }|
    n_labels                = |{ label(o) }|
    n_publications          = |{ pub_id(o) : doi(o)≠"" or pmid(o)≠"" }|
    conflict_score          = 1 - n_win / n_observations
                              (count fraction of observations opposing consensus)
    consensus_support_weight = W_win
    opposing_support_weight  = W_total - W_win
    confidence_margin       = w_1st - w_2nd  (label weights sorted desc)
                              if n_labels = 1: w_1st (no runner-up; unopposed support)
    normalized_margin       = (W_win - W_runner_up) / W_total  if W_total > 0
                              else 0.0
                              W_runner_up = second-largest cumulative label weight,
                              or 0 if n_labels = 1
    method_diversity        = 1 - Σ (count_m / n_m)²   (Simpson's diversity)
                              over non-empty methods;  0.0 if ≤1 distinct method
    provenance_richness     = log2(1+n_sources) + log2(1+n_methods) + log2(1+n_publications)
"""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sha8(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def pub_id_from_obs(attrs: dict) -> str:
    """Reconstruct publication node ID from observation attributes."""
    doi  = attrs.get("doi",  "").strip()
    pmid = attrs.get("pmid", "").strip()
    if doi:
        return f"pub:doi:{sha8(doi)}"
    if pmid:
        return f"pub:pmid:{sha8(pmid)}"
    return ""


def parse_weight(v) -> float:
    if v is None or str(v).strip() in ("", "nan"):
        return 0.0
    try:
        return float(v)
    except (ValueError, TypeError):
        return 0.0


# ---------------------------------------------------------------------------
# Load graph
# ---------------------------------------------------------------------------

def load_graph(nodes_path: str, edges_path: str):
    """
    Returns:
        nodes   : dict  node_id → (node_type, attrs_dict)
        seq_ids : set   of SequenceInstance node IDs
        seq_to_obs : dict  seq_id → list of obs_ids (OBSERVED_AS edges)
    """
    nodes: dict[str, tuple[str, dict]] = {}
    with open(nodes_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            nodes[row["node_id"]] = (
                row["node_type"],
                json.loads(row["attributes_json"]),
            )

    seq_ids = {nid for nid, (nt, _) in nodes.items() if nt == "SequenceInstance"}

    seq_to_obs: dict[str, list[str]] = defaultdict(list)
    with open(edges_path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["relation_type"] == "OBSERVED_AS":
                target = row["target"]
                if target in seq_ids:
                    seq_to_obs[target].append(row["source"])

    return nodes, seq_ids, seq_to_obs


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def compute_metrics(seq_id: str, obs_ids: list[str],
                    nodes: dict, consensus_label: str) -> dict:
    """Compute all metrics for one SequenceInstance."""

    observations = []
    for oid in obs_ids:
        if oid not in nodes:
            continue
        _, attrs = nodes[oid]
        observations.append(attrs)

    n = len(observations)
    if n == 0:
        # Should never happen given validation above, but guard defensively
        return {k: 0 for k in (
            "n_observations n_sources n_methods n_labels n_publications "
            "conflict_score consensus_support_weight opposing_support_weight "
            "confidence_margin normalized_margin method_diversity "
            "provenance_richness"
        ).split()}

    # --- Raw collections ---
    sources      = set()
    methods      = set()
    pub_ids      = set()
    label_counts = defaultdict(int)    # label → observation count
    label_weights = defaultdict(float) # label → cumulative weight

    for a in observations:
        sources.add(a.get("source_db", ""))
        m = a.get("method_universal", "").strip()
        if m:
            methods.add(m)
        pid = pub_id_from_obs(a)
        if pid:
            pub_ids.add(pid)
        lbl = a.get("experimental_label", "")
        w   = parse_weight(a.get("evidence_weight", ""))
        label_counts[lbl]  += 1
        label_weights[lbl] += w

    sources.discard("")   # remove blank source if present

    n_observations = n
    n_sources      = len(sources)
    n_methods      = len(methods)
    n_labels       = len(label_counts)
    n_publications = len(pub_ids)

    # --- conflict_score ---
    # Count fraction of observations that disagree with consensus label
    n_win = label_counts.get(consensus_label, 0)
    conflict_score = round(1.0 - n_win / n, 6) if n > 0 else 0.0

    # --- consensus_support_weight / opposing_support_weight ---
    W_win   = label_weights.get(consensus_label, 0.0)
    W_total = sum(label_weights.values())
    W_opp   = W_total - W_win
    consensus_support_weight = round(W_win,   6)
    opposing_support_weight  = round(W_opp,   6)

    # --- confidence_margin ---
    # Sort labels by cumulative weight descending
    sorted_lw = sorted(label_weights.values(), reverse=True)
    if n_labels == 1:
        # Unopposed: no runner-up exists; margin = full support weight
        confidence_margin = round(sorted_lw[0], 6)
    else:
        # Competitive: winner minus runner-up
        confidence_margin = round(sorted_lw[0] - sorted_lw[1], 6)

    # --- normalized_margin ---
    # (W_win - W_runner_up) / W_total
    W_runner_up = sorted_lw[1] if len(sorted_lw) >= 2 else 0.0
    if W_total > 0.0:
        normalized_margin = round((W_win - W_runner_up) / W_total, 6)
    else:
        normalized_margin = 0.0

    # --- method_diversity (Simpson's diversity index) ---
    # Computed over observations with non-empty method_universal
    obs_with_method = [a for a in observations
                       if a.get("method_universal", "").strip()]
    n_m = len(obs_with_method)
    if n_m <= 1:
        method_diversity = 0.0
    else:
        method_freq: dict[str, int] = defaultdict(int)
        for a in obs_with_method:
            method_freq[a["method_universal"].strip()] += 1
        D = sum((c / n_m) ** 2 for c in method_freq.values())
        method_diversity = round(1.0 - D, 6)

    # --- provenance_richness ---
    # log2(1+n_sources) + log2(1+n_methods) + log2(1+n_publications)
    provenance_richness = round(
        math.log2(1 + n_sources)
        + math.log2(1 + n_methods)
        + math.log2(1 + n_publications),
        6,
    )

    return {
        "n_observations":           n_observations,
        "n_sources":                n_sources,
        "n_methods":                n_methods,
        "n_labels":                 n_labels,
        "n_publications":           n_publications,
        "conflict_score":           conflict_score,
        "consensus_support_weight": consensus_support_weight,
        "opposing_support_weight":  opposing_support_weight,
        "confidence_margin":        confidence_margin,
        "normalized_margin":        normalized_margin,
        "method_diversity":         method_diversity,
        "provenance_richness":      provenance_richness,
    }


# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

METRIC_COLS = [
    "seq_id",
    "sequence",
    "uniprot_id",
    "region_start",
    "region_end",
    "is_amyloid",
    "experimental_label",
    "n_observations",
    "n_sources",
    "n_methods",
    "n_labels",
    "n_publications",
    "conflict_score",
    "consensus_support_weight",
    "opposing_support_weight",
    "confidence_margin",
    "normalized_margin",
    "method_diversity",
    "provenance_richness",
]

METRIC_KEYS = [c for c in METRIC_COLS if c not in (
    "seq_id", "sequence", "uniprot_id", "region_start",
    "region_end", "is_amyloid", "experimental_label",
)]


def write_metrics_csv(rows: list[dict], path: str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=METRIC_COLS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Metrics: {len(rows):,} rows → {path}")


def write_updated_nodes(nodes: dict, metrics_by_seq: dict[str, dict], path: str):
    """
    Re-write graph_nodes.csv with metrics merged into SequenceInstance attrs.
    Other node types are written unchanged.
    """
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "node_type", "attributes_json"])
        for node_id, (node_type, attrs) in sorted(nodes.items()):
            if node_type == "SequenceInstance" and node_id in metrics_by_seq:
                merged = {**attrs, **metrics_by_seq[node_id]}
            else:
                merged = attrs
            writer.writerow([node_id, node_type,
                              json.dumps(merged, ensure_ascii=False)])
    n_si = sum(1 for _, (nt, _) in nodes.items() if nt == "SequenceInstance")
    print(f"  Updated graph_nodes.csv: {len(nodes):,} nodes "
          f"({n_si:,} SequenceInstances enriched) → {path}")


# ---------------------------------------------------------------------------
# Summary stats (printed to stdout)
# ---------------------------------------------------------------------------

def print_summary(rows: list[dict]):
    import statistics

    def stats(col):
        vals = [float(r[col]) for r in rows if r[col] != ""]
        if not vals:
            return "n/a"
        return (f"min={min(vals):.3f}  mean={statistics.mean(vals):.3f}  "
                f"median={statistics.median(vals):.3f}  max={max(vals):.3f}")

    print("\n  Metric summary (over all SequenceInstances):")
    for col in METRIC_KEYS:
        print(f"    {col:<30s}  {stats(col)}")

    # Count non-zero conflict scores
    n_conflict = sum(1 for r in rows if float(r["conflict_score"]) > 0)
    print(f"\n  Entries with conflict_score > 0: {n_conflict}")
    n_multi_src = sum(1 for r in rows if int(r["n_sources"]) > 1)
    print(f"  Entries with n_sources > 1:      {n_multi_src}")
    n_multi_meth = sum(1 for r in rows if int(r["n_methods"]) > 1)
    print(f"  Entries with n_methods > 1:      {n_multi_meth}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    script_dir = Path(__file__).parent

    parser = argparse.ArgumentParser(
        description="Compute sequence-level uncertainty metrics for Amyloid DB graph."
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
    parser.add_argument(
        "--no-update-nodes",
        action="store_true",
        help="Skip updating graph_nodes.csv (only write sequence_metrics.csv)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out)

    # Load
    print(f"Loading graph from {args.nodes} ...")
    nodes, seq_ids, seq_to_obs = load_graph(args.nodes, args.edges)
    print(f"  {len(nodes):,} nodes, {len(seq_ids):,} SequenceInstances")

    # Compute
    print("Computing metrics ...")
    rows = []
    metrics_by_seq: dict[str, dict] = {}

    for seq_id in sorted(seq_ids):
        _, seq_attrs = nodes[seq_id]
        consensus_label = seq_attrs.get("experimental_label", "")
        obs_ids = seq_to_obs.get(seq_id, [])

        m = compute_metrics(seq_id, obs_ids, nodes, consensus_label)
        metrics_by_seq[seq_id] = m

        row = {
            "seq_id":             seq_id,
            "sequence":           seq_attrs.get("sequence", ""),
            "uniprot_id":         seq_attrs.get("uniprot_id", ""),
            "region_start":       seq_attrs.get("region_start", ""),
            "region_end":         seq_attrs.get("region_end", ""),
            "is_amyloid":         seq_attrs.get("is_amyloid", ""),
            "experimental_label": consensus_label,
            **m,
        }
        rows.append(row)

    print_summary(rows)

    # Write outputs
    print("\nWriting outputs ...")
    write_metrics_csv(rows, str(out_dir / "sequence_metrics.csv"))

    if not args.no_update_nodes:
        write_updated_nodes(nodes, metrics_by_seq, str(out_dir / "graph_nodes.csv"))

    print("Done.")


if __name__ == "__main__":
    main()
