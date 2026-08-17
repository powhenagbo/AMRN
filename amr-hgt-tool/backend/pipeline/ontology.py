"""
pipeline/ontology.py
======================
Real evolutionary/functional relationships via CARD's own Antibiotic
Resistance Ontology (ARO) — a proper hierarchical ontology (is_a
relationships), not just a flat "same gene family string" grouping.

RGI's raw output already includes each gene's ARO accession number (the
"ARO" column) — this was being discarded before. Given that accession,
this module provides:
    - the term's parent lineage (walking up is_a relationships)
    - ontology siblings (other ARO terms sharing the same immediate parent)

This is a strictly stronger, more precise relationship than the existing
"same gene family, same genome" grouping — it draws on CARD's actual
curated ontology structure, not a coincidental string match, and isn't
limited to genes present in the current genome.

Downloads and caches the ontology file once (backend/ontology_cache/aro.obo).
Requires the `obonet` package: pip install obonet
"""

from __future__ import annotations

import os

import requests

ARO_OBO_URL = "http://purl.obolibrary.org/obo/aro.obo"
CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "ontology_cache")
OBO_PATH = os.path.join(CACHE_DIR, "aro.obo")

REQUEST_TIMEOUT = 30

_graph_cache = None  # module-level cache — parse once per process, not per request


def _ensure_obo_downloaded():
    if os.path.exists(OBO_PATH):
        return
    os.makedirs(CACHE_DIR, exist_ok=True)
    resp = requests.get(ARO_OBO_URL, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    with open(OBO_PATH, "wb") as f:
        f.write(resp.content)


def _load_graph():
    """
    Downloads (if needed) and parses the ARO ontology, caching the parsed
    graph in memory for the lifetime of the process — parsing ~8,500
    terms takes a couple seconds, not worth repeating per request.
    """
    global _graph_cache
    if _graph_cache is not None:
        return _graph_cache

    try:
        import obonet
    except ImportError:
        raise RuntimeError(
            "The 'obonet' package is required for ontology lookups.\n"
            "Install it with: pip install obonet"
        )

    _ensure_obo_downloaded()
    _graph_cache = obonet.read_obo(OBO_PATH)
    return _graph_cache


def format_aro_id(raw_aro: str | int) -> str:
    """
    RGI's raw 'ARO' column is a bare number (e.g. '3003378'). CARD's
    ontology IDs are formatted 'ARO:3003378'. Normalizes either form.
    """
    s = str(raw_aro).strip()
    if s.startswith("ARO:"):
        return s
    return f"ARO:{s}"


def get_term_info(aro_id: str) -> dict | None:
    """Returns {"name": ..., "definition": ...} for a given ARO term, or None if not found."""
    graph = _load_graph()
    aro_id = format_aro_id(aro_id)
    if aro_id not in graph:
        return None
    data = graph.nodes[aro_id]
    return {
        "aro_id": aro_id,
        "name": data.get("name"),
        "definition": (data.get("def") or "").split('"')[1] if data.get("def") else None,
    }


def get_lineage(aro_id: str, max_levels: int = 5) -> list[dict]:
    """
    Walks up is_a relationships from the given term, returning ancestor
    terms from immediate parent up to `max_levels` levels — e.g. for a
    specific beta-lactamase: "class A beta-lactamase" -> "beta-lactamase"
    -> "antibiotic inactivation enzyme" -> ...

    Returns [] if the term isn't found or has no parents recorded.
    """
    graph = _load_graph()
    aro_id = format_aro_id(aro_id)
    if aro_id not in graph:
        return []

    lineage = []
    current = aro_id
    seen = {current}
    for _ in range(max_levels):
        parents = [
            v for u, v, k in graph.out_edges(current, keys=True)
            if k == "is_a"
        ]
        if not parents:
            break
        parent = parents[0]  # ARO terms are usually single-parent in practice
        if parent in seen:
            break  # defensive — avoid any cycle
        info = get_term_info(parent)
        if info:
            lineage.append(info)
        seen.add(parent)
        current = parent

    return lineage


def get_ontology_siblings(aro_id: str, max_siblings: int = 15) -> list[dict]:
    """
    Returns other ARO terms that share the given term's immediate parent —
    true ontological siblings per CARD's curated hierarchy, not limited to
    genes detected in any particular genome. E.g. for one specific
    class A beta-lactamase, this returns other class A beta-lactamases
    CARD recognizes, whether or not they showed up in this analysis.
    """
    graph = _load_graph()
    aro_id = format_aro_id(aro_id)
    if aro_id not in graph:
        return []

    parents = [v for u, v, k in graph.out_edges(aro_id, keys=True) if k == "is_a"]
    if not parents:
        return []
    parent = parents[0]

    siblings = []
    for u, v, k in graph.in_edges(parent, keys=True):
        if k == "is_a" and u != aro_id:
            info = get_term_info(u)
            if info:
                siblings.append(info)
        if len(siblings) >= max_siblings:
            break

    return siblings
