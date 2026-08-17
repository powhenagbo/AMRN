"""
pipeline/pdb_validation.py
=============================
Real experimental cross-reference — queries RCSB PDB's live public search
API for actual crystal/cryo-EM structures matching a given UniProt
accession, and reports their real bound ligands (if any).

This exists specifically so docking results can be checked against real
experimental data, rather than trusting a Vina score in isolation. It
deliberately does NOT hardcode any specific PDB accession numbers or
ligand names from memory — every result here is a live, verifiable query
against RCSB's actual database, not an assertion I can't independently
confirm.

RCSB Search API docs: https://search.rcsb.org/#search-api
"""

from __future__ import annotations

import requests

SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_API = "https://data.rcsb.org/rest/v1/core/entry"
REQUEST_TIMEOUT = 20


def find_pdb_structures(uniprot_accession: str, max_results: int = 10) -> list[dict]:
    """
    Queries RCSB for PDB entries whose polymer entity is cross-referenced
    to the given UniProt accession — using RCSB's documented attribute
    search on rcsb_polymer_entity_container_identifiers.
    reference_sequence_identifiers.database_accession, not a free-text
    search .

    Returns a list of {"pdb_id": ..., "ligands": [...]} — ligands is the
    list of real bound heteroatom-component names/IDs for that structure,
    fetched from RCSB's data API per matching entry.

    Returns an empty list (not an error) if no structures exist for this
    protein — a normal, common outcome, not every protein has been
    crystallized.
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_accession,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_results}},
    }

    try:
        resp = requests.post(SEARCH_API, json=query, timeout=REQUEST_TIMEOUT)
        if resp.status_code == 204:
            return []  # RCSB's documented response for a query with zero matches — not an error
        if resp.status_code != 200:
            raise RuntimeError(
                f"RCSB search API returned HTTP {resp.status_code}: {resp.text[:300]}"
            )
        data = resp.json()
    except requests.RequestException as e:
        raise RuntimeError(f"RCSB search API request failed: {e}")

    hits = data.get("result_set", [])
    results = []
    for hit in hits:
        pdb_id = hit.get("identifier")
        if not pdb_id:
            continue
        ligands = _fetch_ligands(pdb_id)
        results.append({"pdb_id": pdb_id, "ligands": ligands})

    return results


def _fetch_ligands(pdb_id: str) -> list[str]:
    """
    Fetches the real list of bound non-polymer components (ligands,
    cofactors, ions) for a given PDB entry, filtering out common
    crystallization artifacts (water, common buffer/cryoprotectant
    components) that aren't biologically meaningful bound ligands.
    """
    common_artifacts = {"HOH", "GOL", "EDO", "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA", "K", "ACT", "DMS", "PEG"}
    try:
        resp = requests.get(f"{DATA_API}/{pdb_id}", timeout=REQUEST_TIMEOUT)
        if resp.status_code != 200:
            return []
        data = resp.json()
    except (requests.RequestException, ValueError):
        return []

    components = data.get("rcsb_entry_info", {}).get("nonpolymer_bound_components") or []
    return [c for c in components if c not in common_artifacts]
