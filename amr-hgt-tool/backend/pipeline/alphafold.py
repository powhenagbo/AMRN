"""
pipeline/alphafold.py
=======================
On-demand structural lookup for a single AMR gene, triggered by the user
clicking "Analyze with AlphaFold" on a specific gene in the dashboard —
never run automatically for every detected gene, to keep normal AMR
analysis fast.

Pipeline:
    gene name (+ optional organism)
        -> UniProt REST search (gene name, reviewed entries preferred)
        -> UniProt accession
        -> AlphaFold DB prediction API (https://alphafold.ebi.ac.uk/api/prediction/{accession})
        -> structure metadata (pLDDT, PDB/CIF download URLs)

Matching is done via UniProt's gene-name search rather than raw protein-
sequence BLAST — AlphaFold DB's public API doesn't expose a documented
sequence-search endpoint; UniProt's REST API + AlphaFold's per-accession
API is the stable, documented path. An organism filter and preference for
"reviewed" (Swiss-Prot curated) entries are used to reduce false matches,
since short AMR gene symbols like 'CRP' or 'SHV-11' are reused across
many organisms.
"""

from __future__ import annotations

import json
import os
import re

import requests

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "structure_cache")
MANIFEST_PATH = os.path.join(CACHE_DIR, "manifest.json")

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

REQUEST_TIMEOUT = 15  # seconds — this is a live, on-demand user action, fail fast rather than hang

# Words that mark the end of the actual gene/protein name in CARD's longer
# descriptive ARO terms, e.g. "AcrAB-TolC with MarR mutations conferring
# resistance to..." — everything from "with"/"conferring"/etc. onward is
# resistance-phenotype description, not part of the searchable name.
_CARD_DESCRIPTOR_STOPWORDS = [
    "with", "conferring", "mutant", "mutants", "mutation", "mutations",
    "variant", "variants", "resistant", "resistance",
]

# Matches a leading two-word organism binomial, e.g. "Escherichia coli ",
# "Klebsiella pneumoniae ", "Haemophilus influenzae ".
_ORGANISM_PREFIX_RE = re.compile(r"^([A-Z][a-z]+ [a-z]+)\s+")


def _parse_card_gene_name(raw_name: str) -> dict:
    """
    CARD/RGI gene names range from clean symbols ('SHV-11', 'oqxA') to long
    descriptive phrases ('Escherichia coli AcrAB-TolC with MarR mutations
    conferring resistance to ciprofloxacin and tetracycline'). UniProt's
    gene-name search needs something close to an actual symbol — this
    extracts the best candidate plus any organism named in the string.

    Returns {"candidate": str, "organism": str | None}. For already-clean
    names, candidate == raw_name unchanged.
    """
    name = raw_name.strip()

    organism_match = _ORGANISM_PREFIX_RE.match(name)
    organism = organism_match.group(1) if organism_match else None
    remainder = name[organism_match.end():] if organism_match else name

    # Truncate at the first descriptor stopword (word-boundary, case-insensitive).
    pattern = r"\b(" + "|".join(_CARD_DESCRIPTOR_STOPWORDS) + r")\b"
    stop_match = re.search(pattern, remainder, flags=re.IGNORECASE)
    candidate = remainder[:stop_match.start()] if stop_match else remainder
    candidate = candidate.strip(" ,;:-")

    return {"candidate": candidate or remainder.strip(), "organism": organism}


def _fetch_uniprot_features(accession: str) -> list[dict]:
    """
    Fetches curated active-site and binding-site residue annotations for
    a UniProt entry — this is standard UniProt curation for well-studied
    enzymes, not a new database or external tool; just a field on an
    entry we already have the accession for.

    Returns a list of {"type", "position", "description"} — empty list
    if the entry has none (common for less-studied proteins) or the
    request fails, never raises.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(
            url, params={"fields": "ft_act_site,ft_binding,ft_site"}, timeout=REQUEST_TIMEOUT
        )
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError):
        return []

    relevant = []
    for feature in data.get("features", []):
        ftype = feature.get("type", "")
        if ftype not in ("Active site", "Binding site", "Site"):
            continue
        location = feature.get("location", {})
        start = location.get("start", {}).get("value")
        end = location.get("end", {}).get("value")
        position = str(start) if start == end else f"{start}-{end}"
        relevant.append({
            "type": ftype,
            "position": position,
            "description": feature.get("description") or "",
        })

    return relevant


def _search_uniprot(gene: str, organism: str | None = None, field: str = "gene") -> dict | None:
    """
    Searches UniProt by a specific field ("gene" or "protein_name"),
    preferring reviewed (Swiss-Prot) entries, optionally filtered by
    organism. Returns the top hit's accession + protein name, or None.
    """
    field_query = f'gene:{gene}' if field == "gene" else f'protein_name:"{gene}"'
    query_parts = [field_query]
    if organism:
        query_parts.append(f'organism_name:"{organism}"')

    for reviewed_filter in ["AND reviewed:true", ""]:
        query = " ".join(query_parts) + (" " + reviewed_filter if reviewed_filter else "")
        params = {
            "query": query.strip(),
            "fields": "accession,protein_name,organism_name,gene_names",
            "format": "json",
            "size": 1,
        }
        try:
            resp = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            data = resp.json()
        except (requests.RequestException, ValueError):
            continue

        results = data.get("results", [])
        if results:
            hit = results[0]
            protein_desc = hit.get("proteinDescription", {})
            recommended = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value")
            return {
                "accession": hit.get("primaryAccession"),
                "protein_name": recommended or gene,
                "organism": hit.get("organism", {}).get("scientificName"),
            }

    return None


def _find_uniprot_entry(candidate: str, organism: str | None) -> tuple[dict | None, list[str]]:
    """
    Tries several search strategies in order of specificity/confidence,
    stopping at the first hit. Returns (result_or_None, attempts_tried) —
    the attempts list is surfaced in the "not available" reason so a
    failed lookup is debuggable instead of a bare "not found".
    """
    attempts = []

    strategies = [
        ("gene name + organism", "gene", organism),
        ("gene name only", "gene", None),
        ("protein name + organism", "protein_name", organism),
        ("protein name only", "protein_name", None),
    ]

    for label, field, org in strategies:
        attempts.append(label)
        hit = _search_uniprot(candidate, organism=org, field=field)
        if hit and hit.get("accession"):
            return hit, attempts

    return None, attempts


def _fetch_alphafold_prediction(accession: str) -> dict | None:
    """
    Queries AlphaFold DB's prediction API for a given UniProt accession.
    Returns structure metadata, or None if no prediction exists for it
    (a real, expected outcome — not every UniProt entry has a prediction,
    e.g. very short peptides or entries excluded from AlphaFold DB).
    """
    url = ALPHAFOLD_API_URL.format(accession=accession)
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        data = resp.json()
    except (requests.RequestException, ValueError):
        return None

    if not data:
        return None

    entry = data[0] if isinstance(data, list) else data
    return {
        "pdb_url": entry.get("pdbUrl"),
        "cif_url": entry.get("cifUrl"),
        "mean_plddt": entry.get("globalMetricValue") or entry.get("confidenceAvgLocalScore"),
        "model_created": entry.get("modelCreatedDate"),
    }


def _cache_structure(accession: str, pdb_url: str, gene_name: str, organism: str | None):
    """
    Downloads and caches the structure locally under structure_cache/, and
    records it in manifest.json (accession -> gene/organism metadata).
    This is what lets structural_similarity.py compare structures against
    each other later without needing a huge external reference database —
    the cache just grows naturally as you look up more genes.

    Never raises — a caching failure shouldn't break the main AlphaFold
    lookup response, it just means this gene won't be available for
    similarity comparisons until a retry succeeds.
    """
    try:
        os.makedirs(CACHE_DIR, exist_ok=True)
        pdb_path = os.path.join(CACHE_DIR, f"{accession}.pdb")
        if not os.path.exists(pdb_path):
            resp = requests.get(pdb_url, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            with open(pdb_path, "wb") as f:
                f.write(resp.content)

        manifest = {}
        if os.path.exists(MANIFEST_PATH):
            with open(MANIFEST_PATH) as f:
                manifest = json.load(f)
        manifest[accession] = {"gene_name": gene_name, "organism": organism}
        with open(MANIFEST_PATH, "w") as f:
            json.dump(manifest, f, indent=2)
    except Exception:
        pass  # non-fatal by design — see docstring


def find_alphafold_structure(gene: str, organism: str | None = None) -> dict:
    """
    Full lookup: gene name -> UniProt accession -> AlphaFold structure.
    Always returns a dict with "available" set — callers should check
    that rather than relying on exceptions, since "no structure found"
    is a normal, expected result for many genes, not an error.

    Handles both clean gene symbols ('SHV-11', 'oqxA') and CARD's longer
    descriptive names ('Escherichia coli AcrAB-TolC with MarR mutations
    conferring resistance to...') by extracting a search-friendly
    candidate and organism from the name first (see _parse_card_gene_name).
    """
    parsed = _parse_card_gene_name(gene)
    candidate = parsed["candidate"]
    # An explicitly-passed organism (e.g. from the caller/UI) wins over
    # one parsed out of the gene name itself.
    resolved_organism = organism or parsed["organism"]

    uniprot_hit, attempts = _find_uniprot_entry(candidate, resolved_organism)
    if not uniprot_hit or not uniprot_hit.get("accession"):
        return {
            "gene": gene,
            "available": False,
            "reason": (
                f"No matching UniProt entry found for '{candidate}'"
                + (f" (organism: {resolved_organism})" if resolved_organism else "")
                + f". Tried: {', '.join(attempts)}."
            ),
            "search_candidate": candidate,
        }

    accession = uniprot_hit["accession"]
    structure = _fetch_alphafold_prediction(accession)
    if not structure:
        return {
            "gene": gene,
            "available": False,
            "reason": "UniProt entry found, but no AlphaFold prediction exists for it",
            "uniprot": accession,
            "protein": uniprot_hit.get("protein_name"),
        }

    active_site_features = _fetch_uniprot_features(accession)
    _cache_structure(accession, structure["pdb_url"], gene, uniprot_hit.get("organism"))

    return {
        "gene": gene,
        "available": True,
        "uniprot": accession,
        "protein": uniprot_hit.get("protein_name"),
        "organism": uniprot_hit.get("organism"),
        "structure_url": structure["pdb_url"],
        "pdb_url": structure["pdb_url"],
        "cif_url": structure["cif_url"],
        "mean_plddt": structure["mean_plddt"],
        "active_site_features": active_site_features,
    }
