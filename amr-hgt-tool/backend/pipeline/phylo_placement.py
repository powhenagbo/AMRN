"""
pipeline/phylo_placement.py
=============================
Real phylogenetic placement of related AMR genes, using the SAME
alignment-free k-mer distance approach as the rest of this tool (Bray-
Curtis on k-mer profiles + neighbor-joining tree construction) — applied
to protein sequences instead of whole genomes.

This is a genuine, quantitative phylogeny (built from actual sequence
divergence), not the curated classification hierarchy in ontology.py —
the two are complementary, different kinds of "evolutionary relationship":
ontology.py answers "how does CARD's expert curation classify this gene,"
this module answers "how similar are these genes' actual sequences."

Deliberately alignment-free, consistent with this tool's whole design
philosophy (KALI = K-mer Alignment-free Inference) — no MSA, no external
alignment tool, just the same math already used for genome-level distance
matrices, reused at the protein-sequence scale.

Requires: scikit-bio (already a dependency, used for the genome-level NJ
tree in the panel-build pipeline).
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations

import numpy as np


def protein_kmer_profile(seq: str, k: int = 3) -> dict:
    """
    Same idea as kali_hgt.py's DNA k-mer profiling, applied to a protein
    sequence (20-letter amino acid alphabet instead of 4-letter DNA).
    """
    profile = defaultdict(int)
    seq = seq.upper().replace("*", "")  # strip stop-codon markers if present
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "X" not in kmer:  # X = ambiguous/unknown residue
            profile[kmer] += 1
    return profile


def bray_curtis_distance(p1: dict, p2: dict) -> float:
    """Identical formula to the genome-level version elsewhere in this tool."""
    all_kmers = set(p1) | set(p2)
    num = sum(abs(p1.get(k, 0) - p2.get(k, 0)) for k in all_kmers)
    den = sum(p1.get(k, 0) + p2.get(k, 0) for k in all_kmers)
    return num / den if den > 0 else 1.0


def build_family_tree(genes: list[dict], k: int = 3, min_members: int = 3) -> dict | None:
    """
    Given a list of {"id": gene_name, "protein_seq": sequence} for genes
    in the same CARD gene family (from a single genome's analysis), builds
    a real neighbor-joining tree from their protein k-mer distances.

    Returns a nested dict {"name": ..., "length": ..., "children": [...]}
    suitable for direct JSON serialization and frontend rendering, or
    None if there aren't enough members with usable sequences to build a
    meaningful tree (need at least min_members, default 3 — a tree of 2
    is just a single branch, not informative).
    """
    usable = [g for g in genes if g.get("protein_seq") and len(g["protein_seq"]) >= k]
    if len(usable) < min_members:
        return None

    names = [g["id"] for g in usable]
    profiles = [protein_kmer_profile(g["protein_seq"], k) for g in usable]

    n = len(names)
    matrix = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        d = bray_curtis_distance(profiles[i], profiles[j])
        matrix[i][j] = d
        matrix[j][i] = d

    try:
        from skbio import DistanceMatrix
        from skbio.tree import nj
    except ImportError:
        raise RuntimeError("scikit-bio is required for phylogenetic tree construction.")

    # Deduplicate names for skbio (it requires unique tip labels) — genes
    # detected more than once (rare, but possible with overlapping ORF
    # calls) get a disambiguating suffix.
    seen = {}
    unique_names = []
    for name in names:
        if name not in seen:
            seen[name] = 0
            unique_names.append(name)
        else:
            seen[name] += 1
            unique_names.append(f"{name} ({seen[name]})")

    dm = DistanceMatrix(matrix, unique_names)
    tree = nj(dm)

    return _treenode_to_dict(tree)


def _treenode_to_dict(node) -> dict:
    """Converts an skbio TreeNode into a plain nested dict for JSON serialization."""
    result = {
        "name": node.name or "",
        "length": float(node.length) if node.length is not None else 0.0,
    }
    if node.children:
        result["children"] = [_treenode_to_dict(child) for child in node.children]
    return result
