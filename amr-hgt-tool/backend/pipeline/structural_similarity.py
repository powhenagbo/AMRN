"""
pipeline/structural_similarity.py
====================================
On-demand structural similarity search using Foldseek, comparing a
gene's AlphaFold structure against other structures already fetched and
cached locally during this tool's use (backend/structure_cache/) —
NOT the full AlphaFold DB or PDB, which would require a multi-gigabyte
reference database and separate setup.

This answers "which other AMR proteins I've already looked at have a
similar fold" — genuinely useful and achievable locally. Searching the
entire universe of known structures would need Foldseek's full reference
database as a documented future upgrade, not something bolted on here.

Requires Foldseek installed on the host:
    macOS (Homebrew): brew install foldseek
    or via conda:     conda install -c conda-forge -c bioconda foldseek
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "structure_cache")
MANIFEST_PATH = os.path.join(CACHE_DIR, "manifest.json")


def _check_foldseek():
    if shutil.which("foldseek") is None:
        raise RuntimeError(
            "Foldseek not found on PATH. Install it first:\n"
            "  macOS (Homebrew): brew install foldseek\n"
            "  or via conda:     conda install -c conda-forge -c bioconda foldseek"
        )


def load_manifest() -> dict:
    """{accession: {"gene_name": ..., "organism": ..., "cached_at": ...}}"""
    if not os.path.exists(MANIFEST_PATH):
        return {}
    with open(MANIFEST_PATH) as f:
        return json.load(f)


def save_manifest(manifest: dict):
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)


def find_similar_structures(query_accession: str, max_hits: int = 5) -> list[dict]:
    """
    Compares the cached structure for query_accession against every other
    cached structure (from earlier AlphaFold lookups made through this
    tool), using Foldseek's easy-search. Returns hits ranked by TM-score.

    Returns an empty list (not an error) if fewer than 2 structures are
    cached yet — there's genuinely nothing to compare against, which is
    a normal early-usage state, not a failure.
    """
    _check_foldseek()

    query_path = os.path.join(CACHE_DIR, f"{query_accession}.pdb")
    if not os.path.exists(query_path):
        raise RuntimeError(
            f"No cached structure for {query_accession} — run an AlphaFold "
            "lookup for this gene first (this caches its structure "
            "automatically for future similarity comparisons)."
        )

    other_files = [
        f for f in os.listdir(CACHE_DIR)
        if f.endswith(".pdb") and f != f"{query_accession}.pdb"
    ]
    if not other_files:
        return []

    with tempfile.TemporaryDirectory() as tmpdir:
        target_dir = os.path.join(tmpdir, "targets")
        os.makedirs(target_dir)
        for f in other_files:
            shutil.copyfile(os.path.join(CACHE_DIR, f), os.path.join(target_dir, f))

        out_tsv = os.path.join(tmpdir, "results.tsv")
        work_dir = os.path.join(tmpdir, "foldseek_tmp")
        os.makedirs(work_dir)

        cmd = [
            "foldseek", "easy-search", query_path, target_dir, out_tsv, work_dir,
            "--format-output", "query,target,alntmscore,evalue,bits",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Foldseek failed (exit {result.returncode}).\n"
                f"--- stderr (tail) ---\n{result.stderr[-1500:]}"
            )

        if not os.path.exists(out_tsv):
            return []

        hits = []
        with open(out_tsv) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 5:
                    continue
                _query, target, tm_score, evalue, bits = row[:5]
                hits.append({
                    "target_accession": target.replace(".pdb", ""),
                    "tm_score": float(tm_score) if tm_score else None,
                    "evalue": float(evalue) if evalue else None,
                    "bits": float(bits) if bits else None,
                })

    hits.sort(key=lambda h: h.get("tm_score") or 0, reverse=True)
    return hits[:max_hits]
