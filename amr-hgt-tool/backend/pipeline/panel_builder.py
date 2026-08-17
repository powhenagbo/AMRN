"""
pipeline/panel_builder.py
==========================
Builds (or rebuilds) the reference panel used by wrapper.py's
apply_panel_labels() — i.e. produces transmission_classification.csv and
cooccurrence_classification.csv from a folder of genome FASTAs, the same
way the AMR Pipeline Runbook does by hand, but triggerable from the app
as a background job.

Reuses:
  - wrapper.run_rgi() / wrapper.run_amrfinder()   (per-genome gene calling)
  - amr_transmission.run_transmission_analysis()   (Mantel classification)
  - amr_cooccurrence.run_cooccurrence_analysis()   (co-occurrence scatter)

This does NOT touch the single-upload path in wrapper.py at all — it's a
separate, much longer-running job that only produces the panel artifacts
wrapper.py reads from disk.
"""

from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

pipeline_dir = os.path.dirname(os.path.abspath(__file__))
if pipeline_dir not in sys.path:
    sys.path.insert(0, pipeline_dir)

from wrapper import run_rgi, run_amrfinder  # noqa: reuse existing detector wrappers
import amr_transmission
import amr_cooccurrence

PANEL_ROOT = os.path.join(os.path.dirname(pipeline_dir), "panel")
FASTA_EXTS = {".fasta", ".fa", ".fna", ".fas"}


# ── K-mer distance matrix (Bray-Curtis), same math as run_novel_amr.py ────

def kmer_profile(fasta_path: str, k: int) -> dict:
    profile = defaultdict(int)
    with open(fasta_path) as f:
        seq_chunks = []
        for line in f:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip().upper())
        seq = "".join(seq_chunks)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" not in kmer:
            profile[kmer] += 1
    return profile


def bray_curtis_distance(p1: dict, p2: dict) -> float:
    all_kmers = set(p1) | set(p2)
    num = sum(abs(p1.get(k, 0) - p2.get(k, 0)) for k in all_kmers)
    den = sum(p1.get(k, 0) + p2.get(k, 0) for k in all_kmers)
    return num / den if den > 0 else 1.0


def build_distance_matrix(genome_paths: list, k: int, work_dir: str, report=None) -> pd.DataFrame:
    """
    Builds the pairwise Bray-Curtis distance matrix.

    Memory architecture note (real fix, not theoretical): the original
    version of this function held every genome's full k-mer profile
    dictionary in memory simultaneously before computing any distances —
    at k=21 on ~5Mb bacterial genomes, a single profile can have millions
    of entries; holding 200+ of them at once genuinely exhausted memory
    and crashed the whole backend process on a real 226-genome panel
    build (confirmed twice, including once with the checkpoint/resume
    path, which made things worse by re-pickling that entire in-memory
    dict to disk every 200 pairs).

    Fixed by storing each genome's profile as its own file on disk
    (work_dir/kmer_profiles/<name>.pkl) and loading only the two profiles
    needed for each pairwise comparison, with a small LRU cache since
    itertools.combinations() visits pairs in an order where consecutive
    pairs often share one genome — bounds memory to a handful of
    profiles at a time instead of all 226+.
    """
    import pickle
    from functools import lru_cache

    names = [Path(p).stem for p in genome_paths]
    n = len(names)
    profiles_dir = os.path.join(work_dir, "kmer_profiles")
    os.makedirs(profiles_dir, exist_ok=True)
    checkpoint_file = os.path.join(work_dir, "checkpoint_distances.pkl")

    def profile_path(name):
        return os.path.join(profiles_dir, f"{name}.pkl")

    # Compute (or skip, if already on disk from a prior interrupted run —
    # this stage is now resumable too, which it wasn't before) each
    # genome's profile once, immediately writing it to its own file
    # rather than accumulating them all in one in-memory dict.
    if report:
        report(f"Profiling k-mers (k={k}) for {n} genomes")
    for i, gp in enumerate(genome_paths):
        pp = profile_path(names[i])
        if not os.path.exists(pp):
            profile = kmer_profile(gp, k)
            with open(pp, "wb") as f:
                pickle.dump(profile, f)
        if report and (i + 1) % 10 == 0:
            report(f"Profiled {i + 1}/{n} genomes")

    @lru_cache(maxsize=16)
    def load_profile(name):
        with open(profile_path(name), "rb") as f:
            return pickle.load(f)

    if os.path.exists(checkpoint_file):
        if report:
            report("Resuming distance matrix from checkpoint")
        with open(checkpoint_file, "rb") as f:
            state = pickle.load(f)
        matrix = state["matrix"]
        completed = state["completed"]
    else:
        matrix = np.zeros((n, n))
        completed = set()

    all_pairs = list(combinations(range(n), 2))
    remaining = [(i, j) for i, j in all_pairs if (i, j) not in completed]

    if report:
        report(f"Computing {len(remaining)} of {len(all_pairs)} pairwise distances")

    for idx, (i, j) in enumerate(remaining):
        d = bray_curtis_distance(load_profile(names[i]), load_profile(names[j]))
        matrix[i][j] = d
        matrix[j][i] = d
        completed.add((i, j))

        if (idx + 1) % 200 == 0:
            # Only matrix + completed pairs — profiles already live safely
            # on disk as individual files, no need to duplicate them here.
            with open(checkpoint_file, "wb") as f:
                pickle.dump({"matrix": matrix, "completed": completed}, f)
            if report:
                report(f"Checkpoint saved — {idx + 1}/{len(remaining)} distances this run")

    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    return pd.DataFrame(matrix, index=names, columns=names)


# ── AMR presence/absence matrix across the panel ──────────────────────────

def build_amr_matrix(genome_paths: list, workdir: str, detector: str = "rgi", report=None):
    """
    Runs the chosen detector on every genome in the panel and pivots the
    results into a (sample x gene) presence/absence matrix, matching the
    shape amr_transmission.py and amr_cooccurrence.py expect.
    """
    detect_fn = run_rgi if detector == "rgi" else run_amrfinder

    all_rows = []
    failed = []
    for idx, gp in enumerate(genome_paths):
        sample = Path(gp).stem
        if report:
            report(f"Detecting AMR genes ({detector}) — genome {idx + 1}/{len(genome_paths)}")
        genome_workdir = os.path.join(workdir, "amr_raw", sample)
        os.makedirs(genome_workdir, exist_ok=True)
        try:
            df = detect_fn(gp, genome_workdir)
            if not df.empty:
                df = df.copy()
                df["sample"] = sample
                all_rows.append(df)
        except Exception as e:
            failed.append(sample)
            print(f"  [WARN] {sample}: {e}")

    if failed:
        print(f"[WARN] AMR detection failed for {len(failed)} genomes: {failed}")
    if not all_rows:
        raise RuntimeError("No AMR results generated for any panel genome.")

    combined = pd.concat(all_rows, ignore_index=True)
    matrix = combined.groupby(["sample", "gene"]).size().unstack(fill_value=0)
    matrix = (matrix > 0).astype(int)
    drug_class_map = combined.drop_duplicates("gene").set_index("gene")["drug_class"].to_dict()
    return matrix, drug_class_map


# ── Orchestration ──────────────────────────────────────────────────────────

def run_panel_build(
    genome_dir: str,
    version: str = "v1",
    kmer: int = 21,
    detector: str = "rgi",
    permutations: int = 999,
    clonal_r: float = 0.4,
    hgt_r: float = 0.2,
    phi_threshold: float = 0.3,
    scatter_threshold: float = 0.7,
    taxon_scope: str = "Enterobacteriaceae",
    progress_cb=None,
) -> dict:
    """
    Full panel build: distance matrix -> AMR matrix -> Mantel transmission
    classification -> co-occurrence classification -> write into
    backend/panel/<version>/.

    This can take a long time (30-60+ min for dozens of genomes, per the
    AMR Pipeline Runbook) — meant to be run as a background job, same as
    a single-genome upload, just much longer.
    """
    def report(stage):
        if progress_cb:
            progress_cb(stage)

    genome_paths = sorted([
        str(p) for p in Path(genome_dir).glob("*")
        if p.suffix.lower() in FASTA_EXTS
    ])
    if len(genome_paths) < 3:
        raise ValueError(
            f"Need at least 3 genomes to build a panel, found {len(genome_paths)} in {genome_dir}"
        )

    panel_dir = os.path.join(PANEL_ROOT, version)
    os.makedirs(panel_dir, exist_ok=True)
    work_dir = os.path.join(panel_dir, "_build_work")
    os.makedirs(work_dir, exist_ok=True)

    report(f"Found {len(genome_paths)} genomes")

    # Resume support: if the distance matrix already finished in a prior
    # (interrupted) run of this same version, skip recomputing it entirely.
    dist_csv = os.path.join(work_dir, "distance_matrix.csv")
    if os.path.exists(dist_csv):
        report("Distance matrix already complete — loading from disk")
        dist_df = pd.read_csv(dist_csv, index_col=0)
    else:
        dist_df = build_distance_matrix(genome_paths, k=kmer, work_dir=work_dir, report=report)
        dist_df.to_csv(dist_csv)

    # Same idea for the AMR matrix — but note build_amr_matrix() also
    # caches per-genome (via run_rgi/run_amrfinder's own existence checks),
    # so even a partial prior run benefits without needing the whole
    # matrix to have finished.
    amr_csv = os.path.join(work_dir, "amr_presence_absence.csv")
    if os.path.exists(amr_csv):
        report("AMR matrix already complete — loading from disk")
        amr_matrix = pd.read_csv(amr_csv, index_col=0)
        drug_class_map = {}  # not persisted separately; only needed for co-occurrence plot coloring
    else:
        amr_matrix, drug_class_map = build_amr_matrix(genome_paths, work_dir, detector=detector, report=report)
        amr_matrix.to_csv(amr_csv)
    report(f"AMR matrix ready — {amr_matrix.shape[1]} unique genes across {amr_matrix.shape[0]} isolates")

    report("Running per-gene Mantel transmission classification")
    transmission_dir = os.path.join(work_dir, "transmission")
    transmission_df = amr_transmission.run_transmission_analysis(
        dist_df, amr_matrix, transmission_dir,
        permutations=permutations, clonal_r=clonal_r, hgt_r=hgt_r,
    )

    report("Running co-occurrence scatter classification")
    cooccurrence_dir = os.path.join(work_dir, "cooccurrence")
    _, cooccurrence_df = amr_cooccurrence.run_cooccurrence_analysis(
        dist_df, amr_matrix, drug_class_map, cooccurrence_dir,
        phi_threshold=phi_threshold, scatter_threshold=scatter_threshold,
    )

    report("Publishing panel artifacts")
    # wrapper.py's load_panel_gene_labels() expects these two files directly
    # inside panel/<version>/, not in the nested work subfolders.
    transmission_src = os.path.join(transmission_dir, "transmission_classification.csv")
    if os.path.exists(transmission_src):
        pd.read_csv(transmission_src).to_csv(
            os.path.join(panel_dir, "transmission_classification.csv"), index=False
        )

    cooccurrence_src = os.path.join(cooccurrence_dir, "cooccurrence_classification.csv")
    if os.path.exists(cooccurrence_src):
        pd.read_csv(cooccurrence_src).to_csv(
            os.path.join(panel_dir, "cooccurrence_classification.csv"), index=False
        )

    meta = {
        "version": version,
        "n_genomes": len(genome_paths),
        "taxon_scope": taxon_scope,
        "detector": detector,
        "kmer": kmer,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "genome_dir": genome_dir,
    }
    with open(os.path.join(panel_dir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    n_classified = 0 if transmission_df is None or transmission_df.empty else len(transmission_df)
    n_pairs = 0 if cooccurrence_df is None or cooccurrence_df.empty else len(cooccurrence_df)

    return {
        "panel_built": True,
        "version": version,
        "n_genomes": len(genome_paths),
        "n_genes_classified": n_classified,
        "n_cooccurrence_pairs": n_pairs,
        "panel_dir": panel_dir,
    }
