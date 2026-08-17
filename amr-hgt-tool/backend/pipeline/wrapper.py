"""
pipeline/wrapper.py
====================
Adapts existing standalone scripts (kali_hgt.py, AMRFinderPlus subprocess,
amr_island_overlap.py logic, and the reference-panel lookup) into plain
callables the Flask job runner can invoke in sequence.

Drop-in requirement: copy (or symlink) the following files from your
existing KALI/AMR codebase into this `pipeline/` directory before running:

    kali_hgt.py
    amr_island_overlap.py   (only needed if you want its exact enrichment
                              logic re-used; a light overlap check is
                              reimplemented below so this module works
                              standalone in the meantime)

Nothing here touches the existing KALI Flask app or its job store.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import numpy as np

PANEL_ROOT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "panel")
DEFAULT_PANEL_VERSION = "v1"


def resolve_panel_dir(detector: str | None = None) -> str:
    """
    Picks which panel folder to read from. If a detector-specific panel
    exists (e.g. panel/v1_amrfinder, built via /panel/build with
    detector="amrfinder"), prefer it — its gene names will exactly match
    that detector's output, no name-normalization guessing needed. Falls
    back to the shared default panel (panel/v1) if no detector-specific
    one has been built yet.
    """
    if detector:
        detector_specific = os.path.join(PANEL_ROOT, f"{DEFAULT_PANEL_VERSION}_{detector}")
        if os.path.exists(os.path.join(detector_specific, "transmission_classification.csv")):
            return detector_specific
    return os.path.join(PANEL_ROOT, DEFAULT_PANEL_VERSION)


# Kept for any code still referencing the old constant directly — resolves
# to the shared default panel with no detector preference.
PANEL_DIR = os.path.join(PANEL_ROOT, DEFAULT_PANEL_VERSION)


# ── Step 1: AMR gene detection (AMRFinderPlus) ────────────────────────────

def run_amrfinder(fasta_path: str, workdir: str, organism: str | None = None) -> pd.DataFrame:
    """
    Run AMRFinderPlus on a single genome FASTA and return a normalized
    DataFrame: gene, drug_class, contig, start, stop, pct_identity.

    Requires `amrfinder` on PATH and its database installed
    (amrfinder_update run at least once on the host).
    """
    out_tsv = os.path.join(workdir, "amrfinder_out.tsv")

    # Skip re-running AMRFinderPlus if this genome was already processed —
    # same resumability benefit as RGI's caching, for panel builds.
    if os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0:
        raw = pd.read_csv(out_tsv, sep="\t")
        return _parse_amrfinder_output(raw)

    cmd = ["amrfinder", "--nucleotide", fasta_path, "--output", out_tsv, "--threads", "4"]
    if organism:
        cmd += ["--organism", organism]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"AMRFinderPlus failed (exit {result.returncode}).\n"
            f"--- stdout (tail) ---\n{result.stdout[-2000:]}\n"
            f"--- stderr (tail) ---\n{result.stderr[-2000:]}"
        )

    if not os.path.exists(out_tsv) or os.path.getsize(out_tsv) == 0:
        return pd.DataFrame(columns=["gene", "drug_class", "contig", "start", "stop", "pct_identity"])

    raw = pd.read_csv(out_tsv, sep="\t")
    return _parse_amrfinder_output(raw)


def _parse_amrfinder_output(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Shared column-normalization logic for AMRFinderPlus's raw TSV output —
    used by both a fresh run and a cached one, so the two paths can't
    silently drift apart. Column names vary slightly across AMRFinderPlus
    versions — tries current standard names first, with known alternates
    as fallback.
    """
    rename_candidates = {
        "gene": ["Gene symbol", "Element symbol", "Symbol"],
        "drug_class": ["Class", "Element type"],
        "contig": ["Contig id", "Contig"],
        "start": ["Start"],
        "stop": ["Stop"],
        "pct_identity": ["% Identity to reference sequence", "% Identity"],
        # "POINT" here flags a curated point-mutation-based resistance call
        # (as opposed to "AMR" for acquired-gene presence/absence) —
        # AMRFinderPlus's closest equivalent to RGI's model_type distinction.
        "model_type": ["Element subtype"],
        "gene_family": ["Subclass"],
    }
    rename_map = {}
    for target, candidates in rename_candidates.items():
        for c in candidates:
            if c in raw.columns:
                rename_map[c] = target
                break

    df = raw.rename(columns=rename_map)

    if "gene" not in df.columns:
        raise RuntimeError(
            "AMRFinderPlus ran, but the output's gene-name column couldn't be "
            "identified. This usually means the installed AMRFinderPlus version "
            "uses different column headers than expected.\n"
            f"Actual columns returned: {list(raw.columns)}"
        )

    keep = ["gene", "drug_class", "contig", "start", "stop", "pct_identity",
            "model_type", "gene_family"]
    return df[[c for c in keep if c in df.columns]].copy()


# ── Step 1b: AMR gene detection (CARD-RGI via Docker) ─────────────────────

def run_rgi(fasta_path: str, workdir: str) -> pd.DataFrame:
    """
    Run CARD-RGI on a single genome FASTA via Docker, matching the
    invocation in the AMR Pipeline Runbook (finlaymaguire/rgi:latest).
    Returns a normalized DataFrame: gene, drug_class, contig, start, stop,
    pct_identity — same shape as run_amrfinder() so downstream steps don't
    care which detector produced it.

    Requires Docker running on the host and the image already pulled:
        docker pull finlaymaguire/rgi:latest

    Note: RGI's CARD-derived gene names are what the reference panel's
    Mantel labels (transmission_classification.csv) are keyed on if the
    panel was built per the runbook. Prefer this detector when panel
    lookup accuracy matters more than raw speed.
    """
    fasta_path = os.path.abspath(fasta_path)
    genome_dir = os.path.dirname(fasta_path)
    genome_file = os.path.basename(fasta_path)
    name = Path(genome_file).stem

    out_dir = os.path.join(workdir, "rgi_out")
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = f"rgi_{name}"
    out_txt = os.path.join(out_dir, f"{out_prefix}.txt")

    # Skip re-running RGI/Docker if this genome was already processed —
    # matters most for panel builds, where genome_workdir is deterministic
    # per sample and a resumed build shouldn't redo finished genomes.
    if os.path.exists(out_txt) and os.path.getsize(out_txt) > 0:
        raw = pd.read_csv(out_txt, sep="\t")
        return _parse_rgi_output(raw)

    cmd = [
        "docker", "run", "--rm", "--platform", "linux/amd64",
        "-v", f"{genome_dir}:/genomes",
        "-v", f"{out_dir}:/output",
        "finlaymaguire/rgi:latest",
        "rgi", "main",
        "--input_sequence", f"/genomes/{genome_file}",
        "--output_file", f"/output/{out_prefix}",
        "--input_type", "contig",
        "--alignment_tool", "BLAST",
        "--clean",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"RGI/Docker failed (exit {result.returncode}).\n"
            f"--- stdout (tail) ---\n{result.stdout[-2000:]}\n"
            f"--- stderr (tail) ---\n{result.stderr[-2000:]}"
        )

    if not os.path.exists(out_txt) or os.path.getsize(out_txt) == 0:
        return pd.DataFrame(columns=["gene", "drug_class", "contig", "start", "stop", "pct_identity"])

    raw = pd.read_csv(out_txt, sep="\t")
    return _parse_rgi_output(raw)


def _parse_rgi_output(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Shared column-normalization logic for RGI's raw TSV output — used by
    both a fresh RGI run and a cached one (see run_rgi's existence check),
    so the two paths can never silently drift apart.
    """
    rename_map = {}
    for target, candidates in {
        "gene": ["Best_Hit_ARO"],
        "drug_class": ["Drug Class"],
        "contig": ["Contig"],
        "start": ["Start"],
        "stop": ["Stop"],
        "pct_identity": ["Best_Identities"],
        # RGI already computes these — the pipeline was discarding them.
        # SNPs_in_Best_Hit_ARO: specific point mutations found for CARD
        # "protein variant model" hits (e.g. gyrA S83L for fluoroquinolones).
        "mutations": ["SNPs_in_Best_Hit_ARO"],
        # "protein homolog model" = gene presence/absence matters;
        # "protein variant model" = specific mutations matter, not just presence.
        "model_type": ["Model_type"],
        "mechanism": ["Resistance Mechanism"],
        "gene_family": ["AMR Gene Family"],
        # The ARO accession number itself (e.g. "3003378") — distinct from
        # "Best_Hit_ARO" which is the gene NAME despite the similar column
        # name. This is what lets ontology.py look up the term's real
        # position in CARD's hierarchical ontology.
        "aro_id": ["ARO"],
        # RGI already predicts/translates this — needed for real k-mer-
        # based phylogenetic placement among related genes (see
        # phylo_placement.py), reusing the same alignment-free Bray-Curtis
        # approach as the rest of this tool, rather than a separate
        # sequence-fetching step.
        "protein_seq": ["Predicted_Protein"],
    }.items():
        for c in candidates:
            if c in raw.columns:
                rename_map[c] = target
                break

    df = raw.rename(columns=rename_map)

    if "gene" not in df.columns:
        raise RuntimeError(
            "RGI ran, but the output's gene-name column couldn't be identified.\n"
            f"Actual columns returned: {list(raw.columns)}"
        )

    keep = ["gene", "drug_class", "contig", "start", "stop", "pct_identity",
            "mutations", "model_type", "mechanism", "gene_family", "aro_id", "protein_seq"]
    return df[[c for c in keep if c in df.columns]].copy()


# ── Step 2: Compositional anomaly islands (kali_hgt.py) ───────────────────

def run_kali_islands(fasta_path: str, k_list=(3, 4, 5), bin_size=5000, z_threshold=3.0) -> pd.DataFrame:
    """
    Run kali_hgt.py's detect_islands() across all contigs in the FASTA.
    Returns a DataFrame of islands: contig, start, end, size_bp, max_zscore.

    Import is local so the rest of this module still loads even if
    kali_hgt.py hasn't been copied into pipeline/ yet.
    """
    # kali_hgt.py lives alongside this file, but Python doesn't search a
    # module's own directory for plain "import X" unless it's on sys.path —
    # so add pipeline/ explicitly before importing.
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    if pipeline_dir not in sys.path:
        sys.path.insert(0, pipeline_dir)

    from kali_hgt import read_fasta, kmer_vector, gc_content  # noqa: local drop-in file
    from scipy.spatial import distance as _scipy_distance

    contigs = read_fasta(fasta_path)

    # ── Build bins across ALL contigs first, before any scoring ──
    # This is the key fix: the original per-contig approach scored each
    # contig against its OWN mean composition, which fails in two ways —
    # (1) contigs shorter than bin_size*3 were skipped outright, and more
    # importantly (2) a contig that IS itself a whole horizontally-acquired
    # element (common in fragmented draft assemblies, where assemblers
    # can't merge divergent inserted DNA into the main chromosome) has
    # nothing anomalous-looking *within itself* to contrast against — the
    # anomaly only shows up relative to the genome's true chromosomal
    # background. So: one shared background, computed across every bin
    # from every contig combined, then every bin scored against that.
    bin_records = []  # each: {"contig", "start", "end", "seq"}
    for contig_id, seq in contigs.items():
        n = len(seq)
        if n < bin_size:
            continue  # a contig shorter than one bin has nothing to score
        starts = list(range(0, n - bin_size + 1, bin_size))
        for s in starts:
            bin_records.append({"contig": contig_id, "start": s, "end": s + bin_size - 1, "seq": seq[s:s + bin_size]})

    if len(bin_records) < 3:
        # Not enough bins genome-wide to establish a meaningful background —
        # e.g. an extremely small or heavily fragmented assembly.
        return pd.DataFrame(columns=["contig", "start", "end", "size_bp", "max_zscore"])

    # Feature matrix across all bins, all contigs, concatenated across k values.
    feature_cols = []
    for k in k_list:
        vecs = np.array([kmer_vector(b["seq"], k) for b in bin_records])
        feature_cols.append(vecs)
    X = np.hstack(feature_cols)

    background = X.mean(axis=0)  # genome-wide composition — the correct reference
    scores = np.array([
        _scipy_distance.cosine(X[i], background) if np.linalg.norm(X[i]) > 0 else 0.0
        for i in range(len(bin_records))
    ])
    z_scores = (scores - scores.mean()) / scores.std() if scores.std() > 0 else np.zeros(len(scores))

    for i, b in enumerate(bin_records):
        b["score"] = scores[i]
        b["zscore"] = z_scores[i]
        b["flagged"] = z_scores[i] > z_threshold

    # Merge adjacent flagged bins into islands — adjacency only within the
    # same contig, since "adjacent" across two different contigs is meaningless.
    all_islands = []
    for contig_id in contigs:
        contig_bins = [b for b in bin_records if b["contig"] == contig_id]
        in_island, island_start, island_bins = False, None, []
        for b in contig_bins:
            if b["flagged"]:
                if not in_island:
                    in_island, island_start, island_bins = True, b["start"], []
                island_bins.append(b)
            else:
                if in_island:
                    all_islands.append(_finalize_island(contig_id, island_start, island_bins))
                    in_island = False
        if in_island and island_bins:
            all_islands.append(_finalize_island(contig_id, island_start, island_bins))

    if not all_islands:
        return pd.DataFrame(columns=["contig", "start", "end", "size_bp", "max_zscore"])
    return pd.DataFrame(all_islands)


def _finalize_island(contig_id: str, start: int, island_bins: list) -> dict:
    end = island_bins[-1]["end"]
    return {
        "contig": contig_id,
        "start": start,
        "end": end,
        "size_bp": end - start + 1,
        "max_zscore": max(b["zscore"] for b in island_bins),
    }


# ── Step 3: Overlap check — does each AMR gene fall inside an island? ─────

def overlap_amr_with_islands(amr_df: pd.DataFrame, islands_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each AMR gene hit, flag whether its coordinates fall inside any
    compositional anomaly island on the same contig, and attach that
    island's max Z-score as a continuous anomaly signal.
    """
    if amr_df.empty:
        return amr_df.assign(in_island=False, island_zscore=None)

    records = []
    for _, gene in amr_df.iterrows():
        match = islands_df[
            (islands_df.get("contig") == gene.get("contig"))
            & (islands_df.get("start") <= gene.get("start", -1))
            & (islands_df.get("end") >= gene.get("stop", -1))
        ] if not islands_df.empty else pd.DataFrame()

        in_island = not match.empty
        zscore = float(match["max_zscore"].max()) if in_island else None

        row = gene.to_dict()
        row["in_island"] = in_island
        row["island_zscore"] = zscore
        records.append(row)

    return pd.DataFrame(records)


# ── Step 4: Panel lookup — apply precomputed Mantel labels if available ───

def _normalize_gene_name(name: str) -> str:
    """
    Reduces a gene name to a canonical form so detectors with different
    naming conventions can be matched against the same panel — e.g.
    AMRFinderPlus's 'blaSHV-11' and RGI's 'SHV-11' both normalize to
    'shv11'. Strips common resistance-gene prefixes, case, and punctuation.
    This is intentionally aggressive (case-insensitive, no separators) —
    exact matching is still tried first in apply_panel_labels(), so this
    only kicks in as a fallback and won't override a precise hit.
    """
    if not isinstance(name, str):
        return ""
    n = name.strip().lower()
    # Common acquired beta-lactamase prefix that AMRFinderPlus includes
    # and CARD/RGI typically omits (e.g. "blaSHV-11" vs "SHV-11").
    if n.startswith("bla") and len(n) > 3:
        n = n[3:]
    # Drop everything but letters/digits so punctuation differences
    # (hyphens, underscores, apostrophes, parentheses) don't block a match.
    n = re.sub(r"[^a-z0-9]", "", n)
    return n


def debug_panel_state(detector: str | None = None) -> dict:
    """
    Diagnostic helper — reports exactly what this running process sees
    when it tries to load the panel for a given detector, so mismatches
    between 'the file on disk' and 'what the server is actually reading'
    are immediately visible instead of silently falling back to
    island-only evidence.
    """
    panel_dir = resolve_panel_dir(detector)
    csv_path = os.path.join(panel_dir, "transmission_classification.csv")
    exists = os.path.exists(csv_path)
    labels = load_panel_gene_labels(detector)
    return {
        "wrapper_file": os.path.abspath(__file__),
        "detector_requested": detector,
        "panel_dir_used": os.path.abspath(panel_dir),
        "is_detector_specific_panel": panel_dir != os.path.join(PANEL_ROOT, DEFAULT_PANEL_VERSION),
        "csv_path": os.path.abspath(csv_path),
        "csv_exists": exists,
        "csv_size_bytes": os.path.getsize(csv_path) if exists else 0,
        "n_exact_labels_loaded": len(labels["exact"]),
        "has_CRP": "CRP" in labels["exact"],
        "sample_gene_names": list(labels["exact"].keys())[:10],
    }


def load_panel_meta(detector: str | None = None) -> dict:
    panel_dir = resolve_panel_dir(detector)
    meta_path = os.path.join(panel_dir, "meta.json")
    if not os.path.exists(meta_path):
        return {"version": None, "n_genomes": 0, "taxon_scope": None}
    with open(meta_path) as f:
        return json.load(f)


def load_panel_gene_labels(detector: str | None = None) -> dict:
    """
    Returns {gene_name: {"classification": ..., "mantel_r": ..., "p_value": ...}}
    from transmission_classification.csv, keyed both by the exact gene name
    as it appears in the panel AND by a normalized form (see
    _normalize_gene_name) so lookups can fall back to a fuzzy match when
    the live detector uses different naming conventions than whatever
    built the panel.

    If a detector-specific panel exists (see resolve_panel_dir), its gene
    names should already match the live detector's output exactly, so the
    normalized fallback should rarely be needed in that case.
    """
    panel_dir = resolve_panel_dir(detector)
    csv_path = os.path.join(panel_dir, "transmission_classification.csv")
    if not os.path.exists(csv_path):
        return {"exact": {}, "normalized": {}}

    df = pd.read_csv(csv_path)
    exact = {}
    normalized = {}
    for _, row in df.iterrows():
        entry = {
            "classification": row["classification"],
            "mantel_r": row["mantel_r"],
            "p_value": row["p_value"],
            "n_carriers": row.get("n_carriers"),
            "panel_gene_name": row["gene"],
        }
        exact[row["gene"]] = entry
        norm_key = _normalize_gene_name(row["gene"])
        # If two panel gene names collide after normalization (rare), keep
        # the first — exact matching is always tried first anyway, so this
        # only affects the fallback path.
        normalized.setdefault(norm_key, entry)

    return {"exact": exact, "normalized": normalized}


def apply_panel_labels(merged_df: pd.DataFrame, detector: str | None = None) -> pd.DataFrame:
    """
    Attach panel-derived classification where available; otherwise fall
    back to island-overlap-only evidence and mark confidence accordingly.

    Tries an exact gene-name match first (highest confidence — same
    detector/naming convention as the panel), then a normalized match
    (handles cases like AMRFinderPlus's 'blaSHV-11' vs RGI's 'SHV-11'),
    before giving up and falling back to island-overlap-only evidence.

    detector selects which panel to read from (see resolve_panel_dir) —
    if a detector-specific panel exists, its names should already match
    exactly, so the normalized fallback becomes a rare safety net rather
    than the primary matching path.
    """
    labels = load_panel_gene_labels(detector)
    exact_labels = labels["exact"]
    normalized_labels = labels["normalized"]

    if merged_df.empty:
        # No AMR genes detected in this genome at all — a legitimate,
        # if uncommon, real result, not an error. Without this guard,
        # calling .apply(..., result_type="expand") on an empty frame
        # below produces a shape mismatch, and the resulting
        # results[0]/[1]/[2] column lookups raise a bare KeyError(0) —
        # which renders as the unhelpful literal text "0" wherever the
        # exception message gets surfaced. Return early with the right
        # columns already present instead.
        return merged_df.assign(
            classification=pd.Series(dtype=object),
            mantel_r=pd.Series(dtype=float),
            evidence_source=pd.Series(dtype=object),
        )

    def classify(row):
        gene_name = row.get("gene")

        panel_hit = exact_labels.get(gene_name)
        if panel_hit:
            return (panel_hit["classification"], panel_hit["mantel_r"], "panel",
                    panel_hit.get("p_value"), panel_hit.get("n_carriers"))

        norm_hit = normalized_labels.get(_normalize_gene_name(gene_name))
        if norm_hit:
            return (norm_hit["classification"], norm_hit["mantel_r"], "panel_normalized",
                    norm_hit.get("p_value"), norm_hit.get("n_carriers"))

        if row.get("in_island"):
            return "HGT", None, "island_only", None, None
        return "CLONAL", None, "island_only", None, None

    results = merged_df.apply(classify, axis=1, result_type="expand")
    merged_df = merged_df.copy()
    merged_df["classification"] = results[0]
    merged_df["mantel_r"] = results[1]
    merged_df["evidence_source"] = results[2]
    merged_df["mantel_p_value"] = results[3]
    merged_df["mantel_n_carriers"] = results[4]
    return merged_df


# ── Orchestration ──────────────────────────────────────────────────────────

DETECTORS = {
    "rgi": {
        "label": "CARD-RGI",
        "run": lambda fasta_path, workdir, organism: run_rgi(fasta_path, workdir),
        "note": "Matches the detector the reference panel was built with — "
                "panel Mantel labels apply directly.",
    },
    "amrfinder": {
        "label": "AMRFinderPlus",
        "run": lambda fasta_path, workdir, organism: run_amrfinder(fasta_path, workdir, organism),
        "note": "Faster, NCBI-curated, but gene names may not match the RGI-based "
                "panel labels — expect more genes to fall back to island-only evidence.",
    },
}


def _sanitize(obj):
    """
    Recursively replace NaN/Infinity with None so the result is valid JSON.
    Python's json module happily emits a literal `NaN` token (a non-standard
    extension), but JavaScript's JSON.parse() correctly rejects it — this
    is what breaks the frontend if any pandas/numpy value slips through
    as NaN (e.g. a gene with no computed Mantel r).
    """
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    if isinstance(obj, dict):
        return {k: _sanitize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v) for v in obj]
    return obj


def _normalize_contig_ids(df: pd.DataFrame, fasta_path: str) -> pd.DataFrame:
    """
    RGI's 'Contig' output column isn't the genome's real contig name — it's
    RGI's own internal Prodigal-based ORF identifier, formatted as
    '<real_contig_name>_<ORF_number>' (e.g. 'NZ_CP044158.1_131' means ORF
    #131 on contig 'NZ_CP044158.1', not a contig literally named that).
    This has been silently breaking every AMR-gene-to-island overlap
    check, since the island detector reports real contig names and this
    never matched them.

    Fix: cross-reference against the genome's actual contig names (read
    from the FASTA itself) rather than guessing with a blind regex strip —
    a real contig name could legitimately end in digits/underscores, so
    matching against ground truth avoids false corrections.
    """
    if "contig" not in df.columns or df.empty:
        return df

    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    if pipeline_dir not in sys.path:
        sys.path.insert(0, pipeline_dir)
    from kali_hgt import read_fasta

    real_contigs = set(read_fasta(fasta_path).keys())

    def resolve(raw_contig):
        if pd.isna(raw_contig):
            return raw_contig
        raw_contig = str(raw_contig)
        if raw_contig in real_contigs:
            return raw_contig  # already correct (e.g. AMRFinderPlus, which doesn't have this issue)
        # Progressively strip trailing "_<segment>" pieces until we find
        # a real contig name, or give up and leave it as-is.
        candidate = raw_contig
        for _ in range(3):  # a handful of attempts is plenty — real contig names aren't usually this deep
            if "_" not in candidate:
                break
            candidate = candidate.rsplit("_", 1)[0]
            if candidate in real_contigs:
                return candidate
        return raw_contig  # no match found — leave unchanged rather than guess wrong

    df = df.copy()
    df["contig"] = df["contig"].apply(resolve)
    return df


def run_full_pipeline(fasta_path: str, workdir: str, organism: str | None = None,
                       detector: str = "rgi", progress_cb=None) -> dict:
    """
    Runs the full single-upload pipeline and returns a JSON-serializable
    result matching the shape the frontend dashboard expects.

    detector: "rgi" (default) or "amrfinder" — caller's choice, made once
    per upload. Both produce the same downstream DataFrame shape, so
    nothing past Step 1 needs to know which one ran.

    progress_cb(str) is called with human-readable stage names, so the
    Flask route can forward them over SSE.
    """
    if detector not in DETECTORS:
        raise ValueError(f"Unknown detector '{detector}'. Choose one of: {list(DETECTORS)}")

    def report(stage):
        if progress_cb:
            progress_cb(stage)

    chosen = DETECTORS[detector]
    report(f"Running {chosen['label']}")
    amr_df = chosen["run"](fasta_path, workdir, organism)
    amr_df = _normalize_contig_ids(amr_df, fasta_path)

    report("Scoring compositional islands")
    islands_df = run_kali_islands(fasta_path)

    report("Checking island overlap")
    merged_df = overlap_amr_with_islands(amr_df, islands_df)

    report("Applying reference panel labels")
    final_df = apply_panel_labels(merged_df, detector=detector)

    panel_meta = load_panel_meta(detector=detector)

    genes = []
    for _, row in final_df.iterrows():
        klass = row.get("classification")
        model_type = row.get("model_type") if pd.notna(row.get("model_type")) else None

        # Validation finding (checked against real panel data, not
        # speculative): low Mantel r has two distinct biological causes —
        # true horizontal transfer, or a point mutation that arose
        # independently multiple times across unrelated lineages
        # (convergent evolution / homoplasy) under shared selective
        # pressure — e.g. gyrA/parC fluoroquinolone-resistance mutations,
        # which are essentially never horizontally transferred as a whole
        # gene in Enterobacteriaceae. Both look identical to a presence/
        # absence Mantel test, but mean very different things
        # epidemiologically. Flag this ambiguity whenever it applies,
        # using model_type (RGI's "protein variant model" / AMRFinder's
        # "POINT") — data already captured, just not cross-checked until now.
        interpretation_caveat = None
        if klass == "HGT" and model_type and (
            model_type == "POINT" or "variant" in model_type.lower()
        ):
            interpretation_caveat = (
                "This gene is a point-mutation-based resistance mechanism, not an "
                "acquired gene — its 'HGT' call (low correlation with phylogeny) may "
                "reflect the same mutation arising independently in multiple unrelated "
                "lineages (convergent evolution), rather than true horizontal transfer "
                "of a mobile element. Interpret with this distinction in mind."
            )

        genes.append({
            "id": row.get("gene"),
            "drug": row.get("drug_class"),
            "contig": row.get("contig"),
            "pos": int(row.get("start", 0)) if pd.notna(row.get("start")) else None,
            "stop": int(row.get("stop", 0)) if pd.notna(row.get("stop")) else None,
            "identity": row.get("pct_identity") if pd.notna(row.get("pct_identity")) else None,
            "score": row.get("island_zscore"),
            "mantel_r": row.get("mantel_r"),
            "mantel_p_value": row.get("mantel_p_value") if pd.notna(row.get("mantel_p_value")) else None,
            "mantel_n_carriers": int(row.get("mantel_n_carriers")) if pd.notna(row.get("mantel_n_carriers")) else None,
            "klass": klass,
            "evidence": row.get("evidence_source"),
            # Mutation/mechanism data the detector already computed —
            # previously discarded before it ever reached the frontend.
            "mutations": row.get("mutations") if pd.notna(row.get("mutations")) else None,
            "model_type": model_type,
            "mechanism": row.get("mechanism") if pd.notna(row.get("mechanism")) else None,
            "gene_family": row.get("gene_family") if pd.notna(row.get("gene_family")) else None,
            "interpretation_caveat": interpretation_caveat,
            "protein_seq": row.get("protein_seq") if pd.notna(row.get("protein_seq")) else None,
        })

        aro_id = row.get("aro_id") if pd.notna(row.get("aro_id")) else None
        if aro_id:
            try:
                pipeline_dir = os.path.dirname(os.path.abspath(__file__))
                if pipeline_dir not in sys.path:
                    sys.path.insert(0, pipeline_dir)
                from ontology import format_aro_id, get_lineage, get_ontology_siblings

                genes[-1]["aro_id"] = format_aro_id(aro_id)
                genes[-1]["ontology_lineage"] = get_lineage(aro_id)
                genes[-1]["ontology_siblings"] = get_ontology_siblings(aro_id)
            except Exception as e:
                # Non-fatal by design — this is a new, not-yet-battle-tested
                # feature (first real download of CARD's ontology file,
                # first real parse). A failure here shouldn't break the
                # rest of the analysis; surface it quietly in the gene's
                # own record instead so it's debuggable without crashing.
                raw_id = str(aro_id).strip()
                genes[-1]["aro_id"] = raw_id if raw_id.startswith("ARO:") else f"ARO:{raw_id}"
                genes[-1]["ontology_lineage"] = []
                genes[-1]["ontology_siblings"] = []
                genes[-1]["ontology_error"] = str(e)

    # Real phylogenetic placement (see phylo_placement.py) — grouped by
    # CARD gene family, using protein sequences RGI already computed.
    # Distinct from ontology.py's lineage/siblings: this is an actual
    # quantitative tree built from sequence divergence, not a curated
    # classification. Non-fatal on any failure — a new, first-real-run
    # feature, shouldn't be able to break the rest of the analysis.
    try:
        pipeline_dir = os.path.dirname(os.path.abspath(__file__))
        if pipeline_dir not in sys.path:
            sys.path.insert(0, pipeline_dir)
        from phylo_placement import build_family_tree
        by_family = {}
        for g in genes:
            fam = g.get("gene_family")
            if fam:
                by_family.setdefault(fam, []).append(g)

        for fam, members in by_family.items():
            tree = build_family_tree(members)
            if tree:
                for g in members:
                    g["phylo_tree"] = tree
            else:
                # Distinguish "not enough related genes in THIS genome to
                # build a meaningful tree" from an actual failure — both
                # look identical to the user otherwise (nothing shown).
                usable_count = sum(1 for g in members if g.get("protein_seq"))
                for g in members:
                    g["phylo_tree_status"] = (
                        f"Only {usable_count} gene(s) in this family detected in this "
                        f"genome — need at least 3 to build a meaningful tree."
                    )
    except Exception as e:
        for g in genes:
            g.setdefault("phylo_tree_error", str(e))

    # protein_seq was only needed internally for tree-building — drop it
    # before this reaches the frontend, no reason to ship raw sequences
    # over the wire.
    for g in genes:
        g.pop("protein_seq", None)

    islands_summary = []
    if not islands_df.empty:
        for _, row in islands_df.iterrows():
            islands_summary.append({
                "contig": row.get("contig"),
                "start": int(row.get("start")) if pd.notna(row.get("start")) else None,
                "end": int(row.get("end")) if pd.notna(row.get("end")) else None,
                "size_bp": int(row.get("size_bp")) if pd.notna(row.get("size_bp")) else None,
                "max_zscore": float(row.get("max_zscore")) if pd.notna(row.get("max_zscore")) else None,
            })

    result = {
        "genes": genes,
        "islands": islands_summary,
        "panel": panel_meta,
        "detector": {
            "id": detector,
            "label": chosen["label"],
            "note": chosen["note"],
        },
    }
    return _sanitize(result)
