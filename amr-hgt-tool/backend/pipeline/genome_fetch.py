"""
pipeline/genome_fetch.py
==========================
Downloads N genomes for a given taxon using NCBI's official `datasets`
CLI, and flattens them into a single folder ready to hand straight to
panel_builder.run_panel_build() (or the AMR Pipeline Runbook's manual
steps, if preferred).

Requires the `datasets` CLI to be installed on the host machine:

    macOS (Homebrew):
        brew install ncbi-datasets-cli

    Or via conda:
        conda install -c conda-forge ncbi-datasets-cli


"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path


def _slug(text: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]+", "_", text).strip("_")


def _check_datasets_cli():
    if shutil.which("datasets") is None:
        raise RuntimeError(
            "NCBI 'datasets' CLI not found on PATH. Install it first:\n"
            "  macOS (Homebrew): brew install ncbi-datasets-cli\n"
            "  or via conda:     conda install -c conda-forge ncbi-datasets-cli"
        )


def _list_accessions(taxon: str, limit: int, assembly_level: str) -> list[dict]:
    """
    Runs `datasets summary genome taxon` and returns a list of
    {"accession": ..., "organism": ...} dicts for the top N assemblies.
    """
    cmd = [
        "datasets", "summary", "genome", "taxon", taxon,
        "--assembly-level", assembly_level,
        "--limit", str(limit),
        "--as-json-lines",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"'datasets summary' failed (exit {result.returncode}).\n"
            f"--- stderr ---\n{result.stderr[-2000:]}"
        )

    entries = []
    for line in result.stdout.strip().splitlines():
        if not line.strip():
            continue
        obj = json.loads(line)
        accession = obj.get("accession")
        organism = (
            obj.get("organism", {}).get("organism_name")
            if isinstance(obj.get("organism"), dict)
            else None
        )
        if accession:
            entries.append({"accession": accession, "organism": organism or "unknown"})

    return entries


def _download_accessions(accessions: list[str], outdir: str) -> str:
    """Downloads genome data for the given accessions as a single zip, returns its path."""
    zip_path = os.path.join(outdir, "ncbi_dataset.zip")
    cmd = [
        "datasets", "download", "genome", "accession", *accessions,
        "--include", "genome",
        "--filename", zip_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"'datasets download' failed (exit {result.returncode}).\n"
            f"--- stderr ---\n{result.stderr[-2000:]}"
        )
    return zip_path


def _flatten_fastas(zip_path: str, entries: list[dict], outdir: str) -> list[str]:
    """
    Extracts the zip and copies each genome's .fna file up into outdir
    directly (panel_builder.py globs outdir non-recursively), renamed to
    <accession>_<organism>.fasta for readability.
    """
    extract_dir = os.path.join(outdir, "_extracted")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(extract_dir)

    accession_to_organism = {e["accession"]: e["organism"] for e in entries}
    data_dir = os.path.join(extract_dir, "ncbi_dataset", "data")

    saved_paths = []
    if os.path.isdir(data_dir):
        for accession_dir in sorted(os.listdir(data_dir)):
            src_dir = os.path.join(data_dir, accession_dir)
            if not os.path.isdir(src_dir):
                continue
            fna_files = [f for f in os.listdir(src_dir) if f.endswith(".fna")]
            if not fna_files:
                continue
            organism = accession_to_organism.get(accession_dir, "unknown")
            dest_name = f"{accession_dir}_{_slug(organism)}.fasta"
            dest_path = os.path.join(outdir, dest_name)
            shutil.copyfile(os.path.join(src_dir, fna_files[0]), dest_path)
            saved_paths.append(dest_path)

    shutil.rmtree(extract_dir, ignore_errors=True)
    os.remove(zip_path)
    return saved_paths


def fetch_genomes(
    taxon: str,
    limit: int,
    outdir: str,
    assembly_level: str = "complete",
    progress_cb=None,
) -> dict:
    """
    Downloads up to `limit` genome assemblies for `taxon` into `outdir`
    as flat FASTA files, ready to pass as the `directory` argument to
    panel_builder.run_panel_build().

    assembly_level: one of "complete", "chromosome", "scaffold", "contig"
                     (NCBI datasets CLI values; "complete" is the safest
                     default for building a clean reference panel).
    """
    def report(stage):
        if progress_cb:
            progress_cb(stage)

    _check_datasets_cli()
    os.makedirs(outdir, exist_ok=True)

    report(f"Looking up up to {limit} '{taxon}' assemblies on NCBI")
    entries = _list_accessions(taxon, limit, assembly_level)
    if not entries:
        raise RuntimeError(
            f"No assemblies found for taxon '{taxon}' at assembly-level '{assembly_level}'. "
            "Try a broader taxon name or a lower assembly-level (e.g. 'contig')."
        )

    accessions = [e["accession"] for e in entries]
    report(f"Downloading {len(accessions)} genomes")
    zip_path = _download_accessions(accessions, outdir)

    report("Extracting and organizing FASTA files")
    saved_paths = _flatten_fastas(zip_path, entries, outdir)

    return {
        "taxon": taxon,
        "assembly_level": assembly_level,
        "requested": limit,
        "downloaded": len(saved_paths),
        "genome_dir": outdir,
        "accessions": accessions,
    }
