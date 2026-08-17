"""
pipeline/pockets.py
=====================
On-demand binding-pocket detection for a single protein structure, using
fpocket (open-source CLI) run against the AlphaFold-predicted PDB file.

Requires fpocket installed on the host:
    macOS (Homebrew): brew install fpocket
    or via conda:     conda install -c bioconda fpocket

Never run automatically — triggered only when the user clicks "Detect
binding pockets" after a structure has already been fetched via
alphafold.py, same on-demand philosophy as the rest of the structural
analysis features.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

import requests

REQUEST_TIMEOUT = 20


def _check_fpocket():
    if shutil.which("fpocket") is None:
        raise RuntimeError(
            "fpocket not found on PATH. Install it first:\n"
            "  macOS (Homebrew): brew install fpocket\n"
            "  or via conda:     conda install -c bioconda fpocket"
        )


def detect_pockets(pdb_url: str, max_pockets: int = 5) -> list[dict]:
    """
    Downloads the given PDB file and runs fpocket against it, returning
    up to max_pockets ranked by druggability score (descending).

    Each pocket includes real "center" and "size" (docking box dimensions)
    parsed from fpocket's per-pocket alpha-sphere atom file — not a
    volume-based guess. This matters: docking against an approximated,
    origin-defaulted box can silently dock into empty space for any
    structure whose coordinates aren't centered near (0,0,0), producing
    nonsensical (positive/clashing) affinity scores that look like real
    results but aren't.
    """
    _check_fpocket()

    resp = requests.get(pdb_url, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, "structure.pdb")
        with open(pdb_path, "wb") as f:
            f.write(resp.content)

        result = subprocess.run(
            ["fpocket", "-f", pdb_path],
            capture_output=True, text=True, cwd=tmpdir,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"fpocket failed (exit {result.returncode}).\n"
                f"--- stderr (tail) ---\n{result.stderr[-1500:]}"
            )

        out_dir = os.path.join(tmpdir, "structure_out")
        info_path = os.path.join(out_dir, "structure_info.txt")
        if not os.path.exists(info_path):
            return []

        pockets = _parse_fpocket_info(info_path)

        pockets_dir = os.path.join(out_dir, "pockets")
        for pocket in pockets:
            atm_path = os.path.join(pockets_dir, f"pocket{pocket['pocket_id']}_atm.pdb")
            box = _pocket_box_from_atm_file(atm_path)
            pocket["center"] = box["center"]
            pocket["size"] = box["size"]

        pockets.sort(key=lambda p: p.get("druggability_score") or 0, reverse=True)
        return pockets[:max_pockets]


def _pocket_box_from_atm_file(atm_path: str, padding: float = 5.0) -> dict:
    """
    Parses a fpocket pocketN_atm.pdb file's ATOM coordinates (the real
    alpha-spheres defining this specific pocket) and returns a docking
    box: center = centroid, size = bounding-box extent + padding on each
    axis. Returns center=(0,0,0), size=(20,20,20) as a last-resort
    fallback only if the file is missing/unparseable — should be rare.
    """
    xs, ys, zs = [], [], []
    try:
        with open(atm_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    # Standard PDB fixed-width columns for x/y/z.
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))
    except (FileNotFoundError, ValueError):
        pass

    if not xs:
        return {"center": (0.0, 0.0, 0.0), "size": (20.0, 20.0, 20.0)}

    center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    size = (
        max(xs) - min(xs) + padding * 2,
        max(ys) - min(ys) + padding * 2,
        max(zs) - min(zs) + padding * 2,
    )
    # Keep the box in a sane docking range even for a tiny or huge pocket.
    size = tuple(max(12.0, min(35.0, s)) for s in size)
    return {"center": center, "size": size}


def _parse_fpocket_info(info_path: str) -> list[dict]:
    """
    Parses fpocket's <name>_info.txt summary file. Format is repeated
    blocks like:

        Pocket 1 :
            Score :                     0.734
            Druggability Score :        0.845
            Number of Alpha Spheres :   45
            Volume :                    612.3
            ...
    """
    with open(info_path) as f:
        content = f.read()

    pockets = []
    blocks = re.split(r"\n(?=Pocket \d+ :)", content)
    for block in blocks:
        header = re.match(r"Pocket (\d+) :", block)
        if not header:
            continue
        pocket_id = int(header.group(1))

        def _grab(label):
            m = re.search(rf"{re.escape(label)}\s*:\s*([\-0-9.]+)", block)
            return float(m.group(1)) if m else None

        pockets.append({
            "pocket_id": pocket_id,
            "druggability_score": _grab("Druggability Score"),
            "score": _grab("Score"),
            "volume": _grab("Volume"),
            "n_alpha_spheres": _grab("Number of Alpha Spheres"),
        })

    return pockets
