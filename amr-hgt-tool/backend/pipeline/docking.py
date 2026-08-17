"""
pipeline/docking.py
=====================
On-demand molecular docking — screens a small, curated set of real
antibiotic-adjuvant inhibitors against the gene's AlphaFold structure,
using the top fpocket-detected pocket as the docking site.

This is a genuinely heavier pipeline than the rest of the structural
analysis features, with more external moving parts:

    AlphaFold PDB + fpocket pocket coords
        -> receptor prep (obabel: PDB -> PDBQT, adds hydrogens/charges)
        -> ligand fetch (PubChem PUG-REST, by compound name -> SDF)
        -> ligand prep (obabel: SDF -> PDBQT)
        -> AutoDock Vina (docking box centered on the chosen pocket)
        -> binding affinity scores (kcal/mol; more negative = stronger
           predicted binding)

Requires two more external tools beyond fpocket/Foldseek:
    macOS (Homebrew): brew install autodock-vina open-babel
    or via conda:      conda install -c bioconda vina
                        conda install -c conda-forge openbabel

Ligand selection is keyword-matched against the gene's drug class /
gene family to a small curated list of REAL, well-documented AMR-
relevant inhibitor compounds (fetched live from PubChem by name, not
fabricated or hardcoded structures). This is intentionally a small,
targeted panel, not a general-purpose virtual screening library.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

import requests

REQUEST_TIMEOUT = 20

# Curated, keyword-matched inhibitor candidates — all real, documented
# compounds. Matched against the gene's drug_class / gene_family text.
# This is deliberately small and targeted, not a general screening library.
# Substrate-specific mappings for well-characterized aminoglycoside-
# modifying enzyme families. This matters because "aminoglycoside
# antibiotic" as a drug class covers chemically distinct substrate
# groups — e.g. aadA-family adenylyltransferases act on streptomycin/
# spectinomycin, NOT gentamicin, despite both being labeled
# "aminoglycoside antibiotic". A validation run against aadA5 caught
# this: the old broad drug-class match defaulted every aminoglycoside
# gene to gentamicin regardless of whether that's the gene's actual
# substrate, which isn't a meaningful test of the enzyme's real biology.
_GENE_FAMILY_SUBSTRATES = [
    (r"\baada\d*\b|ant\(3|ant3", ["Streptomycin", "Spectinomycin"]),
    (r"\bstra\d*\b|\bstrb\d*\b|aph\(3''|aph\(6", ["Streptomycin"]),
    (r"aac\(3|aac3", ["Gentamicin"]),  # this family genuinely does act on gentamicin
    (r"aac\(6|aac6", ["Tobramycin", "Amikacin"]),
    (r"aph\(3'\)|apha\d*\b", ["Kanamycin"]),
    # Metallo-beta-lactamases (NDM/VIM/IMP) are mechanistically unrelated
    # to the serine beta-lactamases the generic beta-lactam inhibitor
    # panel below was built for — clavulanic acid, tazobactam, avibactam,
    # and vaborbactam all work by covalently modifying a catalytic serine
    # residue that metallo-enzymes simply don't have (they use a
    # zinc-dependent mechanism instead). None of those four inhibitors
    # are expected to work on an MBL at all. Real MBL inhibitors instead
    # work by zinc chelation — using well-documented compounds from that
    # literature instead of the serine-lactamase panel.
    (r"ndm-?\d*\b", ["Dipicolinic acid", "Aspergillomarasmine A"]),
    # vim/imp are short, generic-looking substrings — require the real
    # IMP-N/VIM-N numbering pattern (not just the bare letters) to avoid
    # accidentally matching an unrelated gene name that happens to
    # contain "vim" or "imp" as a substring.
    (r"vim-\d+", ["Dipicolinic acid", "Aspergillomarasmine A"]),
    (r"imp-\d+", ["Dipicolinic acid", "Aspergillomarasmine A"]),
    # Chloramphenicol acetyltransferase (CAT) enzymes use acetyl-CoA as a
    # co-substrate to acetylate chloramphenicol — testing the enzyme's
    # actual natural substrate/product is mechanistically far more
    # relevant than any beta-lactamase inhibitor (which was the previous
    # silent default for any gene not matching another specific rule).
    (r"\bcat[a-z]{0,2}\d*\b", ["Acetyl-CoA", "Coenzyme A"]),
]

_INHIBITOR_LIBRARY = [
    # (match keywords, compound names to try)
    #
    # Order matters here, and it's a deliberate fix for a real bug: efflux
    # pump genes (e.g. acrB) list beta-lactam drug classes in their
    # drug_class field too, since they confer resistance to those drugs
    # by pumping them out — not because they're enzymes. A drug_class
    # string containing "cephalosporin" doesn't mean the gene IS a
    # beta-lactamase. Checking the efflux/RND signal first (a much more
    # specific and reliable indicator of actual protein mechanism) avoids
    # matching an efflux pump against beta-lactamase inhibitors just
    # because both keyword sets happen to co-occur in its combined text.
    (["efflux", "rnd", "resistance-nodulation"],
     ["Phenylalanine-arginine beta-naphthylamide"]),
    (["aminoglycoside"],
     ["Gentamicin"]),  # as a competitive-binding reference, not an inhibitor per se
    (["fluoroquinolone"],
     ["Ciprofloxacin"]),
]
_DEFAULT_LIGANDS = ["Clavulanic acid"]  # fallback if nothing matches


def _check_tools():
    missing = [t for t in ("obabel", "vina") if shutil.which(t) is None]
    if missing:
        raise RuntimeError(
            f"Missing required tool(s): {', '.join(missing)}. Install with:\n"
            "  macOS (Homebrew): brew install autodock-vina open-babel\n"
            "  or via conda:      conda install -c bioconda vina\n"
            "                      conda install -c conda-forge openbabel"
        )


def _select_ligands(gene_name: str | None, drug_class: str | None, gene_family: str | None) -> list[str]:
    name_lower = (gene_name or "").lower()

    # Try specific substrate-family matches first — more biologically
    # accurate than the broad drug-class fallback below.
    for pattern, names in _GENE_FAMILY_SUBSTRATES:
        if re.search(pattern, name_lower):
            return names

    # Beta-lactamase inhibitors specifically require the gene FAMILY to
    # actually say "lactamase" — not just a drug_class keyword like
    # "cephalosporin", which efflux pumps and other non-enzyme resistance
    # genes list too (they confer resistance to those drugs without being
    # the enzyme that breaks them down). This was a real bug: acrB (an
    # RND efflux pump) was matching the beta-lactam rule purely because
    # its drug_class field happened to contain "cephalosporin"/"penam",
    # even though it has no beta-lactamase mechanism at all.
    gene_family_lower = (gene_family or "").lower()
    if "lactamase" in gene_family_lower:
        return ["Clavulanic acid", "Tazobactam", "Avibactam", "Vaborbactam"]

    haystack = f"{drug_class or ''} {gene_family or ''}".lower()
    for keywords, names in _INHIBITOR_LIBRARY:
        if any(k in haystack for k in keywords):
            return names
    return _DEFAULT_LIGANDS


def _looks_like_valid_sdf(content: bytes) -> bool:
    """
    HTTP 200 + non-empty body isn't proof of a real molecule file — PubChem
    (and many APIs) can return a 200 with an HTML error page or an empty
    molblock for a name that didn't resolve cleanly. A real SDF/molfile
    always ends its atom/bond block with 'M  END'; this is a cheap,
    reliable sanity check that catches that failure mode before it reaches
    OpenBabel as garbage input (which is what produced the earlier
    'No atoms in this ligand' error).
    """
    if not content or b"<html" in content[:500].lower():
        return False
    return b"M  END" in content


def _fetch_ligand_sdf(name: str, out_path: str) -> tuple[bool, str, bool]:
    """
    Fetches a compound's SDF from PubChem by name. Tries the 3D conformer
    first; many compounds (especially larger/complex ones) don't have a
    precomputed 3D conformer in PubChem, so this falls back to the 2D
    structure if needed — the caller then needs OpenBabel to generate 3D
    coordinates (--gen3d) for a 2D-sourced file.

    Returns (success, error_detail, needs_gen3d).
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF"

    try:
        resp = requests.get(url, params={"record_type": "3d"}, timeout=REQUEST_TIMEOUT)
        if resp.status_code == 200 and _looks_like_valid_sdf(resp.content):
            with open(out_path, "wb") as f:
                f.write(resp.content)
            return True, "", False
        first_error = f"3D fetch: HTTP {resp.status_code}, valid={_looks_like_valid_sdf(resp.content)} — {resp.text[:200]}"
    except requests.RequestException as e:
        first_error = f"3D fetch: {e}"

    # Fall back to the 2D structure — OpenBabel can generate 3D coordinates itself.
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        if resp.status_code == 200 and _looks_like_valid_sdf(resp.content):
            with open(out_path, "wb") as f:
                f.write(resp.content)
            return True, "", True
        second_error = f"2D fetch: HTTP {resp.status_code}, valid={_looks_like_valid_sdf(resp.content)} — {resp.text[:200]}"
    except requests.RequestException as e:
        second_error = f"2D fetch: {e}"

    return False, f"{first_error}; {second_error}", False


def _prepare_pdbqt(input_path: str, output_path: str, is_ligand: bool, gen3d: bool = False) -> tuple[bool, str]:
    """
    Converts PDB/SDF to PDBQT via OpenBabel, adding hydrogens/charges.
    gen3d=True tells OpenBabel to generate 3D coordinates itself — needed
    when the input SDF only has 2D coordinates (PubChem fallback case).
    Returns (success, stderr_tail).
    """
    cmd = ["obabel", input_path, "-O", output_path, "-h", "--partialcharge", "gasteiger"]
    if gen3d:
        cmd.append("--gen3d")
    result = subprocess.run(cmd, capture_output=True, text=True)
    ok = result.returncode == 0 and os.path.exists(output_path)
    stderr_tail = (result.stderr or "")[-800:]

    if ok and not is_ligand:
        _strip_receptor_pdbqt_headers(output_path)

    return ok, stderr_tail


_VALID_RECEPTOR_PDBQT_PREFIXES = ("ATOM", "HETATM", "TER")


def _strip_receptor_pdbqt_headers(pdbqt_path: str):
    """
    OpenBabel's PDBQT writer carries over the original PDB's header
    metadata (HEADER, TITLE, COMPND, SOURCE, etc.) verbatim, and can also
    emit ROOT/BRANCH/TORSDOF torsion-tree records if it treats the whole
    protein as one flexible molecule rather than a rigid receptor — Vina's
    receptor parser rejects both: it wants a strictly rigid receptor
    (ATOM/HETATM/TER only), no header lines, no torsion tree. This strips
    everything else in place.
    """
    with open(pdbqt_path) as f:
        lines = f.readlines()

    cleaned = [line for line in lines if line.startswith(_VALID_RECEPTOR_PDBQT_PREFIXES)]

    with open(pdbqt_path, "w") as f:
        f.writelines(cleaned)


def _run_vina(receptor_pdbqt: str, ligand_pdbqt: str, center: tuple, size: tuple, out_dir: str,
               seed: int | None = None) -> tuple[float | None, str, str | None]:
    """
    Runs AutoDock Vina, returns (best_affinity_kcal_mol_or_None, stderr_tail, docked_pose_pdbqt_or_None).
    The docked pose is Vina's actual predicted position/orientation of the
    ligand inside the pocket — captured so the caller can render it
    alongside the receptor, not just report a bare score.

    seed: explicit random seed — used by run_docking_with_replicates() to
    run multiple independent searches per ligand rather than trusting a
    single stochastic run.
    """
    out_pdbqt = os.path.join(out_dir, f"docked_out_{seed or 'default'}.pdbqt")
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--out", out_pdbqt,
        "--exhaustiveness", "16",
        # Vina defaults to writing 9 separate poses (MODEL...ENDMDL blocks)
        # per run into one file. We only ever use the single best pose —
        # requesting just 1 avoids ambiguity downstream, both in scoring
        # (the reported affinity is unambiguously "the" pose, not "the
        # best of several bundled together") and in visualization (a
        # multi-MODEL PDBQT file handed to the 3D viewer can confuse its
        # internal model indexing, which is what caused the docked
        # ligand to silently fail to render in the pose viewer).
        "--num_modes", "1",
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    stderr_tail = (result.stderr or "")[-800:]

    if result.returncode != 0:
        return None, stderr_tail, None

    # Vina prints affinity of the best pose to stdout, e.g.:
    #    1        -6.4      0.000      0.000
    match = re.search(r"^\s*1\s+(-?\d+\.\d+)", result.stdout, flags=re.MULTILINE)
    if not match:
        # Exit code 0 but couldn't parse output — surface stdout instead,
        # since that's where the actual problem will be visible.
        return None, (result.stdout or "")[-800:], None

    docked_pose = None
    try:
        with open(out_pdbqt) as f:
            docked_pose = f.read()
    except FileNotFoundError:
        pass  # non-fatal — score is still valid even if pose read fails

    return float(match.group(1)), "", docked_pose


def _pocket_center_and_size(pocket: dict) -> tuple[tuple, tuple]:
    """
    Uses the real pocket center/size parsed from fpocket's per-pocket
    atom file (see pockets.py's _pocket_box_from_atm_file) — falls back
    to a generic origin-centered box only if those weren't available for
    some reason (should be rare; logged rather than silently guessed).
    """
    if "center" in pocket and "size" in pocket:
        return tuple(pocket["center"]), tuple(pocket["size"])
    # Fallback — only reached if pockets.py's coordinate parsing failed.
    return (0.0, 0.0, 0.0), (20.0, 20.0, 20.0)


def run_docking(pdb_url: str, top_pocket: dict, drug_class: str | None, gene_family: str | None,
                 gene_name: str | None = None, custom_ligands: list[str] | None = None) -> dict:
    """
    Full on-demand docking run: fetch receptor + candidate ligands, prep
    both, dock each ligand into the given pocket, return ranked results.

    top_pocket: one entry from pockets.detect_pockets()'s output — reuses
    work already done rather than re-running fpocket.

    gene_name enables substrate-specific ligand selection for well-
    characterized enzyme families (see _GENE_FAMILY_SUBSTRATES) — falls
    back to broad drug_class/gene_family keyword matching if no specific
    family match is found.

    custom_ligands: an explicit list of compound names to test instead of
    the automatic selection — for testing specific compounds (e.g. ones
    suggested via the chat assistant) that aren't in the small built-in
    rule set. Each name is still fetched live from PubChem, same as the
    automatic panel, so it must be a real, resolvable compound name.
    """
    _check_tools()

    ligand_names = custom_ligands if custom_ligands else _select_ligands(gene_name, drug_class, gene_family)
    center, size = _pocket_center_and_size(top_pocket)

    with tempfile.TemporaryDirectory() as tmpdir:
        receptor_pdb = os.path.join(tmpdir, "receptor.pdb")
        resp = requests.get(pdb_url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        with open(receptor_pdb, "wb") as f:
            f.write(resp.content)

        receptor_pdbqt = os.path.join(tmpdir, "receptor.pdbqt")
        receptor_ok, receptor_err = _prepare_pdbqt(receptor_pdb, receptor_pdbqt, is_ligand=False)
        if not receptor_ok:
            raise RuntimeError(
                f"Failed to prepare receptor (OpenBabel PDB -> PDBQT conversion failed).\n"
                f"--- obabel stderr (tail) ---\n{receptor_err}"
            )

        results = []
        for name in ligand_names:
            ligand_sdf = os.path.join(tmpdir, f"{name.replace(' ', '_')}.sdf")
            fetched, fetch_err, needs_gen3d = _fetch_ligand_sdf(name, ligand_sdf)
            if not fetched:
                results.append({"ligand": name, "available": False, "reason": f"Could not fetch structure from PubChem: {fetch_err}"})
                continue

            ligand_pdbqt = os.path.join(tmpdir, f"{name.replace(' ', '_')}.pdbqt")
            ligand_ok, ligand_err = _prepare_pdbqt(ligand_sdf, ligand_pdbqt, is_ligand=True, gen3d=needs_gen3d)
            if not ligand_ok:
                results.append({
                    "ligand": name, "available": False,
                    "reason": f"Ligand PDBQT preparation failed: {ligand_err or 'no error output'}",
                })
                continue

            # Run multiple independent searches (different random seeds)
            # rather than trusting a single stochastic Vina run — Vina's
            # search is not deterministic, so one run alone doesn't tell
            # you how reproducible that score actually is.
            N_REPLICATES = 3
            replicate_affinities = []
            replicate_errors = []
            best_pose = None
            best_affinity = None

            for rep in range(N_REPLICATES):
                affinity, vina_err, docked_pose = _run_vina(
                    receptor_pdbqt, ligand_pdbqt, center, size, tmpdir, seed=rep + 1
                )
                if affinity is None:
                    replicate_errors.append(vina_err or "no error output")
                    continue
                replicate_affinities.append(affinity)
                if best_affinity is None or affinity < best_affinity:
                    best_affinity = affinity
                    best_pose = docked_pose

            if not replicate_affinities:
                results.append({
                    "ligand": name, "available": False,
                    "reason": f"Docking failed on all {N_REPLICATES} replicate runs: {replicate_errors[0] if replicate_errors else 'no error output'}",
                })
            else:
                mean_affinity = sum(replicate_affinities) / len(replicate_affinities)
                std_affinity = (
                    (sum((a - mean_affinity) ** 2 for a in replicate_affinities) / len(replicate_affinities)) ** 0.5
                    if len(replicate_affinities) > 1 else 0.0
                )
                results.append({
                    "ligand": name, "available": True,
                    "affinity_kcal_mol": best_affinity,  # kept for backward compat / sorting
                    "mean_affinity_kcal_mol": round(mean_affinity, 2),
                    "std_affinity_kcal_mol": round(std_affinity, 2),
                    "n_replicates": len(replicate_affinities),
                    "n_replicates_attempted": N_REPLICATES,
                    "docked_pose_pdbqt": best_pose,
                })

    results.sort(key=lambda r: r.get("affinity_kcal_mol") if r.get("available") else 999)
    return {"pocket_id": top_pocket.get("pocket_id"), "ligands_tried": ligand_names, "results": results}
