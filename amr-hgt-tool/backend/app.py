"""
app.py — AMR-HGT tool backend
==============================
A standalone Flask app, fully separate from the existing KALI web app.
Different codebase, different job store (jobs_store.json, not jobs.json),
different port. Safe to run alongside KALI without touching it.

Run:
    pip install flask flask-cors pandas
    python app.py            # serves on http://localhost:5050

Endpoints:
    POST /upload                  -> {job_id}
    POST /jobs/<job_id>/start     -> {status: "started"}
    GET  /jobs/<job_id>/stream    -> SSE progress stream
    GET  /jobs/<job_id>/results   -> final JSON result
    GET  /panel/meta              -> reference panel version/scope info
"""

from __future__ import annotations

import json
import os
import threading
import time
import uuid

from flask import Flask, request, jsonify, Response
from flask_cors import CORS

import jobs
from pipeline.wrapper import run_full_pipeline, load_panel_meta, debug_panel_state
from pipeline.panel_builder import run_panel_build
from pipeline.genome_fetch import fetch_genomes
from pipeline.alphafold import find_alphafold_structure
from pipeline.pockets import detect_pockets
from pipeline.structural_similarity import find_similar_structures, load_manifest
from pipeline.docking import run_docking
from pipeline.pdb_validation import find_pdb_structures
from pipeline.chat_assistant import ask as chat_ask, list_models as chat_list_models

UPLOAD_DIR = os.path.join(os.path.dirname(__file__), "uploads")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")
GENOMES_DIR = os.path.join(os.path.dirname(__file__), "genomes")
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(GENOMES_DIR, exist_ok=True)

app = Flask(__name__)
CORS(app)  # relax during local dev; scope this down before deploying


@app.route("/upload", methods=["POST"])
def upload():
    if "file" not in request.files:
        return jsonify({"error": "No file provided"}), 400

    f = request.files["file"]
    if not f.filename:
        return jsonify({"error": "Empty filename"}), 400

    if not f.filename.lower().endswith((".fasta", ".fa", ".fna")):
        return jsonify({"error": "Expected a FASTA file (.fasta/.fa/.fna)"}), 400

    saved_name = f"{uuid.uuid4()}_{f.filename}"
    saved_path = os.path.join(UPLOAD_DIR, saved_name)
    f.save(saved_path)

    job_id = jobs.create_job(input_path=saved_path)
    return jsonify({"job_id": job_id})


@app.route("/jobs/<job_id>/start", methods=["POST"])
def start_job(job_id):
    job = jobs.get_job(job_id)
    if not job:
        return jsonify({"error": "Unknown job_id"}), 404

    if job["status"] in ("running", "done"):
        return jsonify({"status": job["status"]})

    body = request.json if request.is_json else {}
    organism = body.get("organism")
    detector = body.get("detector", "rgi")  # "rgi" (default, matches the panel) or "amrfinder"

    thread = threading.Thread(target=_run_job, args=(job_id, organism, detector), daemon=True)
    thread.start()

    jobs.update_job(job_id, status="running", stage="Starting")
    return jsonify({"status": "started", "detector": detector})


def _run_job(job_id: str, organism: str | None, detector: str):
    job = jobs.get_job(job_id)
    workdir = os.path.join(RESULTS_DIR, job_id)
    os.makedirs(workdir, exist_ok=True)

    def progress(stage: str):
        jobs.update_job(job_id, stage=stage)

    try:
        result = run_full_pipeline(job["input_path"], workdir, organism=organism,
                                    detector=detector, progress_cb=progress)
        jobs.update_job(job_id, status="done", stage="Complete", result=result)
    except Exception as e:
        jobs.update_job(job_id, status="error", stage="Failed", error=str(e))


@app.route("/jobs/<job_id>/stream")
def stream(job_id):
    def event_stream():
        last_stage = None
        while True:
            job = jobs.get_job(job_id)
            if not job:
                yield f"data: {json.dumps({'error': 'Unknown job_id'})}\n\n"
                return

            if job["stage"] != last_stage:
                yield f"data: {json.dumps({'status': job['status'], 'stage': job['stage']})}\n\n"
                last_stage = job["stage"]

            if job["status"] in ("done", "error"):
                return

            time.sleep(1)

    return Response(event_stream(), mimetype="text/event-stream")


@app.route("/jobs/<job_id>/results")
def results(job_id):
    job = jobs.get_job(job_id)
    if not job:
        return jsonify({"error": "Unknown job_id"}), 404
    if job["status"] == "error":
        return jsonify({"error": job.get("error") or "Job failed", "stage": job.get("stage")}), 500
    if job["status"] != "done":
        return jsonify({"status": job["status"], "stage": job["stage"]}), 202
    return jsonify(job["result"])


@app.route("/genomes/fetch", methods=["POST"])
def start_genome_fetch():
    """
    Downloads N genomes for a given taxon via NCBI's official 'datasets'
    CLI, saved into backend/genomes/<label>/ as flat FASTA files ready
    to hand to /panel/build's "directory" field.

    Body (JSON):
        taxon           (required) e.g. "Escherichia coli", "Klebsiella pneumoniae"
        limit           (required) how many genomes to fetch
        assembly_level  default "complete" ("complete"|"chromosome"|"scaffold"|"contig")
        label           default: slugified taxon — subfolder name under backend/genomes/
    """
    body = request.json if request.is_json else {}
    taxon = body.get("taxon")
    limit = body.get("limit")
    if not taxon or not limit:
        return jsonify({"error": "'taxon' and 'limit' are required"}), 400

    label = body.get("label") or taxon.replace(" ", "_")
    outdir = os.path.join(GENOMES_DIR, label)

    job_id = jobs.create_job(input_path=outdir)
    jobs.update_job(job_id, status="running", stage="Starting genome fetch")

    thread = threading.Thread(
        target=_run_genome_fetch_job,
        args=(job_id, taxon, int(limit), body.get("assembly_level", "complete"), outdir),
        daemon=True,
    )
    thread.start()

    return jsonify({"job_id": job_id, "status": "started"})


def _run_genome_fetch_job(job_id: str, taxon: str, limit: int, assembly_level: str, outdir: str):
    def progress(stage: str):
        jobs.update_job(job_id, stage=stage)

    try:
        result = fetch_genomes(taxon, limit, outdir, assembly_level=assembly_level, progress_cb=progress)
        jobs.update_job(job_id, status="done", stage="Complete", result=result)
    except Exception as e:
        jobs.update_job(job_id, status="error", stage="Failed", error=str(e))


@app.route("/panel/meta")
def panel_meta():
    return jsonify(load_panel_meta())


@app.route("/panel/debug")
def panel_debug():
    """
    Diagnostic route — shows exactly what this running backend process
    sees when loading the panel for a given detector (file path, whether
    it exists, how many genes loaded, whether a known gene like CRP is
    present). Optional ?detector=rgi|amrfinder query param — matches
    resolve_panel_dir()'s logic for picking a detector-specific panel.
    """
    detector = request.args.get("detector")
    try:
        return jsonify(debug_panel_state(detector=detector))
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/ontology/debug")
def ontology_debug():
    """
    Diagnostic route for the CARD ontology module — confirms the .obo
    file downloads and parses correctly, and does a real lookup on a
    known ARO term, independent of any specific analysis job. Useful for
    isolating "did the ontology feature itself break" from "did this
    particular gene's ARO lookup fail" when troubleshooting.

    Optional ?aro=<accession> to test a specific term (bare number or
    "ARO:" prefixed) — defaults to a well-known CTX-M beta-lactamase term.
    """
    # Placeholder test accession — not independently verified against
    # CARD's current ontology in this session. If this specific ID
    # doesn't resolve, that alone doesn't mean the module is broken —
    # pass a real accession you've confirmed via RGI's own output
    # (the "ARO" column) or via card.mcmaster.ca directly.
    test_aro = request.args.get("aro", "3001864")
    try:
        from pipeline.ontology import format_aro_id, get_term_info, get_lineage, get_ontology_siblings, OBO_PATH
        info = get_term_info(test_aro)
        return jsonify({
            "obo_cached_at": OBO_PATH,
            "obo_file_exists": os.path.exists(OBO_PATH),
            "test_aro_id": format_aro_id(test_aro),
            "test_term_info": info,
            "test_lineage": get_lineage(test_aro),
            "test_siblings_count": len(get_ontology_siblings(test_aro)),
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/alphafold/gene/<gene>")
def alphafold_gene(gene):
    """
    On-demand structural lookup for a single gene — called only when the
    user clicks "Analyze with AlphaFold" on a specific gene in the
    dashboard, never automatically for every detected gene (would slow
    down normal AMR analysis for no benefit until someone actually wants
    the structure).

    Optional ?organism=<name> query param narrows the UniProt search —
    recommended, since short AMR gene symbols (e.g. 'CRP', 'SHV-11') are
    reused across many organisms and an unfiltered search can return the
    wrong species' protein.
    """
    organism = request.args.get("organism")
    try:
        result = find_alphafold_structure(gene, organism=organism)
        return jsonify(result)
    except Exception as e:
        return jsonify({"gene": gene, "available": False, "error": str(e)}), 500


@app.route("/pockets/gene/<gene>")
def pockets_gene(gene):
    """
    On-demand binding-pocket detection via fpocket — only called after
    the user has already fetched a structure (needs its pdb_url passed
    in), and only when they click "Detect binding pockets" specifically.
    Requires fpocket installed on the host; see pockets.py's docstring.
    """
    pdb_url = request.args.get("pdb_url")
    if not pdb_url:
        return jsonify({"gene": gene, "available": False, "error": "pdb_url query param required"}), 400
    try:
        pockets = detect_pockets(pdb_url)
        return jsonify({"gene": gene, "available": True, "pockets": pockets})
    except Exception as e:
        return jsonify({"gene": gene, "available": False, "error": str(e)}), 500


@app.route("/similarity/gene/<gene>")
def similarity_gene(gene):
    """
    On-demand structural similarity search via Foldseek, comparing this
    gene's cached structure against every other structure previously
    fetched through this tool (backend/structure_cache/) — not the full
    AlphaFold DB / PDB. Requires ?uniprot=<accession> (the frontend
    already has this after a successful AlphaFold lookup). Requires
    Foldseek installed on the host; see structural_similarity.py's
    docstring.
    """
    uniprot = request.args.get("uniprot")
    if not uniprot:
        return jsonify({"gene": gene, "available": False, "error": "uniprot query param required"}), 400
    try:
        hits = find_similar_structures(uniprot)
        manifest = load_manifest()
        for h in hits:
            entry = manifest.get(h["target_accession"], {})
            h["gene_name"] = entry.get("gene_name")
            h["organism"] = entry.get("organism")
        return jsonify({"gene": gene, "available": True, "hits": hits})
    except Exception as e:
        return jsonify({"gene": gene, "available": False, "error": str(e)}), 500


@app.route("/docking/gene/<gene>")
def docking_gene(gene):
    """
    On-demand molecular docking — screens a small, keyword-matched set of
    real inhibitor compounds (fetched live from PubChem) against this
    gene's top detected binding pocket, using AutoDock Vina.

    This is the heaviest structural feature — requires fpocket already
    having been run (re-runs it here to get pocket geometry), plus vina
    and obabel installed on the host. See docking.py's docstring for the
    full pipeline and installation instructions.

    Query params:
        pdb_url         (required) from a prior /alphafold/gene lookup
        drug_class      (optional) used to pick a relevant inhibitor panel
        gene_family     (optional) same purpose
        custom_ligands  (optional) comma-separated compound names to test
                         instead of automatic selection, e.g.
                         "Acetyl-CoA,Coenzyme A,Malonyl-CoA" — each is
                         still fetched live from PubChem, so must be a
                         real, resolvable name
    """
    pdb_url = request.args.get("pdb_url")
    drug_class = request.args.get("drug_class")
    gene_family = request.args.get("gene_family")
    custom_ligands_raw = request.args.get("custom_ligands")
    custom_ligands = [n.strip() for n in custom_ligands_raw.split(",") if n.strip()] if custom_ligands_raw else None
    if not pdb_url:
        return jsonify({"gene": gene, "available": False, "error": "pdb_url query param required"}), 400

    try:
        pockets = detect_pockets(pdb_url, max_pockets=1)
        if not pockets:
            return jsonify({"gene": gene, "available": False, "error": "No binding pocket detected to dock against"})

        result = run_docking(pdb_url, pockets[0], drug_class, gene_family, gene_name=gene, custom_ligands=custom_ligands)
        return jsonify({"gene": gene, "available": True, **result})
    except Exception as e:
        return jsonify({"gene": gene, "available": False, "error": str(e)}), 500


@app.route("/validation/pdb/<gene>")
def pdb_validation_gene(gene):
    """
    Real experimental cross-reference — queries RCSB PDB's live database
    for actual crystal structures of this protein and their real bound
    ligands, so a Vina docking score can be checked against genuine
    experimental data rather than trusted in isolation.

    Query params:
        uniprot   (required) accession from a prior /alphafold/gene lookup
    """
    uniprot = request.args.get("uniprot")
    if not uniprot:
        return jsonify({"gene": gene, "available": False, "error": "uniprot query param required"}), 400
    try:
        structures = find_pdb_structures(uniprot)
        return jsonify({"gene": gene, "available": True, "structures": structures})
    except Exception as e:
        return jsonify({"gene": gene, "available": False, "error": str(e)}), 500


@app.route("/chat/models")
def chat_models():
    """
    Real, live list of OpenRouter's current model catalog — used by the
    frontend's model picker so it never shows a stale or guessed list.
    """
    try:
        return jsonify({"available": True, "models": chat_list_models()})
    except Exception as e:
        return jsonify({"available": False, "error": str(e)}), 500


@app.route("/chat", methods=["POST"])
def chat():
    """
    Chat assistant grounded in the current genome analysis. Expects a
    JSON body:
        {
          "question": "which genes are HGT?",
          "genes": [...],              # the current analysis's gene list
          "conversation_history": [    # prior turns this session, optional
            {"role": "user", "content": "..."},
            {"role": "assistant", "content": "..."}
          ],
          "active_gene_id": "CTX-M-15", # optional, whichever gene is selected
          "model": "anthropic/claude-3.5-sonnet",  # optional, from the model picker
          "structural_context": {      # optional — whatever's been fetched for the active gene
            "structure": {...},          # AlphaFold result, if fetched
            "pockets": {...},            # pocket detection result, if run
            "docking": {...},            # docking result, if run
            "pdb_validation": {...}      # real PDB cross-reference, if checked
          }
        }

    Requires OPENROUTER_API_KEY to be set before the backend starts.
    If "model" isn't given, falls back to the OPENROUTER_MODEL
    environment variable — see chat_assistant.py's docstring.
    """
    body = request.get_json(force=True) or {}
    question = body.get("question")
    genes = body.get("genes", [])
    conversation_history = body.get("conversation_history", [])
    active_gene_id = body.get("active_gene_id")
    model = body.get("model")
    structural_context = body.get("structural_context")

    if not question:
        return jsonify({"available": False, "error": "question is required"}), 400

    try:
        answer = chat_ask(question, genes, conversation_history, active_gene_id=active_gene_id,
                           model=model, structural_context=structural_context)
        return jsonify({"available": True, "answer": answer})
    except Exception as e:
        return jsonify({"available": False, "error": str(e)}), 500


@app.route("/panel/build", methods=["POST"])
def start_panel_build():
    """
    Kicks off a full panel rebuild from a folder of genomes already on
    disk (e.g. /Users/pauloa/Desktop/Virus/AMR/Ecoli1). This is a long
    job (30-60+ min for dozens of genomes) — reuses the same job store,
    so progress can be watched at /jobs/<job_id>/stream and the final
    summary read at /jobs/<job_id>/results, exactly like a single upload.

    Body (JSON):
        directory        (required) absolute path to a folder of FASTA files
        version           default "v1" — writes to backend/panel/<version>/
        detector          default "rgi"
        kmer              default 21
        permutations      default 999  (Mantel test permutations)
        clonal_r          default 0.4
        hgt_r             default 0.2
        phi_threshold     default 0.3
        scatter_threshold default 0.7
        taxon_scope       default "Enterobacteriaceae"
    """
    body = request.json if request.is_json else {}
    directory = body.get("directory")
    if not directory or not os.path.isdir(directory):
        return jsonify({"error": f"'directory' must be a valid path on this machine. Got: {directory!r}"}), 400

    job_id = jobs.create_job(input_path=directory)
    jobs.update_job(job_id, status="running", stage="Starting panel build")

    thread = threading.Thread(target=_run_panel_build_job, args=(job_id, body), daemon=True)
    thread.start()

    return jsonify({"job_id": job_id, "status": "started"})


def _run_panel_build_job(job_id: str, body: dict):
    def progress(stage: str):
        jobs.update_job(job_id, stage=stage)

    try:
        result = run_panel_build(
            genome_dir=body["directory"],
            version=body.get("version", "v1"),
            kmer=body.get("kmer", 21),
            detector=body.get("detector", "rgi"),
            permutations=body.get("permutations", 999),
            clonal_r=body.get("clonal_r", 0.4),
            hgt_r=body.get("hgt_r", 0.2),
            phi_threshold=body.get("phi_threshold", 0.3),
            scatter_threshold=body.get("scatter_threshold", 0.7),
            taxon_scope=body.get("taxon_scope", "Enterobacteriaceae"),
            progress_cb=progress,
        )
        jobs.update_job(job_id, status="done", stage="Complete", result=result)
    except Exception as e:
        jobs.update_job(job_id, status="error", stage="Failed", error=str(e))


if __name__ == "__main__":
    # Deliberately a different port from the existing KALI app.
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5050)), debug=os.environ.get("FLASK_DEBUG") == "1")
