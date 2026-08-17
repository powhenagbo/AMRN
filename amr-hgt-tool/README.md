# AMR-HGT Tool — Standalone Scaffold

A separate web app from the existing KALI Flask app. Own codebase, own job
store, own port (5050 by default). Nothing here imports from or writes to
the existing KALI app's files.

## What's here

```
amr-hgt-tool/
├── backend/
│   ├── app.py                 # Flask app: upload, job start, SSE stream, results
│   ├── jobs.py                 # File-backed job store (jobs_store.json)
│   ├── pipeline/
│   │   └── wrapper.py           # Calls AMRFinderPlus + kali_hgt.py + panel lookup
│   ├── panel/v1/
│   │   └── meta.json             # Reference panel metadata (placeholder)
│   ├── uploads/                 # Uploaded FASTA files land here
│   ├── results/                 # Per-job working directories
│   └── requirements.txt
└── frontend/
    └── amr-hgt-dashboard.jsx    # The dashboard UI (currently uses hardcoded demo data)
```

## What still needs to happen before this is real (in order)

1. **Copy `kali_hgt.py` into `backend/pipeline/`**
   `pipeline/wrapper.py` imports `read_fasta` and `detect_islands` from it
   directly — it's a local import, so the file needs to physically sit next
   to `wrapper.py`.

2. **Install one or both AMR detectors** on whatever machine runs this
   backend — the caller picks per upload (`detector: "rgi"` or `"amrfinder"`
   in the `POST /jobs/{id}/start` body; defaults to `"rgi"`):
   - **CARD-RGI** (recommended default): requires Docker running, with
     `docker pull finlaymaguire/rgi:latest` done once. `wrapper.py`'s
     `run_rgi()` shells out to it exactly as in the AMR Pipeline Runbook.
     Use this when panel lookup accuracy matters — the reference panel's
     Mantel labels are keyed on RGI/CARD gene names.
   - **AMRFinderPlus**: requires the `amrfinder` binary + database
     (`amrfinder_update` run once). Faster, but gene names may not match
     the RGI-based panel, so more genes will fall back to island-only
     evidence (`evidence_source: "island_only"`) rather than a panel-derived
     Mantel label. Every result includes a `detector` block noting which
     one ran and why that matters.

3. **Fetch genomes to build the panel from** — optional if you already
   have a genome folder (e.g. from `download_o104h4_genomes.py` or
   elsewhere). To have the app download them instead, install NCBI's
   official `datasets` CLI:
   ```bash
   brew install ncbi-datasets-cli
   # or: conda install -c conda-forge ncbi-datasets-cli
   ```
   Then fetch N genomes for a taxon:
   ```bash
   curl -X POST http://127.0.0.1:5050/genomes/fetch \
     -H "Content-Type: application/json" \
     -d '{"taxon": "Klebsiella pneumoniae", "limit": 50, "assembly_level": "complete"}'
   # -> {"job_id": "...", "status": "started"}
   ```
   Watch it and check results the same way as any other job:
   ```bash
   curl http://127.0.0.1:5050/jobs/<job_id>/stream
   curl http://127.0.0.1:5050/jobs/<job_id>/results
   # -> {"downloaded": 50, "genome_dir": "backend/genomes/Klebsiella_pneumoniae", ...}
   ```
   The `genome_dir` it returns is exactly what you feed into `/panel/build`'s
   `directory` field next. This uses NCBI's official CLI rather than custom
   scraping — `download_o104h4_genomes.py` remains the right tool if you
   specifically want that curated 32-genome O104:H4 outbreak set instead.

4. **Build the reference panel** — now triggerable from the app itself,
   no manual CLI runbook needed. `amr_transmission.py` and
   `amr_cooccurrence.py` (your real modules) are already copied into
   `backend/pipeline/`, wired up by `backend/pipeline/panel_builder.py`.

   With the backend running, trigger a build via `curl` (or Postman) —
   this is a long job (30-60+ min for dozens of genomes), so it reuses
   the same job/SSE machinery as a single upload:

   ```bash
   curl -X POST http://127.0.0.1:5050/panel/build \
     -H "Content-Type: application/json" \
     -d '{
       "directory": "/Users/pauloa/Desktop/Virus/AMR/Ecoli1",
       "version": "v1",
       "detector": "rgi",
       "kmer": 21
     }'
   # -> {"job_id": "...", "status": "started"}
   ```

   Watch progress:
   ```bash
   curl http://127.0.0.1:5050/jobs/<job_id>/stream
   ```

   Check the final summary once done:
   ```bash
   curl http://127.0.0.1:5050/jobs/<job_id>/results
   # -> {"panel_built": true, "version": "v1", "n_genomes": 47,
   #     "n_genes_classified": 95, "n_cooccurrence_pairs": 1481, ...}
   ```

   This writes `transmission_classification.csv`, `cooccurrence_classification.csv`,
   and an updated `meta.json` straight into `backend/panel/v1/` — the exact
   files `wrapper.py`'s `apply_panel_labels()` reads for every single-genome
   upload afterward. No separate manual runbook execution needed anymore,
   though the original CLI scripts still work fine if you prefer them.

   **Note:** `directory` must be a path that exists on the same machine
   running the Flask backend — this doesn't accept browser uploads for the
   panel build, by design (see earlier discussion: pointing at a local
   folder is far faster than re-uploading dozens of large FASTA files).

   Without a built panel, `wrapper.py` still runs fine for single uploads —
   every gene just falls back to island-overlap-only evidence
   (`evidence_source: "island_only"`) instead of a panel-derived Mantel label.

5. **Wire the frontend to the backend.** In `amr-hgt-dashboard.jsx`, replace
   the hardcoded `const GENES = [...]` with a `fetch` against
   `/jobs/{job_id}/results`, and wire the "Upload FASTA" button to
   `POST /upload` → `POST /jobs/{job_id}/start` → poll or subscribe to
   `/jobs/{job_id}/stream`.

## Running locally

```bash
cd backend
pip install -r requirements.txt
python app.py
# serves on http://localhost:5050 — separate from the existing KALI app's port
```

## Notes

- `evidence_source` in each gene result tells you whether a classification
  came from the precomputed panel (`"panel"`, high confidence) or from
  island-overlap alone (`"island_only"`, lower confidence — this is the
  case for genes not yet seen in the reference panel, or genomes outside
  its taxonomic scope). Surface this distinction in the UI rather than
  showing a bare CLONAL/HGT badge with no confidence signal.
- This scaffold does not yet compute distance-to-panel for taxonomic-scope
  checking (Case 2 from earlier discussion — genome outside Enterobacteriaceae).
  That's the next layer to add once the panel itself exists.
