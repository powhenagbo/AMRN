# Getting RGI + AMRFinderPlus working on Render

## What changed and why

Locally, `run_rgi()` shells out to `docker run finlaymaguire/rgi:latest` —
that works because Docker is running on your Mac and you've already pulled
the image. Render's containers don't expose a Docker daemon to run a nested
container, so that call always fails in production, regardless of the
Render runtime you pick.

The fix has two parts:
1. **`pipeline/wrapper.py`** — `run_rgi()` now checks an `RGI_MODE` env var.
   - Unset (default) → behaves exactly as before: shells out to
     `docker run finlaymaguire/rgi:latest`. **Your local setup is untouched.**
   - `RGI_MODE=native` → calls the `rgi` CLI directly, no Docker involved.
     This is the path the deployed container uses.
2. **`Dockerfile`** — a Docker-based image for the Render service itself
   (this is Render running *your app* in a container, not your app running
   *nested* Docker containers). It installs `rgi` and `amrfinder` directly
   via conda/bioconda, loads CARD reference data at build time, and sets
   `RGI_MODE=native` so `run_rgi()` uses the direct-call path.

`run_amrfinder()` needed no code changes — it already just calls the bare
`amrfinder` binary, which the Dockerfile installs.

## Files to drop in

| File | Destination |
|---|---|
| `wrapper.py` | replaces `amr-hgt-tool/backend/pipeline/wrapper.py` |
| `Dockerfile` | new file at `amr-hgt-tool/backend/Dockerfile` |
| `render.yaml` | replaces (or adds, if you don't have one) the root `render.yaml` |
| `backend.gitignore` | rename to `.gitignore`, place at `amr-hgt-tool/backend/.gitignore` |

## Steps

1. **Copy the files in** as above.

2. **Test the Docker build locally before pushing** — this is the step
   most likely to need iteration (CARD data URLs and package versions
   shift over time):
   ```bash
   cd amr-hgt-tool/backend
   docker build -t amrn-api .
   docker run -p 10000:10000 -e OPENROUTER_API_KEY=your-key amrn-api
   ```
   Then from another terminal, upload a real genome and start a job
   against `http://localhost:10000` (same requests your frontend sends)
   and confirm it completes with `detector: "rgi"` and separately with
   `detector: "amrfinder"`.

3. **Trim the tool list if you don't need it all.** The Dockerfile installs
   `fpocket`, `foldseek`, `autodock-vina`, `openbabel`, and
   `ncbi-datasets-cli` too, since the pipeline's structural-analysis
   modules (`docking.py`, `pockets.py`, `structural_similarity.py`,
   `genome_fetch.py`) need them. If you're not deploying those endpoints
   yet, comment those lines out of the Dockerfile — smaller image, faster
   builds.

4. **Push and connect the repo to Render as a Blueprint** (New → Blueprint,
   point it at this repo). Render will read the new `render.yaml` and
   build from the Dockerfile instead of the old native-Python buildpack.

5. **Set `OPENROUTER_API_KEY`** in the Render dashboard's environment
   variables for the service — don't put the real value in `render.yaml`
   itself (it's marked `sync: false` there on purpose, so Render prompts
   you for it instead of it living in git).

6. **Watch the first build closely.** The CARD data download
   (`https://card.mcmaster.ca/latest/data`) and `rgi load` step are the
   most likely things to break — CARD occasionally changes its release
   URL structure. If the build fails there, check
   https://card.mcmaster.ca/download for the current data URL and update
   the `wget` line in the Dockerfile.

## What stays exactly the same

- Local development: unset `RGI_MODE`, keep using Docker + your existing
  `finlaymaguire/rgi:latest` image and local `amrfinder` install. Nothing
  about your local workflow changes.
- `run_amrfinder()`, the Flask routes, the frontend — untouched.
