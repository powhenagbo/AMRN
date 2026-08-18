# Place this at amr-hgt-tool/backend/Dockerfile
#
# Installs the external CLI tools the pipeline needs so they exist inside
# the deployed container (previously only true on your local dev machine
# via Docker + `docker pull finlaymaguire/rgi:latest` and a locally
# installed `amrfinder`).
#
# RGI is installed directly via conda/bioconda here — NOT run via nested
# `docker run`, because Render's containers don't expose a Docker daemon
# to run a second, nested container. RGI_MODE=native (set below) tells
# the patched pipeline/wrapper.py::run_rgi() to call the `rgi` binary
# directly instead of shelling out to `docker run`. Local dev is
# untouched: without RGI_MODE set, run_rgi() still defaults to the
# Docker-based path you're already using.

FROM condaforge/miniforge3:latest

WORKDIR /app

# System deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Bioinformatics tools via bioconda/conda-forge.
# Trim this list to whatever your deployed feature set actually needs —
# fpocket/foldseek/vina/openbabel/datasets are only required if you're
# also deploying the structural-analysis endpoints (pockets, docking,
# structural similarity, genome_fetch). Comment them out if you're only
# deploying the core AMR/HGT detection pipeline, to keep the image and
# build time smaller.
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.11 \
    rgi \
    ncbi-amrfinderplus \
    fpocket \
    foldseek \
    autodock-vina \
    openbabel \
    ncbi-datasets-cli \
    && mamba clean -afy

# AMRFinderPlus needs its database initialized once at build time
RUN amrfinder_update || true

# RGI needs CARD reference data loaded before `rgi main` will run.
# Downloads and loads the current CARD database at build time so it's
# baked into the image (avoids doing this on every cold start).
RUN wget -q https://card.mcmaster.ca/latest/data -O /tmp/card-data.tar.bz2 \
    && mkdir -p /tmp/card-data \
    && tar -xjf /tmp/card-data.tar.bz2 -C /tmp/card-data \
    && rgi load --card_json /tmp/card-data/card.json --local \
    && rm -rf /tmp/card-data /tmp/card-data.tar.bz2

# Tell run_rgi() to call the rgi binary directly (no nested Docker)
ENV RGI_MODE=native

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV PORT=10000
EXPOSE 10000

# --workers 1 --threads 8: a single sync worker holds the whole process
# hostage for the entire duration of any /jobs/<id>/stream SSE
# connection (it blocks in a time.sleep loop for as long as a job runs,
# often several minutes on a real genome) — during that time NOTHING
# else can be served, not other users' requests, not even Render's
# health checks, and gunicorn's own --timeout eventually kills the
# worker mid-stream (seen as "WORKER TIMEOUT" / SIGKILL in the logs).
# Threads let the same process interleave multiple requests — Python
# releases the GIL during time.sleep()/I/O, so one job's progress
# stream no longer blocks everyone else. --timeout raised to give a
# real end-to-end pipeline run room to finish without getting killed.
CMD gunicorn --bind 0.0.0.0:$PORT --workers 1 --threads 8 --timeout 1800 app:app
