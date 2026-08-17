"""
jobs.py
=======
Minimal file-backed job store for the AMR-HGT tool.

Deliberately separate from KALI's existing jobs.json — this tool has its
own store, own file, own schema, so nothing about the existing KALI app
is touched or shared.

Job record shape:
{
  "id": "uuid",
  "status": "queued" | "running" | "done" | "error",
  "stage": "human-readable current step",
  "created_at": iso timestamp,
  "updated_at": iso timestamp,
  "input_path": path to uploaded FASTA,
  "result": { ... } | null,
  "error": str | null
}
"""

import json
import os
import threading
import uuid
from datetime import datetime, timezone

JOBS_FILE = os.path.join(os.path.dirname(__file__), "jobs_store.json")
_lock = threading.Lock()


def _now():
    return datetime.now(timezone.utc).isoformat()


def _read_all() -> dict:
    if not os.path.exists(JOBS_FILE):
        return {}
    with open(JOBS_FILE, "r") as f:
        try:
            return json.load(f)
        except json.JSONDecodeError:
            return {}


def _write_all(data: dict):
    tmp = JOBS_FILE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(tmp, JOBS_FILE)


def create_job(input_path: str) -> str:
    job_id = str(uuid.uuid4())
    with _lock:
        data = _read_all()
        data[job_id] = {
            "id": job_id,
            "status": "queued",
            "stage": "Queued",
            "created_at": _now(),
            "updated_at": _now(),
            "input_path": input_path,
            "result": None,
            "error": None,
        }
        _write_all(data)
    return job_id


def update_job(job_id: str, **fields):
    with _lock:
        data = _read_all()
        if job_id not in data:
            return
        data[job_id].update(fields)
        data[job_id]["updated_at"] = _now()
        _write_all(data)


def get_job(job_id: str) -> dict | None:
    with _lock:
        data = _read_all()
        return data.get(job_id)
