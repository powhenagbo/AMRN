"""
pipeline/chat_assistant.py
=============================
A chat assistant embedded in the dashboard, grounded in the user's
current genome analysis — answers questions like "which genes are HGT?"
or "explain this Mantel r value" using the actual detected genes, not
just general knowledge.

Uses OpenRouter (openrouter.ai) — an OpenAI-compatible API gateway
giving access to many providers/models through one endpoint. Requires
an OpenRouter API key, set as an environment variable before starting
the backend:

    export OPENROUTER_API_KEY=sk-or-...
    export OPENROUTER_MODEL=anthropic/claude-3.5-sonnet   # or whichever model you use
    python app.py

Get a key at https://openrouter.ai/keys — model IDs are provider-
prefixed (e.g. "openai/gpt-4o", "anthropic/claude-3.5-sonnet",
"meta-llama/llama-3.1-70b-instruct"); check openrouter.ai/models for the
current catalog rather than assuming a specific string here, since
available models and exact IDs change over time.
"""

from __future__ import annotations

import os

import requests

OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"
OPENROUTER_MODELS_URL = "https://openrouter.ai/api/v1/models"
REQUEST_TIMEOUT = 30

MAX_GENES_IN_CONTEXT = 100  # keeps the prompt bounded even on large genomes


def _get_api_key():
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        raise RuntimeError(
            "No OpenRouter API key configured. Set it before starting the "
            "backend:\n"
            "  export OPENROUTER_API_KEY=sk-or-...\n"
            "  python app.py\n"
            "Get a key at https://openrouter.ai/keys"
        )
    return api_key


def list_models() -> list[dict]:
    """
    Fetches OpenRouter's real, current model catalog — used so the
    frontend can offer an actual live list to pick from, rather than a
    hardcoded guess that could be stale or wrong. This endpoint is public
    (no API key needed to list models, only to actually use one).

    Returns [{"id": "...", "name": "..."}], sorted alphabetically by id.
    """
    try:
        resp = requests.get(OPENROUTER_MODELS_URL, timeout=REQUEST_TIMEOUT)
    except requests.RequestException as e:
        raise RuntimeError(f"Could not reach OpenRouter's model list: {e}")

    if resp.status_code != 200:
        raise RuntimeError(f"OpenRouter model list returned HTTP {resp.status_code}: {resp.text[:300]}")

    data = resp.json()
    models = data.get("data", [])
    result = [{"id": m.get("id"), "name": m.get("name") or m.get("id")} for m in models if m.get("id")]
    result.sort(key=lambda m: m["id"])
    return result


def _build_context_summary(genes: list[dict], active_gene_id: str | None) -> str:
    """
    A condensed, token-efficient summary of the current analysis — not
    the full raw gene JSON (which includes mutation text, phylogenetic
    trees, etc.). Keeps each chat request fast and inexpensive.
    """
    lines = []
    for g in genes[:MAX_GENES_IN_CONTEXT]:
        lines.append(
            f"- {g.get('id')}: {g.get('klass')} (evidence: {g.get('evidence')}), "
            f"drug_class={g.get('drug')}, mantel_r={g.get('mantel_r')}, "
            f"island_zscore={g.get('score')}, gene_family={g.get('gene_family')}"
        )
    summary = "\n".join(lines) if lines else "(no genes detected in this analysis)"

    truncation_note = ""
    if len(genes) > MAX_GENES_IN_CONTEXT:
        truncation_note = f"\n(showing first {MAX_GENES_IN_CONTEXT} of {len(genes)} genes)"

    active_note = f"\n\nThe user is currently viewing gene: {active_gene_id}" if active_gene_id else ""

    return f"Current genome analysis \u2014 {len(genes)} genes detected:\n{summary}{truncation_note}{active_note}"


def _build_structural_context(structural_context: dict | None) -> str:
    """
    Summarizes whatever structural/docking results the user has already
    generated for the currently-viewed gene (AlphaFold structure, binding
    pockets, docking scores, real PDB cross-reference) — this is what
    lets the assistant actually discuss "this structure" or "this docking
    result" rather than only knowing the bare gene list. Only includes
    what's actually been run; doesn't trigger new analyses itself.
    """
    if not structural_context:
        return ""

    parts = []

    struct = structural_context.get("structure")
    if struct and struct.get("available"):
        parts.append(
            f"AlphaFold structure: protein={struct.get('protein')}, "
            f"uniprot={struct.get('uniprot')}, mean_pLDDT={struct.get('mean_plddt')}"
        )

    pockets = structural_context.get("pockets")
    if pockets and pockets.get("available"):
        pocket_list = pockets.get("pockets", [])
        if pocket_list:
            top = sorted(pocket_list, key=lambda p: p.get("druggability_score") or 0, reverse=True)[:3]
            pocket_summary = "; ".join(
                f"pocket #{p.get('pocket_id')} druggability={p.get('druggability_score')}" for p in top
            )
            parts.append(f"Binding pockets detected: {pocket_summary}")
        else:
            parts.append("Binding pockets: none detected")

    docking = structural_context.get("docking")
    if docking and docking.get("available"):
        results = [r for r in docking.get("results", []) if r.get("available")]
        if results:
            docking_summary = "; ".join(
                f"{r.get('ligand')}: {r.get('affinity_kcal_mol')} kcal/mol "
                f"(mean {r.get('mean_affinity_kcal_mol')}\u00b1{r.get('std_affinity_kcal_mol')}, "
                f"{r.get('n_replicates')}/{r.get('n_replicates_attempted')} runs)"
                for r in results
            )
            parts.append(f"Docking results: {docking_summary}")

    pdb = structural_context.get("pdb_validation")
    if pdb and pdb.get("available"):
        structures = pdb.get("structures", [])
        if structures:
            pdb_summary = "; ".join(
                f"{s.get('pdb_id')} (ligands: {', '.join(s.get('ligands', [])) or 'none'})" for s in structures
            )
            parts.append(f"Real experimental structures (RCSB PDB): {pdb_summary}")
        else:
            parts.append("Real experimental structures (RCSB PDB): none found")

    if not parts:
        return ""

    return "\n\nStructural analysis run so far for the currently-viewed gene:\n" + "\n".join(f"- {p}" for p in parts)


def ask(question: str, genes: list[dict], conversation_history: list[dict],
        active_gene_id: str | None = None, model: str | None = None,
        structural_context: dict | None = None) -> str:
    """
    Sends a question, grounded in the current analysis, to OpenRouter.

    model: a specific OpenRouter model ID to use for this request (e.g.
    from the frontend's model picker) — falls back to the OPENROUTER_MODEL
    environment variable if not given, and raises a clear error if
    neither is set.

    conversation_history: list of {"role": "user"|"assistant", "content": str}
    from prior turns in this chat session — passed through so the
    assistant has conversational memory within a session (not persisted
    across page reloads).

    structural_context: whatever structural/docking results the frontend
    has already fetched for the active gene (structure, pockets, docking,
    pdb_validation) — lets the assistant discuss those results directly.
    """
    api_key = _get_api_key()
    model = model or os.environ.get("OPENROUTER_MODEL")
    if not model:
        raise RuntimeError(
            "No model specified. Either pick one in the chat widget, or set "
            "a default before starting the backend:\n"
            "  export OPENROUTER_MODEL=anthropic/claude-3.5-sonnet\n"
            "  python app.py"
        )

    context = _build_context_summary(genes, active_gene_id)
    struct_context = _build_structural_context(structural_context)

    system_prompt = (
        "You are KALI AI, an assistant embedded in KALI-AMR-HGT, a tool that "
        "detects antimicrobial resistance genes and classifies each as "
        "clonally inherited or horizontally transferred using alignment-free "
        "compositional analysis and population-level phylogenetic "
        "correlation. Answer the user's questions about their current "
        "genome analysis using the data provided below. Be concise and "
        "accurate. If something isn't in the provided data, say so rather "
        "than guessing or inventing values.\n\n" + context + struct_context
    )

    messages = [{"role": "system", "content": system_prompt}] + conversation_history + [
        {"role": "user", "content": question}
    ]

    try:
        resp = requests.post(
            OPENROUTER_URL,
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            json={"model": model, "messages": messages, "max_tokens": 1000},
            timeout=REQUEST_TIMEOUT,
        )
    except requests.RequestException as e:
        raise RuntimeError(f"OpenRouter request failed: {e}")

    if resp.status_code != 200:
        raise RuntimeError(f"OpenRouter returned HTTP {resp.status_code}: {resp.text[:400]}")

    data = resp.json()
    try:
        return data["choices"][0]["message"]["content"]
    except (KeyError, IndexError):
        raise RuntimeError(f"Unexpected OpenRouter response shape: {str(data)[:400]}")
