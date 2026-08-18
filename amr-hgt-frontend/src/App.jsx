import React, { useState, useRef, useEffect } from "react";
import {
  Dna, Search, ChevronRight, Info, Upload, Loader2, AlertCircle,
  Download, Layers, ChevronDown, ChevronUp, CheckCircle2, Atom,
} from "lucide-react";

// Point this at wherever the Flask backend is actually running.
// Reads VITE_API_BASE_URL at build time (set in Vercel's Environment
// Variables for production), falling back to localhost for local dev
// when that var isn't set.
const API_BASE = import.meta.env.VITE_API_BASE_URL || "http://127.0.0.1:5050";

// ── Shared job-runner hook ──────────────────────────────────────────────
// Every long-running backend action (single-genome analysis, genome fetch,
// panel build) follows the same shape: POST to start it, get a job_id,
// stream progress over SSE, then read the final result. This hook covers
// that pattern once instead of three times.
function useJobRunner() {
  const [status, setStatus] = useState("idle"); // idle | starting | running | done | error
  const [stage, setStage] = useState("");
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const esRef = useRef(null);

  function reset() {
    setStatus("idle");
    setStage("");
    setError(null);
    setResult(null);
  }

  async function fetchResult(jobId) {
    try {
      const res = await fetch(`${API_BASE}/jobs/${jobId}/results`);
      const data = await res.json();

      if (res.status === 202) {
        setTimeout(() => fetchResult(jobId), 1000);
        return;
      }
      if (!res.ok) {
        setError(data.error || "Job failed");
        setStatus("error");
        return;
      }
      setResult(data);
      setStatus("done");
      setStage("Complete");
    } catch (err) {
      setError(err.message);
      setStatus("error");
    }
  }

  // starter: async () => jobId  — caller does whatever POST(s) it needs
  // and returns the resulting job_id; this hook takes it from there.
  async function run(starter) {
    reset();
    setStatus("starting");
    try {
      const jobId = await starter();
      setStatus("running");

      const es = new EventSource(`${API_BASE}/jobs/${jobId}/stream`);
      esRef.current = es;

      es.onmessage = (event) => {
        const data = JSON.parse(event.data);
        if (data.error) {
          setError(data.error);
          setStatus("error");
          es.close();
          return;
        }
        setStage(data.stage);
        if (data.status === "done" || data.status === "error") {
          es.close();
          fetchResult(jobId);
        }
      };
      es.onerror = () => es.close();
    } catch (err) {
      setError(err.message);
      setStatus("error");
    }
  }

  useEffect(() => {
    return () => {
      if (esRef.current) esRef.current.close();
    };
  }, []);

  return { status, stage, error, result, run, reset };
}

function Badge({ klass }) {
  const isHgt = klass === "HGT";
  const isAmbiguous = klass !== "HGT" && klass !== "Clonal" && klass !== "CLONAL";
  const color = isHgt ? "#5FA777" : isAmbiguous ? "#8592A6" : "#A8825E";
  return (
    <span
      style={{
        display: "inline-flex",
        alignItems: "center",
        gap: 6,
        padding: "3px 10px",
        borderRadius: 999,
        fontSize: 12,
        fontFamily: "'IBM Plex Mono', monospace",
        fontWeight: 500,
        letterSpacing: "0.02em",
        background: `${color}1F`,
        color,
        border: `1px solid ${color}59`,
      }}
    >
      <span style={{ width: 6, height: 6, borderRadius: "50%", background: color }} />
      {klass}
    </span>
  );
}

function EvidenceTag({ evidence }) {
  if (!evidence) return null;
  const config = {
    panel: { label: "panel-derived", color: "#3FA796" },
    panel_normalized: { label: "panel-derived (name-matched)", color: "#3FA796" },
    island_only: { label: "island-only", color: "#8592A6" },
  }[evidence] || { label: evidence, color: "#8592A6" };

  return (
    <span
      style={{
        fontSize: 10,
        fontFamily: "'IBM Plex Mono', monospace",
        color: config.color,
        textTransform: "uppercase",
        letterSpacing: "0.04em",
      }}
    >
      {config.label}
    </span>
  );
}

function Row({ label, value, mono }) {
  return (
    <div style={{ display: "flex", justifyContent: "space-between" }}>
      <span style={{ color: "#6B7A93" }}>{label}</span>
      <span style={{ fontFamily: mono ? "'IBM Plex Mono', monospace" : "inherit", color: "#E7ECF3" }}>
        {value ?? "—"}
      </span>
    </div>
  );
}

function ProgressLine({ status, stage, error }) {
  if (status === "starting" || status === "running") {
    return (
      <div style={{ marginTop: 14, fontSize: 13, color: "#8592A6", fontFamily: "'IBM Plex Mono', monospace", display: "flex", alignItems: "center", gap: 8 }}>
        <Loader2 size={14} className="spin" color="#E8A33D" />
        {stage || "Working..."}
      </div>
    );
  }
  if (status === "error") {
    return (
      <div style={{ marginTop: 14, fontSize: 13, color: "#E85D3D", display: "flex", alignItems: "flex-start", gap: 8, whiteSpace: "pre-wrap" }}>
        <AlertCircle size={14} style={{ flexShrink: 0, marginTop: 2 }} />
        {error || "Something went wrong."}
      </div>
    );
  }
  return null;
}

// ── Embedded 3D structure viewer ────────────────────────────────────────
// Uses 3Dmol.js, loaded via a <script> tag in index.html (window.$3Dmol) —
// simpler and more bundler-friendly than the npm package for a viewer
// like this. Colors by pLDDT confidence using AlphaFold's own convention
// (blue = very high, cyan = high, yellow = low, orange = very low).
function plddtColor(bfactor) {
  if (bfactor > 90) return "#1E6FEB";
  if (bfactor > 70) return "#65CBF3";
  if (bfactor > 50) return "#FFDB13";
  return "#FF7D45";
}

// Plain-language interpretation of classification confidence — NOT a
// blended single score across Mantel r / pLDDT / docking (those measure
// fundamentally different things; averaging them would misrepresent
// what's actually known). This interprets the Mantel-test statistics
// specifically, honestly reflecting sample size and significance rather
// than treating any nonzero r as equally trustworthy.
function classificationConfidence(evidence, pValue, nCarriers) {
  if (evidence !== "panel" && evidence !== "panel_normalized") {
    return { label: "No statistical test", color: "#6B7A93", detail: "Classification falls back to compositional island overlap only — no Mantel test was run for this gene." };
  }
  if (pValue == null) {
    return { label: "Panel match, no p-value recorded", color: "#6B7A93", detail: "" };
  }
  if (pValue <= 0.01 && nCarriers >= 20) {
    return { label: "Strong", color: "#3FA796", detail: `p=${pValue.toFixed(4)}, n=${nCarriers} carriers — high statistical power.` };
  }
  if (pValue <= 0.05) {
    return { label: "Moderate", color: "#E8A33D", detail: `p=${pValue.toFixed(4)}, n=${nCarriers ?? "?"} carriers.` };
  }
  return { label: "Weak / not significant", color: "#E85D3D", detail: `p=${pValue.toFixed(4)} — treat this classification cautiously.` };
}

function StructureViewer({ pdbUrl }) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);
  const [status, setStatus] = useState("loading"); // loading | ready | error

  useEffect(() => {
    if (!pdbUrl || !containerRef.current) return;
    if (!window.$3Dmol) {
      setStatus("error");
      return;
    }

    setStatus("loading");
    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: "#0B1220",
    });
    viewerRef.current = viewer;

    fetch(pdbUrl)
      .then((r) => {
        if (!r.ok) throw new Error("Failed to fetch structure file");
        return r.text();
      })
      .then((pdbText) => {
        viewer.addModel(pdbText, "pdb");
        viewer.setStyle({}, { cartoon: { colorfunc: (atom) => plddtColor(atom.b) } });
        viewer.zoomTo();
        viewer.render();
        setStatus("ready");
      })
      .catch(() => setStatus("error"));

    return () => {
      try {
        viewer.clear();
      } catch (e) {
        // viewer may already be torn down — safe to ignore
      }
    };
  }, [pdbUrl]);

  return (
    <div style={{ position: "relative", marginTop: 12 }}>
      <div
        ref={containerRef}
        style={{
          width: "100%",
          height: 280,
          borderRadius: 8,
          overflow: "hidden",
          background: "#0B1220",
          border: "1px solid #1B2740",
        }}
      />
      {status === "loading" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#8592A6", fontFamily: "'IBM Plex Mono', monospace" }}>
          <Loader2 size={14} className="spin" color="#E8A33D" />
          Loading structure...
        </div>
      )}
      {status === "error" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#E85D3D" }}>
          <AlertCircle size={14} />
          Couldn't load the 3D viewer
        </div>
      )}
      {status === "ready" && (
        <div style={{ display: "flex", gap: 12, marginTop: 8, fontSize: 10, fontFamily: "'IBM Plex Mono', monospace", color: "#6B7A93" }}>
          <span><span style={{ color: "#1E6FEB" }}>■</span> Very high</span>
          <span><span style={{ color: "#65CBF3" }}>■</span> High</span>
          <span><span style={{ color: "#FFDB13" }}>■</span> Low</span>
          <span><span style={{ color: "#FF7D45" }}>■</span> Very low</span>
        </div>
      )}
    </div>
  );
}

// ── Docked pose viewer ────────────────────────────────────────────────
// Shows the receptor (dimmed cartoon, pLDDT-colored) with the docked
// ligand overlaid as bright magenta sticks at Vina's actual predicted
// position/orientation — not just a bare affinity number.
function DockingViewer({ pdbUrl, ligandPdbqt, ligandName }) {
  const containerRef = useRef(null);
  const [status, setStatus] = useState("loading");

  useEffect(() => {
    if (!pdbUrl || !ligandPdbqt || !containerRef.current) return;
    if (!window.$3Dmol) {
      setStatus("error");
      return;
    }

    setStatus("loading");
    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: "#0B1220",
    });

    fetch(pdbUrl)
      .then((r) => {
        if (!r.ok) throw new Error("Failed to fetch receptor structure");
        return r.text();
      })
      .then((receptorText) => {
        viewer.addModel(receptorText, "pdb");
        viewer.setStyle({ model: 0 }, { cartoon: { colorfunc: (atom) => plddtColor(atom.b), opacity: 0.85 } });

        const ligandModel = viewer.addModel(ligandPdbqt, "pdbqt");
        viewer.setStyle({ model: 1 }, { stick: { color: "#FF3EC9", radius: 0.22 } });

        // Defensive check — a silently-empty ligand model (e.g. from a
        // parsing mismatch) would otherwise just render nothing, with no
        // indication anything went wrong. Verify it actually has atoms
        // before declaring success.
        const ligandAtomCount = ligandModel.selectedAtoms({}).length;
        if (ligandAtomCount === 0) {
          setStatus("error");
          return;
        }

        viewer.zoomTo({ model: 1 });
        viewer.zoom(0.4); // pull back slightly from the ligand alone to show pocket context
        viewer.render();
        setStatus("ready");
      })
      .catch(() => setStatus("error"));

    return () => {
      try {
        viewer.clear();
      } catch (e) {
        // viewer may already be torn down — safe to ignore
      }
    };
  }, [pdbUrl, ligandPdbqt]);

  return (
    <div style={{ position: "relative", marginTop: 12 }}>
      <div
        ref={containerRef}
        style={{
          width: "100%",
          height: 280,
          borderRadius: 8,
          overflow: "hidden",
          background: "#0B1220",
          border: "1px solid #1B2740",
        }}
      />
      {status === "loading" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#8592A6", fontFamily: "'IBM Plex Mono', monospace" }}>
          <Loader2 size={14} className="spin" color="#E8A33D" />
          Loading docked pose...
        </div>
      )}
      {status === "error" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#E85D3D" }}>
          <AlertCircle size={14} />
          Couldn't load the docked pose
        </div>
      )}
      {status === "ready" && (
        <div style={{ fontSize: 10, color: "#6B7A93", marginTop: 8, fontFamily: "'IBM Plex Mono', monospace" }}>
          <span style={{ color: "#FF3EC9" }}>■</span> {ligandName} (predicted pose) — protein colored by pLDDT confidence
        </div>
      )}
    </div>
  );
}

// ── Real experimental structure viewer ──────────────────────────────────
// Renders an actual RCSB crystal structure — NOT colored by pLDDT (real
// structures store crystallographic B-factors in that column, a
// completely different measurement; reusing the confidence color scheme
// here would be misleading). Uses a standard spectrum/rainbow coloring
// by chain position instead, the conventional way to color a structure
// with no confidence data.
function RealStructureViewer({ pdbId }) {
  const containerRef = useRef(null);
  const [status, setStatus] = useState("loading");

  useEffect(() => {
    if (!pdbId || !containerRef.current) return;
    if (!window.$3Dmol) {
      setStatus("error");
      return;
    }

    setStatus("loading");
    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: "#0B1220",
    });

    fetch(`https://files.rcsb.org/download/${pdbId}.pdb`)
      .then((r) => {
        if (!r.ok) throw new Error("Failed to fetch structure from RCSB");
        return r.text();
      })
      .then((pdbText) => {
        viewer.addModel(pdbText, "pdb");
        viewer.setStyle({}, { cartoon: { color: "spectrum" } });
        viewer.zoomTo();
        viewer.render();
        setStatus("ready");
      })
      .catch(() => setStatus("error"));

    return () => {
      try {
        viewer.clear();
      } catch (e) {
        // viewer may already be torn down — safe to ignore
      }
    };
  }, [pdbId]);

  return (
    <div style={{ position: "relative", marginTop: 10 }}>
      <div
        ref={containerRef}
        style={{
          width: "100%",
          height: 260,
          borderRadius: 8,
          overflow: "hidden",
          background: "#0B1220",
          border: "1px solid #1B2740",
        }}
      />
      {status === "loading" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#8592A6", fontFamily: "'IBM Plex Mono', monospace" }}>
          <Loader2 size={14} className="spin" color="#E8A33D" />
          Loading {pdbId} from RCSB...
        </div>
      )}
      {status === "error" && (
        <div style={{ position: "absolute", inset: 0, display: "flex", alignItems: "center", justifyContent: "center", gap: 8, fontSize: 12, color: "#E85D3D" }}>
          <AlertCircle size={14} />
          Couldn't load {pdbId}
        </div>
      )}
      {status === "ready" && (
        <div style={{ fontSize: 10, color: "#6B7A93", marginTop: 6, fontFamily: "'IBM Plex Mono', monospace" }}>
          {pdbId} — real experimental structure, colored by chain position (no confidence data applies here)
        </div>
      )}
    </div>
  );
}

// ── Phylogenetic tree renderer ───────────────────────────────────────────
// Renders the real neighbor-joining tree (see phylo_placement.py) as a
// simple SVG dendrogram — no charting library, consistent with the rest
// of this app's minimal-dependency approach. Layout: leaves get evenly
// spaced y-positions (in traversal order), x-position is cumulative
// branch length from the root, internal nodes connect their children
// with a vertical bar plus horizontal lines out to each child.
function PhyloTree({ tree, activeGeneName }) {
  const ROW_HEIGHT = 22;
  const X_SCALE = 900; // branch-length units are small (k-mer distances ~0-1), scale up for visibility
  const LABEL_PAD = 8;
  const LEFT_MARGIN = 12;

  const leaves = [];
  const edges = [];
  let leafCounter = 0;

  // Two-pass layout: first assign x by cumulative branch length (in tree units),
  // then assign y by leaf order, then propagate y up to internal nodes.
  function assignX(node, cumulative) {
    const x = cumulative + (node.length || 0);
    node._x = x;
    if (node.children) {
      node.children.forEach((c) => assignX(c, x));
    } else {
      node._leafIndex = leafCounter++;
      leaves.push(node);
    }
  }
  assignX(tree, 0);

  function assignY(node) {
    if (!node.children) {
      node._y = node._leafIndex * ROW_HEIGHT + ROW_HEIGHT / 2;
      return node._y;
    }
    const childYs = node.children.map(assignY);
    node._y = (Math.min(...childYs) + Math.max(...childYs)) / 2;
    return node._y;
  }
  assignY(tree);

  function collectEdges(node) {
    if (!node.children) return;
    const childXs = node.children.map((c) => c._x);
    const childYs = node.children.map((c) => c._y);
    edges.push({ type: "v", x: LEFT_MARGIN + node._x * X_SCALE, y1: Math.min(...childYs), y2: Math.max(...childYs) });
    node.children.forEach((c) => {
      edges.push({ type: "h", x1: LEFT_MARGIN + node._x * X_SCALE, x2: LEFT_MARGIN + c._x * X_SCALE, y: c._y });
      collectEdges(c);
    });
  }
  collectEdges(tree);

  const maxX = Math.max(...leaves.map((l) => l._x)) * X_SCALE + LEFT_MARGIN;
  const width = maxX + 220;
  const height = leaves.length * ROW_HEIGHT + 10;

  return (
    <svg viewBox={`0 0 ${width} ${height}`} width="100%" height={height} style={{ overflow: "visible" }}>
      {edges.map((e, i) =>
        e.type === "v" ? (
          <line key={i} x1={e.x} x2={e.x} y1={e.y1} y2={e.y2} stroke="#2A3B5C" strokeWidth={1.5} />
        ) : (
          <line key={i} x1={e.x1} x2={e.x2} y1={e.y} y2={e.y} stroke="#2A3B5C" strokeWidth={1.5} />
        )
      )}
      {leaves.map((leaf, i) => {
        const isActive = leaf.name === activeGeneName;
        return (
          <g key={i}>
            <circle cx={LEFT_MARGIN + leaf._x * X_SCALE} cy={leaf._y} r={3} fill={isActive ? "#E8A33D" : "#3FA796"} />
            <text
              x={LEFT_MARGIN + leaf._x * X_SCALE + LABEL_PAD}
              y={leaf._y + 4}
              fontSize={11}
              fontFamily="'IBM Plex Mono', monospace"
              fill={isActive ? "#E8A33D" : "#A8B3C4"}
              fontWeight={isActive ? 600 : 400}
            >
              {leaf.name}
            </text>
          </g>
        );
      })}
    </svg>
  );
}

function TextField({ label, value, onChange, placeholder, disabled, width }) {
  return (
    <div style={{ display: "flex", flexDirection: "column", gap: 5, width }}>
      <span style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em" }}>
        {label}
      </span>
      <input
        type="text"
        value={value}
        onChange={(e) => onChange(e.target.value)}
        placeholder={placeholder}
        disabled={disabled}
        style={{
          padding: "8px 10px",
          borderRadius: 6,
          border: "1px solid #2A3B5C",
          background: "#0B1220",
          color: "#E7ECF3",
          fontSize: 13,
          fontFamily: "'IBM Plex Mono', monospace",
          opacity: disabled ? 0.5 : 1,
        }}
      />
    </div>
  );
}

function ActionButton({ onClick, disabled, busy, icon: Icon, children }) {
  return (
    <button
      onClick={onClick}
      disabled={disabled}
      style={{
        display: "flex",
        alignItems: "center",
        gap: 8,
        padding: "9px 18px",
        borderRadius: 6,
        border: "none",
        background: disabled ? "#2A3B5C" : "#E8A33D",
        color: disabled ? "#6B7A93" : "#0B1220",
        fontSize: 13,
        fontWeight: 600,
        fontFamily: "'IBM Plex Mono', monospace",
        cursor: disabled ? "not-allowed" : "pointer",
        whiteSpace: "nowrap",
        height: 38,
      }}
    >
      {busy ? <Loader2 size={14} className="spin" /> : <Icon size={14} />}
      {children}
    </button>
  );
}

export default function AMRDashboard() {
  // ── Reference panel setup state ──
  const [panelOpen, setPanelOpen] = useState(false);
  const [taxon, setTaxon] = useState("");
  const [limit, setLimit] = useState("30");
  const [assemblyLevel, setAssemblyLevel] = useState("complete");
  const fetchJob = useJobRunner();

  const [buildDirectory, setBuildDirectory] = useState("");
  const [buildDetector, setBuildDetector] = useState("rgi");
  const [buildKmer, setBuildKmer] = useState("21");
  const buildJob = useJobRunner();

  // Auto-fill the build directory once genomes have been fetched.
  useEffect(() => {
    if (fetchJob.status === "done" && fetchJob.result?.genome_dir) {
      setBuildDirectory(fetchJob.result.genome_dir);
    }
  }, [fetchJob.status, fetchJob.result]);

  // ── Single-genome analysis state ──
  const [file, setFile] = useState(null);
  const [detector, setDetector] = useState("rgi");
  const analyzeJob = useJobRunner();

  // Optional: a pre-computed RGI/AMRFinder output file, generated locally
  // (e.g. `rgi main` run on your own machine where memory isn't
  // constrained). When set, the server skips running the detector itself
  // — no RGI/AMRFinder subprocess, no memory spike on Render — and just
  // parses this file, then runs the lightweight downstream steps.
  const [precomputedFile, setPrecomputedFile] = useState(null);

  const [panelMeta, setPanelMeta] = useState(null);
  const [selected, setSelected] = useState(null);
  const [hovered, setHovered] = useState(null);

  // ── Structural analysis state (on-demand, per selected gene) ──
  const [structureResult, setStructureResult] = useState(null);
  const [structureLoading, setStructureLoading] = useState(false);
  const [structureError, setStructureError] = useState(null);

  const [pocketsResult, setPocketsResult] = useState(null);
  const [pocketsLoading, setPocketsLoading] = useState(false);
  const [pocketsError, setPocketsError] = useState(null);

  const [similarityResult, setSimilarityResult] = useState(null);
  const [similarityLoading, setSimilarityLoading] = useState(false);
  const [similarityError, setSimilarityError] = useState(null);

  const [dockingResult, setDockingResult] = useState(null);
  const [customLigandInput, setCustomLigandInput] = useState("");
  const [dockingLoading, setDockingLoading] = useState(false);
  const [dockingError, setDockingError] = useState(null);
  const [viewingPose, setViewingPose] = useState(null); // { ligand, pdbqt } | null

  const [pdbValidationResult, setPdbValidationResult] = useState(null);
  const [pdbValidationLoading, setPdbValidationLoading] = useState(false);
  const [pdbValidationError, setPdbValidationError] = useState(null);
  const [viewingPdbId, setViewingPdbId] = useState(null);

  const [chatOpen, setChatOpen] = useState(false);
  const [chatMessages, setChatMessages] = useState([]); // {role, content}
  const [chatInput, setChatInput] = useState("");
  const [chatLoading, setChatLoading] = useState(false);
  const [chatError, setChatError] = useState(null);
  const [availableModels, setAvailableModels] = useState([]);
  const [selectedModel, setSelectedModel] = useState("");
  const [modelsLoading, setModelsLoading] = useState(false);

  const fileInputRef = useRef(null);
  const precomputedFileInputRef = useRef(null);

  const genes = analyzeJob.result?.genes || [];
  const islands = analyzeJob.result?.islands || [];
  const detectorMeta = analyzeJob.result?.detector || null;
  const active = hovered || selected || genes[0] || null;

  // Other detected genes sharing the same CARD AMR gene family — a real
  // evolutionary/functional grouping (e.g. "RND efflux pump family"),
  // computed entirely client-side from data already in `genes`. No new
  // API call, no new backend dependency.
  const familySiblings =
    active?.gene_family
      ? genes.filter(
          (g) => g.gene_family === active.gene_family && !(g.id === active.id && g.pos === active.pos)
        )
      : [];

  const genomeLen =
    genes.length > 0 ? Math.max(...genes.map((g) => g.pos || 0)) * 1.15 : null;

  // Genome-level predicted resistance profile — aggregates every detected
  // gene's curated drug-class association into a per-drug-class summary.
  // This is standard genotype-based phenotype inference (same approach
  // used by ResFinder/PATRIC-style tools), NOT a black-box ML prediction
  // of MIC values — it answers "which drug classes have gene-level
  // evidence of resistance" using data already fetched, no new API calls.
  const drugProfile = (() => {
    const byClass = {};
    for (const g of genes) {
      if (!g.drug) continue;
      const classes = g.drug.split(";").map((d) => d.trim()).filter(Boolean);
      for (const cls of classes) {
        if (!byClass[cls]) byClass[cls] = { genes: [], hgtCount: 0, clonalCount: 0 };
        byClass[cls].genes.push(g);
        if (g.klass === "HGT") byClass[cls].hgtCount += 1;
        else if (g.klass === "CLONAL") byClass[cls].clonalCount += 1;
      }
    }
    return Object.entries(byClass)
      .map(([drugClass, data]) => ({ drugClass, ...data }))
      .sort((a, b) => b.genes.length - a.genes.length);
  })();

  // Fresh gene selection should not carry over the previous gene's structure.
  useEffect(() => {
    setStructureResult(null);
    setStructureError(null);
    setStructureLoading(false);
    setPocketsResult(null);
    setPocketsError(null);
    setPocketsLoading(false);
    setSimilarityResult(null);
    setSimilarityError(null);
    setSimilarityLoading(false);
    setDockingResult(null);
    setDockingError(null);
    setDockingLoading(false);
    setViewingPose(null);
    setPdbValidationResult(null);
    setPdbValidationError(null);
    setPdbValidationLoading(false);
    setViewingPdbId(null);
  }, [active?.id, active?.pos]);

  function handleAnalyzeStructure() {
    if (!active) return;
    setStructureLoading(true);
    setStructureError(null);
    setStructureResult(null);
    fetch(`${API_BASE}/alphafold/gene/${encodeURIComponent(active.id)}`)
      .then((r) => r.json())
      .then((data) => {
        setStructureResult(data);
        setStructureLoading(false);
      })
      .catch((err) => {
        setStructureError(err.message);
        setStructureLoading(false);
      });
  }

  function handleDetectPockets() {
    if (!structureResult?.available || !active) return;
    setPocketsLoading(true);
    setPocketsError(null);
    setPocketsResult(null);
    const url = `${API_BASE}/pockets/gene/${encodeURIComponent(active.id)}?pdb_url=${encodeURIComponent(structureResult.pdb_url)}`;
    fetch(url)
      .then((r) => r.json())
      .then((data) => {
        setPocketsResult(data);
        setPocketsLoading(false);
      })
      .catch((err) => {
        setPocketsError(err.message);
        setPocketsLoading(false);
      });
  }

  function handleFindSimilar() {
    if (!structureResult?.available || !active) return;
    setSimilarityLoading(true);
    setSimilarityError(null);
    setSimilarityResult(null);
    const url = `${API_BASE}/similarity/gene/${encodeURIComponent(active.id)}?uniprot=${encodeURIComponent(structureResult.uniprot)}`;
    fetch(url)
      .then((r) => r.json())
      .then((data) => {
        setSimilarityResult(data);
        setSimilarityLoading(false);
      })
      .catch((err) => {
        setSimilarityError(err.message);
        setSimilarityLoading(false);
      });
  }

  function handleRunDocking() {
    if (!structureResult?.available || !active) return;
    setDockingLoading(true);
    setDockingError(null);
    setDockingResult(null);
    const params = new URLSearchParams({
      pdb_url: structureResult.pdb_url,
      drug_class: active.drug || "",
      gene_family: active.gene_family || "",
    });
    if (customLigandInput.trim()) {
      params.set("custom_ligands", customLigandInput.trim());
    }
    fetch(`${API_BASE}/docking/gene/${encodeURIComponent(active.id)}?${params}`)
      .then((r) => r.json())
      .then((data) => {
        setDockingResult(data);
        setDockingLoading(false);
      })
      .catch((err) => {
        setDockingError(err.message);
        setDockingLoading(false);
      });
  }

  function handleCheckPdb() {
    if (!structureResult?.available || !active) return;
    setPdbValidationLoading(true);
    setPdbValidationError(null);
    setPdbValidationResult(null);
    const url = `${API_BASE}/validation/pdb/${encodeURIComponent(active.id)}?uniprot=${encodeURIComponent(structureResult.uniprot)}`;
    fetch(url)
      .then((r) => r.json())
      .then((data) => {
        setPdbValidationResult(data);
        setPdbValidationLoading(false);
      })
      .catch((err) => {
        setPdbValidationError(err.message);
        setPdbValidationLoading(false);
      });
  }

  function fetchAvailableModels() {
    if (availableModels.length > 0 || modelsLoading) return; // only fetch once per session
    setModelsLoading(true);
    fetch(`${API_BASE}/chat/models`)
      .then((r) => r.json())
      .then((data) => {
        if (data.available && data.models?.length > 0) {
          setAvailableModels(data.models);
          if (!selectedModel) setSelectedModel(data.models[0].id);
        }
        setModelsLoading(false);
      })
      .catch(() => setModelsLoading(false));
  }

  function handleSendChat() {
    const question = chatInput.trim();
    if (!question || chatLoading) return;

    const newMessages = [...chatMessages, { role: "user", content: question }];
    setChatMessages(newMessages);
    setChatInput("");
    setChatLoading(true);
    setChatError(null);

    fetch(`${API_BASE}/chat`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        question,
        genes,
        conversation_history: chatMessages, // prior turns, not including this new question yet
        active_gene_id: active?.id || null,
        model: selectedModel || undefined,
        structural_context: {
          structure: structureResult || null,
          pockets: pocketsResult || null,
          docking: dockingResult || null,
          pdb_validation: pdbValidationResult || null,
        },
      }),
    })
      .then((r) => r.json())
      .then((data) => {
        if (data.available) {
          setChatMessages([...newMessages, { role: "assistant", content: data.answer }]);
        } else {
          setChatError(data.error || "Something went wrong");
        }
        setChatLoading(false);
      })
      .catch((err) => {
        setChatError(err.message);
        setChatLoading(false);
      });
  }

  function refreshPanelMeta() {
    fetch(`${API_BASE}/panel/meta`)
      .then((r) => r.json())
      .then(setPanelMeta)
      .catch(() => setPanelMeta(null));
  }

  useEffect(() => {
    refreshPanelMeta();
  }, []);

  // Once a panel build finishes, pull the fresh meta so the header updates.
  useEffect(() => {
    if (buildJob.status === "done") refreshPanelMeta();
  }, [buildJob.status]);

  function handleFetchGenomes() {
    if (!taxon || !limit) return;
    fetchJob.run(async () => {
      const res = await fetch(`${API_BASE}/genomes/fetch`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ taxon, limit: Number(limit), assembly_level: assemblyLevel }),
      });
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || "Failed to start genome fetch");
      return data.job_id;
    });
  }

  function handleBuildPanel() {
    if (!buildDirectory) return;
    buildJob.run(async () => {
      const res = await fetch(`${API_BASE}/panel/build`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          directory: buildDirectory,
          detector: buildDetector,
          kmer: Number(buildKmer),
        }),
      });
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || "Failed to start panel build");
      return data.job_id;
    });
  }

  function handleAnalyze() {
    if (!file) return;
    setSelected(null);
    analyzeJob.run(async () => {
      const formData = new FormData();
      formData.append("file", file);
      const uploadRes = await fetch(`${API_BASE}/upload`, { method: "POST", body: formData });
      const uploadData = await uploadRes.json();
      if (!uploadRes.ok) throw new Error(uploadData.error || "Upload failed");

      // If a pre-computed RGI/AMRFinder output was provided, upload it
      // too — the backend will detect it's attached to this job_id and
      // skip running the detector itself.
      if (precomputedFile) {
        const detectorFormData = new FormData();
        detectorFormData.append("file", precomputedFile);
        const detectorUploadRes = await fetch(
          `${API_BASE}/jobs/${uploadData.job_id}/upload-detector-output`,
          { method: "POST", body: detectorFormData }
        );
        const detectorUploadData = await detectorUploadRes.json();
        if (!detectorUploadRes.ok) {
          throw new Error(detectorUploadData.error || "Failed to upload detector output");
        }
      }

      const startRes = await fetch(`${API_BASE}/jobs/${uploadData.job_id}/start`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ detector }),
      });
      if (!startRes.ok) throw new Error("Failed to start job");
      return uploadData.job_id;
    });
  }

  const hgtCount = genes.filter((g) => g.klass === "HGT").length;
  const isAnalyzeBusy = analyzeJob.status === "starting" || analyzeJob.status === "running";
  const isFetchBusy = fetchJob.status === "starting" || fetchJob.status === "running";
  const isBuildBusy = buildJob.status === "starting" || buildJob.status === "running";

  return (
    <div
      style={{
        minHeight: "100vh",
        background: "#0B1220",
        color: "#E7ECF3",
        fontFamily: "'IBM Plex Sans', sans-serif",
      }}
    >
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500;600&display=swap');
        * { box-sizing: border-box; }
        body { margin: 0; }
        .row:hover { background: rgba(255,255,255,0.03) !important; }
        .track-tick { transition: transform 0.15s ease, filter 0.15s ease; cursor: pointer; }
        .track-tick:hover { transform: scaleY(1.4); filter: brightness(1.3); }
        .spin { animation: spin 1s linear infinite; }
        @keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }
        .detector-btn { cursor: pointer; transition: all 0.15s ease; }
        select.dark-select { appearance: none; }
      `}</style>

      {/* Header */}
      <div style={{ borderBottom: "1px solid #1B2740", padding: "20px 32px" }}>
        <div style={{ maxWidth: 1100, margin: "0 auto", display: "flex", alignItems: "center", justifyContent: "space-between" }}>
          <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
            <Dna size={20} color="#E8A33D" strokeWidth={1.75} />
            <span style={{ fontFamily: "'IBM Plex Mono', monospace", fontSize: 15, fontWeight: 600, letterSpacing: "0.01em" }}>
              KALI · AMR-HGT
            </span>
          </div>
          <div style={{ fontSize: 13, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace" }}>
            {panelMeta?.taxon_scope
              ? `${panelMeta.taxon_scope} · panel ${panelMeta.version || "—"} (${panelMeta.n_genomes ?? 0} genomes)`
              : "Panel not loaded"}
          </div>
        </div>
      </div>

      <div style={{ maxWidth: 1100, margin: "0 auto", padding: "40px 32px 80px" }}>
        {/* ── Reference Panel Setup (collapsible) ── */}
        <div
          style={{
            border: "1px solid #1B2740",
            borderRadius: 8,
            background: "#0F1830",
            marginBottom: 36,
            overflow: "hidden",
          }}
        >
          <button
            onClick={() => setPanelOpen((v) => !v)}
            style={{
              width: "100%",
              display: "flex",
              alignItems: "center",
              justifyContent: "space-between",
              padding: "16px 20px",
              background: "none",
              border: "none",
              cursor: "pointer",
              color: "#E7ECF3",
            }}
          >
            <div style={{ display: "flex", alignItems: "center", gap: 10 }}>
              <Layers size={16} color="#E8A33D" />
              <span style={{ fontFamily: "'IBM Plex Mono', monospace", fontSize: 13, fontWeight: 600 }}>
                Reference Panel Setup
              </span>
              <span style={{ fontSize: 12, color: "#6B7A93" }}>
                — download genomes and (re)build the panel, no terminal needed
              </span>
            </div>
            {panelOpen ? <ChevronUp size={16} color="#6B7A93" /> : <ChevronDown size={16} color="#6B7A93" />}
          </button>

          {panelOpen && (
            <div style={{ padding: "0 20px 24px", display: "flex", flexDirection: "column", gap: 24 }}>
              {/* Step 1: Fetch genomes */}
              <div style={{ borderTop: "1px solid #1B2740", paddingTop: 20 }}>
                <div style={{ fontSize: 12, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", marginBottom: 12, textTransform: "uppercase", letterSpacing: "0.04em" }}>
                  Step 1 — Download genomes
                </div>
                <div style={{ display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end" }}>
                  <TextField label="Taxon" value={taxon} onChange={setTaxon} placeholder="Klebsiella pneumoniae" disabled={isFetchBusy} width={260} />
                  <TextField label="How many" value={limit} onChange={setLimit} placeholder="30" disabled={isFetchBusy} width={100} />
                  <div style={{ display: "flex", flexDirection: "column", gap: 5 }}>
                    <span style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em" }}>
                      Assembly level
                    </span>
                    <select
                      className="dark-select"
                      value={assemblyLevel}
                      onChange={(e) => setAssemblyLevel(e.target.value)}
                      disabled={isFetchBusy}
                      style={{
                        padding: "8px 10px",
                        borderRadius: 6,
                        border: "1px solid #2A3B5C",
                        background: "#0B1220",
                        color: "#E7ECF3",
                        fontSize: 13,
                        fontFamily: "'IBM Plex Mono', monospace",
                      }}
                    >
                      <option value="complete">complete</option>
                      <option value="chromosome">chromosome</option>
                      <option value="scaffold">scaffold</option>
                      <option value="contig">contig</option>
                    </select>
                  </div>
                  <ActionButton onClick={handleFetchGenomes} disabled={!taxon || !limit || isFetchBusy} busy={isFetchBusy} icon={Download}>
                    {isFetchBusy ? "Fetching" : "Fetch genomes"}
                  </ActionButton>
                </div>
                <ProgressLine status={fetchJob.status} stage={fetchJob.stage} error={fetchJob.error} />
                {fetchJob.status === "done" && fetchJob.result && (
                  <div style={{ marginTop: 14, fontSize: 13, color: "#3FA796", display: "flex", alignItems: "center", gap: 8 }}>
                    <CheckCircle2 size={14} />
                    Downloaded {fetchJob.result.downloaded} genome{fetchJob.result.downloaded !== 1 ? "s" : ""} to{" "}
                    <span style={{ fontFamily: "'IBM Plex Mono', monospace", color: "#A8B3C4" }}>{fetchJob.result.genome_dir}</span>
                  </div>
                )}
              </div>

              {/* Step 2: Build panel */}
              <div style={{ borderTop: "1px solid #1B2740", paddingTop: 20 }}>
                <div style={{ fontSize: 12, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", marginBottom: 12, textTransform: "uppercase", letterSpacing: "0.04em" }}>
                  Step 2 — Build panel from those genomes
                </div>
                <div style={{ display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end" }}>
                  <TextField
                    label="Genome folder"
                    value={buildDirectory}
                    onChange={setBuildDirectory}
                    placeholder="Auto-filled after Step 1, or paste a path"
                    disabled={isBuildBusy}
                    width={340}
                  />
                  <div style={{ display: "flex", alignItems: "flex-end", gap: 8 }}>
                    {["rgi", "amrfinder"].map((opt) => (
                      <button
                        key={opt}
                        className="detector-btn"
                        onClick={() => !isBuildBusy && setBuildDetector(opt)}
                        disabled={isBuildBusy}
                        style={{
                          padding: "8px 12px",
                          borderRadius: 6,
                          fontSize: 12,
                          fontFamily: "'IBM Plex Mono', monospace",
                          border: `1px solid ${buildDetector === opt ? "#E8A33D" : "#2A3B5C"}`,
                          background: buildDetector === opt ? "rgba(232,163,61,0.12)" : "transparent",
                          color: buildDetector === opt ? "#E8A33D" : "#8592A6",
                          height: 38,
                        }}
                      >
                        {opt === "rgi" ? "CARD-RGI" : "AMRFinderPlus"}
                      </button>
                    ))}
                  </div>
                  <TextField label="K-mer" value={buildKmer} onChange={setBuildKmer} disabled={isBuildBusy} width={70} />
                  <ActionButton onClick={handleBuildPanel} disabled={!buildDirectory || isBuildBusy} busy={isBuildBusy} icon={Layers}>
                    {isBuildBusy ? "Building" : "Build panel"}
                  </ActionButton>
                </div>
                <p style={{ fontSize: 12, color: "#4B5A73", marginTop: 10, marginBottom: 0 }}>
                  This can take 30–60+ minutes for dozens of genomes — safe to leave running in the background.
                </p>
                <ProgressLine status={buildJob.status} stage={buildJob.stage} error={buildJob.error} />
                {buildJob.status === "done" && buildJob.result && (
                  <div style={{ marginTop: 14, fontSize: 13, color: "#3FA796", display: "flex", alignItems: "center", gap: 8 }}>
                    <CheckCircle2 size={14} />
                    Panel {buildJob.result.version} built — {buildJob.result.n_genomes} genomes,{" "}
                    {buildJob.result.n_genes_classified} genes classified,{" "}
                    {buildJob.result.n_cooccurrence_pairs} co-occurrence pairs.
                  </div>
                )}
              </div>
            </div>
          )}
        </div>

        {/* ── Single-genome upload panel ── */}
        <div
          style={{
            border: "1px solid #1B2740",
            borderRadius: 8,
            padding: 24,
            background: "#0F1830",
            marginBottom: 36,
          }}
        >
          <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", gap: 20, flexWrap: "wrap" }}>
            <div style={{ display: "flex", alignItems: "center", gap: 14, flex: 1, minWidth: 260 }}>
              <input
                ref={fileInputRef}
                type="file"
                accept=".fasta,.fa,.fna"
                onChange={(e) => setFile(e.target.files?.[0] || null)}
                style={{ display: "none" }}
              />
              <button
                onClick={() => fileInputRef.current?.click()}
                disabled={isAnalyzeBusy}
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: 8,
                  padding: "9px 16px",
                  borderRadius: 6,
                  border: "1px solid #2A3B5C",
                  background: "#131C2E",
                  color: "#E7ECF3",
                  fontSize: 13,
                  fontFamily: "'IBM Plex Mono', monospace",
                  cursor: isAnalyzeBusy ? "not-allowed" : "pointer",
                  opacity: isAnalyzeBusy ? 0.5 : 1,
                }}
              >
                <Upload size={14} />
                Choose FASTA
              </button>
              <span style={{ fontSize: 13, color: file ? "#E7ECF3" : "#4B5A73", fontFamily: "'IBM Plex Mono', monospace" }}>
                {file ? file.name : "No file selected"}
              </span>
            </div>

            <div style={{ display: "flex", alignItems: "center", gap: 14 }}>
              <input
                ref={precomputedFileInputRef}
                type="file"
                accept=".txt,.tsv"
                onChange={(e) => setPrecomputedFile(e.target.files?.[0] || null)}
                style={{ display: "none" }}
              />
              <button
                onClick={() => precomputedFileInputRef.current?.click()}
                disabled={isAnalyzeBusy}
                title="Optional: upload an already-computed RGI/AMRFinder output file (e.g. run locally) to skip running the detector on the server"
                style={{
                  display: "flex",
                  alignItems: "center",
                  gap: 8,
                  padding: "9px 16px",
                  borderRadius: 6,
                  border: "1px solid #2A3B5C",
                  background: precomputedFile ? "rgba(232,163,61,0.12)" : "#131C2E",
                  color: precomputedFile ? "#E8A33D" : "#E7ECF3",
                  fontSize: 13,
                  fontFamily: "'IBM Plex Mono', monospace",
                  cursor: isAnalyzeBusy ? "not-allowed" : "pointer",
                  opacity: isAnalyzeBusy ? 0.5 : 1,
                }}
              >
                <Upload size={14} />
                Detector output (optional)
              </button>
              <span style={{ fontSize: 13, color: precomputedFile ? "#E7ECF3" : "#4B5A73", fontFamily: "'IBM Plex Mono', monospace" }}>
                {precomputedFile ? precomputedFile.name : "Run RGI/AMRFinder server-side"}
              </span>
              {precomputedFile && (
                <button
                  onClick={() => { setPrecomputedFile(null); if (precomputedFileInputRef.current) precomputedFileInputRef.current.value = ""; }}
                  disabled={isAnalyzeBusy}
                  style={{
                    fontSize: 12,
                    fontFamily: "'IBM Plex Mono', monospace",
                    color: "#8592A6",
                    background: "transparent",
                    border: "none",
                    cursor: isAnalyzeBusy ? "not-allowed" : "pointer",
                    textDecoration: "underline",
                  }}
                >
                  Clear
                </button>
              )}
            </div>

            <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
              <span style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em" }}>
                Detector
              </span>
              {["rgi", "amrfinder"].map((opt) => (
                <button
                  key={opt}
                  className="detector-btn"
                  onClick={() => !isAnalyzeBusy && setDetector(opt)}
                  disabled={isAnalyzeBusy}
                  style={{
                    padding: "6px 12px",
                    borderRadius: 6,
                    fontSize: 12,
                    fontFamily: "'IBM Plex Mono', monospace",
                    border: `1px solid ${detector === opt ? "#E8A33D" : "#2A3B5C"}`,
                    background: detector === opt ? "rgba(232,163,61,0.12)" : "transparent",
                    color: detector === opt ? "#E8A33D" : "#8592A6",
                    opacity: isAnalyzeBusy ? 0.5 : 1,
                  }}
                >
                  {opt === "rgi" ? "CARD-RGI" : "AMRFinderPlus"}
                </button>
              ))}
            </div>

            <ActionButton onClick={handleAnalyze} disabled={!file || isAnalyzeBusy} busy={isAnalyzeBusy} icon={Search}>
              {isAnalyzeBusy ? "Running" : "Analyze"}
            </ActionButton>
          </div>

          <ProgressLine status={analyzeJob.status} stage={analyzeJob.stage} error={analyzeJob.error} />
        </div>

        {/* Empty state — before any upload */}
        {analyzeJob.status === "idle" && genes.length === 0 && (
          <div style={{ textAlign: "center", padding: "60px 20px", color: "#4B5A73" }}>
            <Dna size={32} strokeWidth={1.2} style={{ marginBottom: 12, opacity: 0.6 }} />
            <p style={{ fontSize: 14 }}>Upload a genome to detect AMR genes and classify their origin.</p>
          </div>
        )}

        {/* A completed analysis that genuinely found zero AMR genes — a real
            result, not an error, distinct from the pre-upload empty state above. */}
        {analyzeJob.status === "done" && genes.length === 0 && (
          <div style={{ textAlign: "center", padding: "60px 20px", color: "#4B5A73" }}>
            <Dna size={32} strokeWidth={1.2} style={{ marginBottom: 12, opacity: 0.6 }} />
            <p style={{ fontSize: 14 }}>No AMR genes detected in this genome.</p>
          </div>
        )}

        {/* Results */}
        {genes.length > 0 && (
          <>
            <div style={{ marginBottom: 36 }}>
              <h1 style={{ fontFamily: "'IBM Plex Mono', monospace", fontSize: 26, fontWeight: 600, margin: 0, lineHeight: 1.3 }}>
                {genes.length} resistance gene{genes.length !== 1 ? "s" : ""} found.{" "}
                <span style={{ color: "#5FA777" }}>{hgtCount} arrived by horizontal transfer.</span>
              </h1>
              <p style={{ color: "#8592A6", fontSize: 14, marginTop: 8, maxWidth: 620 }}>
                Every gene below is scored for compositional anomaly against its local genomic
                background, then classified as clonally inherited or horizontally acquired.
                {detectorMeta && (
                  <> Detected with <strong style={{ color: "#A8B3C4" }}>{detectorMeta.label}</strong>.</>
                )}
                {" "}
                <span style={{ color: islands.length > 0 ? "#3FA796" : "#6B7A93" }}>
                  {islands.length} compositional island{islands.length !== 1 ? "s" : ""} detected genome-wide
                  {islands.length > 0 && ` (up to Z=${Math.max(...islands.map((i) => i.max_zscore || 0)).toFixed(1)})`}.
                </span>
              </p>
            </div>

            {drugProfile.length > 0 && (
              <div style={{ marginBottom: 44 }}>
                <div style={{ fontSize: 12, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", marginBottom: 10, textTransform: "uppercase", letterSpacing: "0.06em" }}>
                  Predicted resistance profile — {drugProfile.length} drug class{drugProfile.length !== 1 ? "es" : ""} implicated
                </div>
                <p style={{ fontSize: 12, color: "#4B5A73", marginBottom: 12, maxWidth: 620 }}>
                  Genotype-based inference from detected genes' curated drug-class associations —
                  not a predicted MIC or clinical outcome. HGT-acquired genes are a stronger
                  resistance signal than intrinsic/regulatory ones.
                </p>
                <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                  {drugProfile.map(({ drugClass, genes: classGenes, hgtCount: h, clonalCount: c }) => (
                    <div
                      key={drugClass}
                      style={{
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "space-between",
                        padding: "8px 12px",
                        borderRadius: 6,
                        border: "1px solid #1B2740",
                        background: "#0F1830",
                        fontSize: 13,
                      }}
                    >
                      <span style={{ color: "#E7ECF3" }}>{drugClass}</span>
                      <span style={{ display: "flex", gap: 10, fontSize: 11, fontFamily: "'IBM Plex Mono', monospace" }}>
                        <span style={{ color: "#8592A6" }}>{classGenes.length} gene{classGenes.length !== 1 ? "s" : ""}</span>
                        {h > 0 && <span style={{ color: "#5FA777" }}>{h} HGT</span>}
                        {c > 0 && <span style={{ color: "#A8825E" }}>{c} clonal</span>}
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            )}

            {genomeLen && (
              <div style={{ marginBottom: 44 }}>
                <div style={{ fontSize: 12, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", marginBottom: 10, textTransform: "uppercase", letterSpacing: "0.06em" }}>
                  Genome track ({(genomeLen / 1e6).toFixed(2)} Mb, approx.)
                </div>
                <div style={{ position: "relative", height: 64, background: "#0F1830", border: "1px solid #1B2740", borderRadius: 6 }}>
                  <div style={{ position: "absolute", top: "50%", left: 0, right: 0, height: 1, background: "#1B2740" }} />
                  {islands.map((isl, i) => {
                    if (isl.start == null || isl.end == null) return null;
                    const left = (isl.start / genomeLen) * 100;
                    const width = Math.max(((isl.end - isl.start) / genomeLen) * 100, 0.3);
                    return (
                      <div
                        key={`island-${i}`}
                        title={`Island: ${isl.start.toLocaleString()}\u2013${isl.end.toLocaleString()} bp (Z=${isl.max_zscore?.toFixed(2)})`}
                        style={{
                          position: "absolute",
                          left: `${left}%`,
                          width: `${width}%`,
                          top: 0,
                          bottom: 0,
                          background: "rgba(63,167,150,0.14)",
                          borderLeft: "1px solid rgba(63,167,150,0.4)",
                          borderRight: "1px solid rgba(63,167,150,0.4)",
                        }}
                      />
                    );
                  })}
                  {genes.map((g, i) => {
                    if (!g.pos) return null;
                    const left = (g.pos / genomeLen) * 100;
                    const isHgt = g.klass === "HGT";
                    const scoreForHeight = g.score != null ? Math.min(g.score / 5, 1) : 0.3;
                    const height = 10 + scoreForHeight * 34;
                    return (
                      <div
                        key={`${g.id}-${i}`}
                        className="track-tick"
                        onMouseEnter={() => setHovered(g)}
                        onMouseLeave={() => setHovered(null)}
                        onClick={() => setSelected(g)}
                        style={{
                          position: "absolute",
                          left: `${left}%`,
                          top: "50%",
                          transform: "translate(-50%, -50%)",
                          width: 3,
                          height,
                          background: isHgt ? "#5FA777" : "#A8825E",
                          borderRadius: 2,
                          boxShadow: isHgt ? "0 0 8px rgba(95,167,119,0.5)" : "none",
                        }}
                        title={g.id}
                      />
                    );
                  })}
                </div>
                <div style={{ display: "flex", justifyContent: "space-between", marginTop: 6, fontSize: 11, color: "#4B5A73", fontFamily: "'IBM Plex Mono', monospace" }}>
                  <span>0 bp</span>
                  <span>{Math.round(genomeLen).toLocaleString()} bp</span>
                </div>
                {islands.length > 0 && (
                  <div style={{ display: "flex", alignItems: "center", gap: 6, marginTop: 8, fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace" }}>
                    <span style={{ display: "inline-block", width: 12, height: 10, background: "rgba(63,167,150,0.25)", border: "1px solid rgba(63,167,150,0.5)" }} />
                    Compositional anomaly island
                  </div>
                )}
              </div>
            )}

            <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
              <div>
                <div style={{ fontSize: 12, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", marginBottom: 10, textTransform: "uppercase", letterSpacing: "0.06em" }}>
                  Detected genes ({genes.length})
                </div>
                <select
                  value={genes.findIndex((g) => g.id === active?.id && g.pos === active?.pos)}
                  onChange={(e) => setSelected(genes[Number(e.target.value)])}
                  style={{
                    width: "100%",
                    padding: "12px 14px",
                    fontSize: 14,
                    fontFamily: "'IBM Plex Mono', monospace",
                    background: "#0F1830",
                    color: "#E7ECF3",
                    border: "1px solid #1B2740",
                    borderRadius: 8,
                    cursor: "pointer",
                  }}
                >
                  <option value={-1} disabled style={{ color: "#4B5A73" }}>
                    Select a gene\u2026
                  </option>
                  {genes.map((g, i) => (
                    <option
                      key={`${g.id}-${i}`}
                      value={i}
                      style={{
                        color: g.klass === "HGT" ? "#5FA777" : g.klass === "CLONAL" ? "#A8825E" : "#6B7A93",
                        background: "#0F1830",
                      }}
                    >
                      {g.klass === "HGT" ? "\u25cf HGT" : g.klass === "CLONAL" ? "\u25cf CLONAL" : "\u25cb AMBIGUOUS"} \u2014 {g.id}
                    </option>
                  ))}
                </select>
                <p style={{ fontSize: 11, color: "#4B5A73", marginTop: 8 }}>
                  Select a gene for drug class, score, and full origin evidence.
                </p>
              </div>

              {active && (
                <div
                  style={{
                    border: "1px solid #1B2740",
                    borderRadius: 8,
                    padding: 20,
                    background: "#0F1830",
                  }}
                >
                  <div style={{ display: "flex", justifyContent: "space-between", alignItems: "flex-start", marginBottom: 6 }}>
                    <div>
                      <div style={{ fontFamily: "'IBM Plex Mono', monospace", fontSize: 16, fontWeight: 600 }}>{active.id}</div>
                      <div style={{ fontSize: 12, color: "#6B7A93", marginTop: 2 }}>{active.drug || "Drug class unknown"}</div>
                    </div>
                    <Badge klass={active.klass} />
                  </div>
                  <div style={{ marginBottom: 14 }}>
                    <EvidenceTag evidence={active.evidence} />
                  </div>

                  <div style={{ display: "flex", flexDirection: "column", gap: 10, fontSize: 13 }}>
                    <Row label="Position" value={active.pos ? `${active.pos.toLocaleString()} bp` : null} />
                    <Row label="Contig" value={active.contig} mono />
                    <Row label="Island Z-score" value={active.score != null ? active.score.toFixed(3) : null} />
                    <Row label="Mantel r" value={active.mantel_r != null ? active.mantel_r.toFixed(3) : null} />
                    {active.gene_family && <Row label="Gene family" value={active.gene_family} />}
                    {active.mechanism && <Row label="Mechanism" value={active.mechanism} />}
                    {active.model_type && <Row label="Detection model" value={active.model_type} />}
                  </div>

                  {/* Evidence & Confidence Summary — every independent
                      signal in one place, each labeled on its own terms
                      rather than blended into one misleading number. */}
                  <div style={{ marginTop: 14, padding: "10px 12px", borderRadius: 6, background: "#0F1830", border: "1px solid #1B2740" }}>
                    <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 8 }}>
                      Evidence &amp; confidence summary
                    </div>

                    {(() => {
                      const cls = classificationConfidence(active.evidence, active.mantel_p_value, active.mantel_n_carriers);
                      return (
                        <div style={{ display: "flex", alignItems: "flex-start", gap: 8, marginBottom: 8 }}>
                          <span style={{ fontSize: 11, color: "#6B7A93", minWidth: 90 }}>Classification</span>
                          <div>
                            <span style={{ fontSize: 12, color: cls.color, fontWeight: 600 }}>{cls.label}</span>
                            {cls.detail && <div style={{ fontSize: 11, color: "#6B7A93", marginTop: 2 }}>{cls.detail}</div>}
                          </div>
                        </div>
                      );
                    })()}

                    {structureResult?.available && (
                      <div style={{ display: "flex", alignItems: "flex-start", gap: 8, marginBottom: 8 }}>
                        <span style={{ fontSize: 11, color: "#6B7A93", minWidth: 90 }}>Structure</span>
                        <div>
                          <span
                            style={{
                              fontSize: 12,
                              fontWeight: 600,
                              color: structureResult.mean_plddt > 90 ? "#1E6FEB" : structureResult.mean_plddt > 70 ? "#65CBF3" : structureResult.mean_plddt > 50 ? "#FFDB13" : "#FF7D45",
                            }}
                          >
                            {structureResult.mean_plddt > 90 ? "Very high" : structureResult.mean_plddt > 70 ? "High" : structureResult.mean_plddt > 50 ? "Low" : "Very low"}
                          </span>
                          <div style={{ fontSize: 11, color: "#6B7A93", marginTop: 2 }}>
                            pLDDT {typeof structureResult.mean_plddt === "number" ? structureResult.mean_plddt.toFixed(1) : structureResult.mean_plddt} — AlphaFold's own confidence in this predicted fold.
                          </div>
                        </div>
                      </div>
                    )}

                    {dockingResult?.available && dockingResult.results?.some((r) => r.available) && (
                      <div style={{ display: "flex", alignItems: "flex-start", gap: 8 }}>
                        <span style={{ fontSize: 11, color: "#6B7A93", minWidth: 90 }}>Docking</span>
                        <div style={{ fontSize: 11, color: "#6B7A93" }}>
                          {dockingResult.results.filter((r) => r.available).map((r, i) => (
                            <div key={i}>
                              {r.ligand}: {r.n_replicates != null ? `${r.n_replicates}/${r.n_replicates_attempted} replicate runs, std ±${r.std_affinity_kcal_mol?.toFixed(1)}` : "single run only"}
                              {" "}— small keyword-matched ligand panel, not a comprehensive screen.
                            </div>
                          ))}
                        </div>
                      </div>
                    )}

                    {!structureResult?.available && !dockingResult?.available && (
                      <div style={{ fontSize: 11, color: "#4B5A73" }}>
                        Structural and docking confidence will appear here once you've run those analyses below.
                      </div>
                    )}
                  </div>

                  {active.mutations && (
                    <div
                      style={{
                        marginTop: 14,
                        padding: "10px 12px",
                        borderRadius: 6,
                        background: "rgba(232,163,61,0.08)",
                        border: "1px solid rgba(232,163,61,0.25)",
                      }}
                    >
                      <div style={{ fontSize: 10, color: "#E8A33D", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 4 }}>
                        Specific mutations detected
                      </div>
                      <div style={{ fontSize: 13, fontFamily: "'IBM Plex Mono', monospace", color: "#E7ECF3" }}>
                        {active.mutations}
                      </div>
                    </div>
                  )}

                  {!active.mutations &&
                    active.model_type &&
                    (active.model_type === "POINT" || active.model_type.toLowerCase().includes("variant")) && (
                      <div style={{ marginTop: 14, fontSize: 12, color: "#8592A6" }}>
                        Flagged as a mutation-based resistance call by the detector, but specific
                        mutation details weren't included in its output.
                      </div>
                    )}

                  {familySiblings.length > 0 && (
                    <div style={{ marginTop: 14, paddingTop: 14, borderTop: "1px solid #1B2740" }}>
                      <div style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 8 }}>
                        Same gene family — {familySiblings.length} other gene{familySiblings.length !== 1 ? "s" : ""} this genome
                      </div>
                      <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
                        {familySiblings.map((g, i) => (
                          <button
                            key={`${g.id}-${i}`}
                            onClick={() => setSelected(g)}
                            style={{
                              padding: "4px 10px",
                              borderRadius: 999,
                              fontSize: 11,
                              fontFamily: "'IBM Plex Mono', monospace",
                              border: "1px solid #2A3B5C",
                              background: "#131C2E",
                              color: "#A8B3C4",
                              cursor: "pointer",
                            }}
                          >
                            {g.id}
                          </button>
                        ))}
                      </div>
                      <p style={{ fontSize: 11, color: "#4B5A73", marginTop: 6, marginBottom: 0 }}>
                        Grouped by CARD's curated gene family ({active.gene_family}) — click to jump to that gene.
                      </p>
                    </div>
                  )}

                  {active.aro_id && (active.ontology_lineage?.length > 0 || active.ontology_siblings?.length > 0) && (
                    <div style={{ marginTop: 14, paddingTop: 14, borderTop: "1px solid #1B2740" }}>
                      <div style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 8 }}>
                        Evolutionary lineage — CARD ontology
                      </div>

                      {active.ontology_lineage?.length > 0 && (
                        <div style={{ display: "flex", flexWrap: "wrap", alignItems: "center", gap: 6, marginBottom: active.ontology_siblings?.length > 0 ? 12 : 0 }}>
                          {[...active.ontology_lineage].reverse().map((term, i) => (
                            <React.Fragment key={term.aro_id || i}>
                              {i > 0 && <ChevronRight size={12} color="#4B5A73" />}
                              <a
                                href={`https://card.mcmaster.ca/aro/${(term.aro_id || "").replace("ARO:", "")}`}
                                target="_blank"
                                rel="noreferrer"
                                title={term.definition || term.name}
                                style={{
                                  fontSize: 12,
                                  color: "#3FA796",
                                  textDecoration: "none",
                                  borderBottom: "1px dotted rgba(63,167,150,0.4)",
                                }}
                              >
                                {term.name}
                              </a>
                            </React.Fragment>
                          ))}
                          <ChevronRight size={12} color="#4B5A73" />
                          <span style={{ fontSize: 12, color: "#E7ECF3", fontWeight: 500 }}>{active.id}</span>
                        </div>
                      )}

                      {active.ontology_siblings?.length > 0 && (
                        <div>
                          <div style={{ fontSize: 10, color: "#4B5A73", marginBottom: 6 }}>
                            Other terms CARD classifies under the same immediate parent
                            (not limited to genes detected in this genome):
                          </div>
                          <div style={{ display: "flex", flexWrap: "wrap", gap: 6 }}>
                            {active.ontology_siblings.map((term, i) => (
                              <a
                                key={term.aro_id || i}
                                href={`https://card.mcmaster.ca/aro/${(term.aro_id || "").replace("ARO:", "")}`}
                                target="_blank"
                                rel="noreferrer"
                                title={term.definition || ""}
                                style={{
                                  padding: "4px 10px",
                                  borderRadius: 999,
                                  fontSize: 11,
                                  fontFamily: "'IBM Plex Mono', monospace",
                                  border: "1px solid #2A3B5C",
                                  background: "#131C2E",
                                  color: "#8592A6",
                                  textDecoration: "none",
                                }}
                              >
                                {term.name}
                              </a>
                            ))}
                          </div>
                        </div>
                      )}
                    </div>
                  )}

                  {active.ontology_error && (
                    <div style={{ marginTop: 10, fontSize: 11, color: "#4B5A73" }}>
                      Ontology lookup unavailable for this gene: {active.ontology_error}
                    </div>
                  )}

                  {active.phylo_tree && (
                    <div style={{ marginTop: 14, paddingTop: 14, borderTop: "1px solid #1B2740" }}>
                      <div style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 4 }}>
                        Phylogenetic placement — protein sequence
                      </div>
                      <p style={{ fontSize: 10, color: "#4B5A73", marginBottom: 10, marginTop: 0 }}>
                        Real tree from actual k-mer sequence distance (alignment-free,
                        same method as this tool's genome-level analysis) — not CARD's
                        curated classification above. Genes here are related by measured
                        sequence similarity, computed only from what's detected in this genome.
                      </p>
                      <div style={{ background: "#0B1220", border: "1px solid #1B2740", borderRadius: 6, padding: "10px 8px", overflowX: "auto" }}>
                        <PhyloTree tree={active.phylo_tree} activeGeneName={active.id} />
                      </div>
                    </div>
                  )}

                  {active.phylo_tree_status && (
                    <div style={{ marginTop: 14, paddingTop: 14, borderTop: "1px solid #1B2740" }}>
                      <div style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                        Phylogenetic placement — protein sequence
                      </div>
                      <div style={{ fontSize: 12, color: "#8592A6" }}>{active.phylo_tree_status}</div>
                    </div>
                  )}

                  {active.phylo_tree_error && (
                    <div style={{ marginTop: 10, fontSize: 11, color: "#4B5A73" }}>
                      Phylogenetic placement unavailable: {active.phylo_tree_error}
                    </div>
                  )}

                  <div
                    style={{
                      marginTop: 18,
                      paddingTop: 14,
                      borderTop: "1px solid #1B2740",
                      display: "flex",
                      gap: 8,
                      fontSize: 12,
                      color: "#6B7A93",
                      lineHeight: 1.5,
                    }}
                  >
                    <Info size={14} style={{ flexShrink: 0, marginTop: 2 }} />
                    <span>
                      {active.evidence === "panel"
                        ? "Classification drawn from the reference panel's Mantel test — this gene has been seen and correlated against known isolate phylogeny."
                        : active.evidence === "panel_normalized"
                        ? "Classification drawn from the reference panel's Mantel test, matched by normalized gene name (e.g. a naming difference like 'blaSHV-11' vs 'SHV-11') rather than an exact match."
                        : "No panel-derived label for this gene yet — classification falls back to this genome's own compositional island overlap only."}
                    </span>
                  </div>

                  {active.interpretation_caveat && (
                    <div
                      style={{
                        marginTop: 10,
                        padding: "10px 12px",
                        borderRadius: 6,
                        background: "rgba(232,93,61,0.08)",
                        border: "1px solid rgba(232,93,61,0.3)",
                        display: "flex",
                        gap: 8,
                      }}
                    >
                      <AlertCircle size={14} style={{ flexShrink: 0, marginTop: 2, color: "#E85D3D" }} />
                      <span style={{ fontSize: 12, color: "#E7ECF3", lineHeight: 1.5 }}>
                        {active.interpretation_caveat}
                      </span>
                    </div>
                  )}

                  {/* Structural Analysis — on-demand, never runs automatically */}
                  <div style={{ marginTop: 18, paddingTop: 14, borderTop: "1px solid #1B2740" }}>
                    <div style={{ fontSize: 11, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 10 }}>
                      Structural analysis
                    </div>

                    {!structureResult && !structureLoading && !structureError && (
                      <button
                        onClick={handleAnalyzeStructure}
                        style={{
                          display: "flex",
                          alignItems: "center",
                          gap: 8,
                          padding: "8px 14px",
                          borderRadius: 6,
                          border: "1px solid #2A3B5C",
                          background: "#131C2E",
                          color: "#E7ECF3",
                          fontSize: 12,
                          fontFamily: "'IBM Plex Mono', monospace",
                          cursor: "pointer",
                        }}
                      >
                        <Atom size={14} color="#E8A33D" />
                        Analyze with AlphaFold
                      </button>
                    )}

                    {structureLoading && (
                      <div style={{ display: "flex", alignItems: "center", gap: 8, fontSize: 12, color: "#8592A6", fontFamily: "'IBM Plex Mono', monospace" }}>
                        <Loader2 size={14} className="spin" color="#E8A33D" />
                        Querying AlphaFold DB...
                      </div>
                    )}

                    {structureError && (
                      <div style={{ display: "flex", alignItems: "center", gap: 8, fontSize: 12, color: "#E85D3D" }}>
                        <AlertCircle size={14} />
                        {structureError}
                      </div>
                    )}

                    {structureResult && !structureResult.available && (
                      <div style={{ fontSize: 12, color: "#8592A6" }}>
                        No structure available — {structureResult.reason || "not found in AlphaFold DB"}.
                      </div>
                    )}

                    {structureResult && structureResult.available && (
                      <div>
                        <div style={{ display: "flex", flexDirection: "column", gap: 8, fontSize: 13, marginBottom: 4 }}>
                          <Row label="Protein" value={structureResult.protein} />
                          <Row label="UniProt" value={structureResult.uniprot} mono />
                          <Row
                            label="pLDDT"
                            value={
                              typeof structureResult.mean_plddt === "number"
                                ? structureResult.mean_plddt.toFixed(1)
                                : structureResult.mean_plddt
                            }
                          />
                        </div>

                        <StructureViewer pdbUrl={structureResult.pdb_url} />

                        {structureResult.active_site_features && structureResult.active_site_features.length > 0 && (
                          <div style={{ marginTop: 12 }}>
                            <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                              Active / binding site residues (UniProt curated)
                            </div>
                            <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                              {structureResult.active_site_features.map((f, i) => (
                                <div
                                  key={i}
                                  style={{
                                    fontSize: 12,
                                    display: "flex",
                                    gap: 8,
                                    color: "#A8B3C4",
                                  }}
                                >
                                  <span style={{ fontFamily: "'IBM Plex Mono', monospace", color: "#E8A33D", minWidth: 42 }}>
                                    {f.position}
                                  </span>
                                  <span>
                                    {f.type}
                                    {f.description ? ` — ${f.description}` : ""}
                                  </span>
                                </div>
                              ))}
                            </div>
                          </div>
                        )}

                        <div style={{ display: "flex", gap: 16, marginTop: 10 }}>
                          <a
                            href={structureResult.pdb_url}
                            target="_blank"
                            rel="noreferrer"
                            style={{ display: "flex", alignItems: "center", gap: 6, fontSize: 12, color: "#E8A33D", fontFamily: "'IBM Plex Mono', monospace", textDecoration: "none" }}
                          >
                            <Download size={12} /> Download PDB
                          </a>
                          {structureResult.cif_url && (
                            <a
                              href={structureResult.cif_url}
                              target="_blank"
                              rel="noreferrer"
                              style={{ display: "flex", alignItems: "center", gap: 6, fontSize: 12, color: "#E8A33D", fontFamily: "'IBM Plex Mono', monospace", textDecoration: "none" }}
                            >
                              <Download size={12} /> Download CIF
                            </a>
                          )}
                        </div>

                        {/* Binding pockets + structural similarity — both on-demand.
                            Buttons stay visible after running (checkmark instead of
                            vanishing) — still clickable to re-run. */}
                        <div style={{ display: "flex", gap: 10, marginTop: 14, flexWrap: "wrap" }}>
                          <button
                            onClick={handleDetectPockets}
                            disabled={pocketsLoading}
                            style={{
                              padding: "6px 12px", borderRadius: 6, border: "1px solid #2A3B5C",
                              background: "#131C2E", color: pocketsResult ? "#3FA796" : "#A8B3C4", fontSize: 11,
                              fontFamily: "'IBM Plex Mono', monospace", cursor: "pointer",
                              display: "flex", alignItems: "center", gap: 6,
                            }}
                          >
                            {pocketsResult && <CheckCircle2 size={12} />}
                            Detect binding pockets
                          </button>
                          <button
                            onClick={handleFindSimilar}
                            disabled={similarityLoading}
                            style={{
                              padding: "6px 12px", borderRadius: 6, border: "1px solid #2A3B5C",
                              background: "#131C2E", color: similarityResult ? "#3FA796" : "#A8B3C4", fontSize: 11,
                              fontFamily: "'IBM Plex Mono', monospace", cursor: "pointer",
                              display: "flex", alignItems: "center", gap: 6,
                            }}
                          >
                            {similarityResult && <CheckCircle2 size={12} />}
                            Find similar structures
                          </button>
                          <button
                            onClick={handleRunDocking}
                            disabled={dockingLoading}
                            style={{
                              padding: "6px 12px", borderRadius: 6, border: "1px solid #2A3B5C",
                              background: "#131C2E", color: dockingResult ? "#3FA796" : "#A8B3C4", fontSize: 11,
                              fontFamily: "'IBM Plex Mono', monospace", cursor: "pointer",
                              display: "flex", alignItems: "center", gap: 6,
                            }}
                          >
                            {dockingResult && <CheckCircle2 size={12} />}
                            Screen inhibitors
                          </button>
                        </div>

                        <div style={{ marginTop: 8, display: "flex", gap: 6, alignItems: "center" }}>
                          <input
                            value={customLigandInput}
                            onChange={(e) => setCustomLigandInput(e.target.value)}
                            placeholder="Custom compounds (comma-separated), e.g. Acetyl-CoA, Coenzyme A"
                            style={{
                              flex: 1, background: "#131C2E", border: "1px solid #2A3B5C", borderRadius: 6,
                              padding: "5px 10px", color: "#E7ECF3", fontSize: 11, outline: "none",
                              fontFamily: "'IBM Plex Mono', monospace",
                            }}
                          />
                        </div>
                        <p style={{ fontSize: 10, color: "#4B5A73", marginTop: 4, marginBottom: 0 }}>
                          Leave blank for automatic selection, or type specific compound names (e.g. ones
                          suggested by KALI AI) to test instead — each must be a real, PubChem-resolvable name.
                        </p>

                        {pocketsLoading && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6", display: "flex", alignItems: "center", gap: 6 }}>
                            <Loader2 size={12} className="spin" color="#E8A33D" /> Running fpocket...
                          </div>
                        )}
                        {pocketsError && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#E85D3D" }}>{pocketsError}</div>
                        )}
                        {pocketsResult && !pocketsResult.available && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6" }}>{pocketsResult.error}</div>
                        )}
                        {pocketsResult?.available && (
                          <div style={{ marginTop: 10 }}>
                            <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                              Predicted binding pockets ({pocketsResult.pockets.length})
                            </div>
                            {pocketsResult.pockets.length === 0 ? (
                              <div style={{ fontSize: 12, color: "#8592A6" }}>No pockets detected.</div>
                            ) : (
                              <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                                {pocketsResult.pockets.map((p) => (
                                  <div key={p.pocket_id} style={{ fontSize: 12, color: "#A8B3C4", display: "flex", gap: 10 }}>
                                    <span style={{ fontFamily: "'IBM Plex Mono', monospace", color: "#E8A33D" }}>#{p.pocket_id}</span>
                                    <span>druggability {p.druggability_score?.toFixed(2) ?? "—"}</span>
                                    <span>volume {p.volume?.toFixed(0) ?? "—"} Å³</span>
                                  </div>
                                ))}
                              </div>
                            )}
                          </div>
                        )}

                        {similarityLoading && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6", display: "flex", alignItems: "center", gap: 6 }}>
                            <Loader2 size={12} className="spin" color="#E8A33D" /> Running Foldseek...
                          </div>
                        )}
                        {similarityError && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#E85D3D" }}>{similarityError}</div>
                        )}
                        {similarityResult && !similarityResult.available && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6" }}>{similarityResult.error}</div>
                        )}
                        {similarityResult?.available && (
                          <div style={{ marginTop: 10 }}>
                            <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                              Similar structures already looked up ({similarityResult.hits.length})
                            </div>
                            {similarityResult.hits.length === 0 ? (
                              <div style={{ fontSize: 12, color: "#8592A6" }}>
                                Nothing to compare against yet — analyze another gene's structure first, then retry.
                              </div>
                            ) : (
                              <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                                {similarityResult.hits.map((h, i) => (
                                  <div key={i} style={{ fontSize: 12, color: "#A8B3C4", display: "flex", gap: 10 }}>
                                    <span style={{ fontFamily: "'IBM Plex Mono', monospace", color: "#E8A33D" }}>
                                      TM {h.tm_score?.toFixed(3) ?? "—"}
                                    </span>
                                    <span>{h.gene_name || h.target_accession}</span>
                                    {h.organism && <span style={{ color: "#6B7A93" }}>({h.organism})</span>}
                                  </div>
                                ))}
                              </div>
                            )}
                          </div>
                        )}

                        {dockingLoading && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6", display: "flex", alignItems: "center", gap: 6 }}>
                            <Loader2 size={12} className="spin" color="#E8A33D" /> Running docking (3 replicate runs per ligand, fetching ligands — may take several minutes)...
                          </div>
                        )}
                        {dockingError && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#E85D3D" }}>{dockingError}</div>
                        )}
                        {dockingResult && !dockingResult.available && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6" }}>{dockingResult.error}</div>
                        )}
                        {dockingResult?.available && (
                          <div style={{ marginTop: 10 }}>
                            <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                              Inhibitor screen — pocket #{dockingResult.pocket_id}
                            </div>
                            <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                              {dockingResult.results.map((r, i) => {
                                const isSuspect = r.available && r.affinity_kcal_mol > 0;
                                return (
                                  <div key={i} style={{ fontSize: 12, color: "#A8B3C4", display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap" }}>
                                    {r.available ? (
                                      <>
                                        <span
                                          style={{
                                            fontFamily: "'IBM Plex Mono', monospace",
                                            color: isSuspect ? "#E85D3D" : "#E8A33D",
                                          }}
                                        >
                                          {r.affinity_kcal_mol?.toFixed(1)} kcal/mol (best)
                                        </span>
                                        {typeof r.mean_affinity_kcal_mol === "number" && (
                                          <span style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace" }}>
                                            mean {r.mean_affinity_kcal_mol.toFixed(1)} ± {r.std_affinity_kcal_mol.toFixed(1)}
                                            {" "}({r.n_replicates}/{r.n_replicates_attempted} runs)
                                          </span>
                                        )}
                                        <span>{r.ligand}</span>
                                        {isSuspect && (
                                          <span style={{ fontSize: 10, color: "#E85D3D" }}>
                                            (positive = clash/no real binding, not a valid result)
                                          </span>
                                        )}
                                        {!isSuspect && r.docked_pose_pdbqt && (
                                          <button
                                            onClick={() => setViewingPose({ ligand: r.ligand, pdbqt: r.docked_pose_pdbqt })}
                                            style={{
                                              padding: "2px 8px", borderRadius: 999, fontSize: 10,
                                              border: "1px solid #2A3B5C", background: "transparent",
                                              color: viewingPose?.ligand === r.ligand ? "#FF3EC9" : "#6B7A93",
                                              cursor: "pointer", fontFamily: "'IBM Plex Mono', monospace",
                                            }}
                                          >
                                            View pose
                                          </button>
                                        )}
                                      </>
                                    ) : (
                                      <span style={{ color: "#4B5A73", whiteSpace: "pre-wrap", fontSize: 11 }}>
                                        {r.ligand} — {r.reason}
                                      </span>
                                    )}
                                  </div>
                                );
                              })}
                            </div>

                            {viewingPose && (
                              <DockingViewer
                                pdbUrl={structureResult.pdb_url}
                                ligandPdbqt={viewingPose.pdbqt}
                                ligandName={viewingPose.ligand}
                              />
                            )}

                            <p style={{ fontSize: 11, color: "#4B5A73", marginTop: 8, marginBottom: 0 }}>
                              More negative = stronger predicted binding. Small, keyword-matched panel —
                              not a comprehensive virtual screen.
                            </p>
                          </div>
                        )}

                        <div style={{ marginTop: 12 }}>
                          <button
                            onClick={handleCheckPdb}
                            disabled={pdbValidationLoading}
                            style={{
                              padding: "6px 12px", borderRadius: 6, border: "1px solid #2A3B5C",
                              background: "#131C2E", color: pdbValidationResult ? "#3FA796" : "#A8B3C4", fontSize: 11,
                              fontFamily: "'IBM Plex Mono', monospace", cursor: "pointer",
                              display: "flex", alignItems: "center", gap: 6,
                            }}
                          >
                            {pdbValidationResult && <CheckCircle2 size={12} />}
                            Check real PDB structures
                          </button>
                        </div>

                        {pdbValidationLoading && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6", display: "flex", alignItems: "center", gap: 6 }}>
                            <Loader2 size={12} className="spin" color="#E8A33D" /> Querying RCSB PDB...
                          </div>
                        )}
                        {pdbValidationError && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#E85D3D" }}>{pdbValidationError}</div>
                        )}
                        {pdbValidationResult && !pdbValidationResult.available && (
                          <div style={{ marginTop: 10, fontSize: 11, color: "#8592A6" }}>{pdbValidationResult.error}</div>
                        )}
                        {pdbValidationResult?.available && (
                          <div style={{ marginTop: 10 }}>
                            <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginBottom: 6 }}>
                              Real crystal structures ({pdbValidationResult.structures.length})
                            </div>
                            {pdbValidationResult.structures.length === 0 ? (
                              <div style={{ fontSize: 12, color: "#8592A6" }}>
                                No experimental structures found for this protein in RCSB PDB — this is
                                normal, not every protein has been crystallized. AlphaFold's prediction
                                is the only structural data available for this gene.
                              </div>
                            ) : (
                              <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
                                {pdbValidationResult.structures.map((s, i) => (
                                  <div key={i} style={{ fontSize: 12, display: "flex", alignItems: "center", gap: 8, flexWrap: "wrap" }}>
                                    <a
                                      href={`https://www.rcsb.org/structure/${s.pdb_id}`}
                                      target="_blank"
                                      rel="noreferrer"
                                      style={{ fontFamily: "'IBM Plex Mono', monospace", color: "#E8A33D", textDecoration: "none" }}
                                    >
                                      {s.pdb_id}
                                    </a>
                                    {s.ligands.length > 0 ? (
                                      s.ligands.map((lig, j) => (
                                        <span
                                          key={j}
                                          style={{
                                            padding: "2px 8px", borderRadius: 999, fontSize: 10,
                                            border: "1px solid #2A3B5C", background: "#131C2E", color: "#8592A6",
                                          }}
                                        >
                                          {lig}
                                        </span>
                                      ))
                                    ) : (
                                      <span style={{ fontSize: 11, color: "#4B5A73" }}>no bound ligand recorded</span>
                                    )}
                                    <button
                                      onClick={() => setViewingPdbId(viewingPdbId === s.pdb_id ? null : s.pdb_id)}
                                      style={{
                                        padding: "2px 8px", borderRadius: 999, fontSize: 10,
                                        border: "1px solid #2A3B5C", background: "transparent",
                                        color: viewingPdbId === s.pdb_id ? "#3FA796" : "#6B7A93",
                                        cursor: "pointer", fontFamily: "'IBM Plex Mono', monospace",
                                      }}
                                    >
                                      {viewingPdbId === s.pdb_id ? "Hide" : "View"}
                                    </button>
                                  </div>
                                ))}
                              </div>
                            )}

                            {viewingPdbId && (
                              <div>
                                <div style={{ fontSize: 10, color: "#6B7A93", fontFamily: "'IBM Plex Mono', monospace", textTransform: "uppercase", letterSpacing: "0.04em", marginTop: 12, marginBottom: 4 }}>
                                  Real structure vs. AlphaFold prediction — scroll up to compare directly
                                </div>
                                <RealStructureViewer pdbId={viewingPdbId} />
                              </div>
                            )}

                            <p style={{ fontSize: 11, color: "#4B5A73", marginTop: 8, marginBottom: 0 }}>
                              Live query against RCSB PDB — real structures, real bound ligand codes.
                              Compare against the docking panel above: a ligand actually crystallized
                              with this protein is a much stronger validation signal than a Vina score alone.
                            </p>
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                </div>
              )}
            </div>
          </>
        )}
      </div>

      {/* Floating chat assistant — grounded in the current analysis */}
      <div style={{ position: "fixed", bottom: 20, right: 20, zIndex: 50 }}>
        {!chatOpen && (
          <button
            onClick={() => { setChatOpen(true); fetchAvailableModels(); }}
            style={{
              width: 52, height: 52, borderRadius: "50%",
              background: "#E8A33D", border: "none", cursor: "pointer",
              display: "flex", alignItems: "center", justifyContent: "center",
              boxShadow: "0 4px 16px rgba(0,0,0,0.4)",
            }}
            title="KALI AI"
          >
            <Info size={22} color="#0B1220" />
          </button>
        )}

        {chatOpen && (
          <div
            style={{
              width: 340, height: 500, background: "#0F1830", border: "1px solid #1B2740",
              borderRadius: 10, display: "flex", flexDirection: "column",
              boxShadow: "0 8px 32px rgba(0,0,0,0.5)", overflow: "hidden",
            }}
          >
            <div style={{ padding: "10px 14px", borderBottom: "1px solid #1B2740", display: "flex", justifyContent: "space-between", alignItems: "center" }}>
              <span style={{ fontSize: 12, fontFamily: "'IBM Plex Mono', monospace", color: "#A8B3C4", fontWeight: 600 }}>
                KALI AI
              </span>
              <button onClick={() => setChatOpen(false)} style={{ background: "none", border: "none", color: "#6B7A93", cursor: "pointer", fontSize: 16 }}>×</button>
            </div>

            <div style={{ padding: "8px 14px", borderBottom: "1px solid #1B2740" }}>
              {modelsLoading && (
                <span style={{ fontSize: 11, color: "#6B7A93" }}>Loading models...</span>
              )}
              {!modelsLoading && availableModels.length > 0 && (
                <select
                  value={selectedModel}
                  onChange={(e) => setSelectedModel(e.target.value)}
                  style={{
                    width: "100%", background: "#131C2E", border: "1px solid #2A3B5C", borderRadius: 6,
                    padding: "4px 8px", color: "#A8B3C4", fontSize: 11, fontFamily: "'IBM Plex Mono', monospace",
                  }}
                >
                  {availableModels.map((m) => (
                    <option key={m.id} value={m.id}>{m.name}</option>
                  ))}
                </select>
              )}
              {!modelsLoading && availableModels.length === 0 && (
                <span style={{ fontSize: 11, color: "#4B5A73" }}>
                  Couldn't load model list \u2014 using OPENROUTER_MODEL default from the backend.
                </span>
              )}
            </div>

            <div style={{ flex: 1, overflowY: "auto", padding: "10px 14px", display: "flex", flexDirection: "column", gap: 10 }}>
              {chatMessages.length === 0 && (
                <div style={{ fontSize: 12, color: "#4B5A73" }}>
                  Ask things like "which genes are HGT?" or "explain this Mantel r value" — answers are grounded in the {genes.length} gene{genes.length !== 1 ? "s" : ""} in your current analysis.
                </div>
              )}
              {chatMessages.map((m, i) => (
                <div
                  key={i}
                  style={{
                    alignSelf: m.role === "user" ? "flex-end" : "flex-start",
                    maxWidth: "85%",
                    padding: "6px 10px",
                    borderRadius: 8,
                    background: m.role === "user" ? "#E8A33D" : "#131C2E",
                    color: m.role === "user" ? "#0B1220" : "#E7ECF3",
                    fontSize: 12,
                    lineHeight: 1.5,
                    whiteSpace: "pre-wrap",
                  }}
                >
                  {m.content}
                </div>
              ))}
              {chatLoading && (
                <div style={{ display: "flex", alignItems: "center", gap: 6, fontSize: 11, color: "#8592A6" }}>
                  <Loader2 size={12} className="spin" color="#E8A33D" /> Thinking...
                </div>
              )}
              {chatError && (
                <div style={{ fontSize: 11, color: "#E85D3D" }}>{chatError}</div>
              )}
            </div>

            <div style={{ padding: 10, borderTop: "1px solid #1B2740", display: "flex", gap: 6 }}>
              <input
                value={chatInput}
                onChange={(e) => setChatInput(e.target.value)}
                onKeyDown={(e) => e.key === "Enter" && handleSendChat()}
                placeholder="Ask a question..."
                style={{
                  flex: 1, background: "#131C2E", border: "1px solid #2A3B5C", borderRadius: 6,
                  padding: "6px 10px", color: "#E7ECF3", fontSize: 12, outline: "none",
                }}
              />
              <button
                onClick={handleSendChat}
                disabled={chatLoading}
                style={{
                  background: "#E8A33D", border: "none", borderRadius: 6, padding: "6px 12px",
                  color: "#0B1220", fontSize: 12, fontWeight: 600, cursor: "pointer",
                }}
              >
                Send
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
