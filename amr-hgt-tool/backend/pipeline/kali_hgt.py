#!/usr/bin/env python3
"""
kali_hgt.py — Phase 3 of the KALI platform

Genomic island and horizontal gene transfer (HGT) detector.

Uses sliding-window k-mer frequency vectors to find genomic regions
whose compositional profile deviates significantly from the genome
background. These are candidate genomic islands, HGT insertions,
pathogenicity islands, or plasmid integrations.

METHOD
------
1. Divide genome into non-overlapping bins of size B bp
2. Compute k-mer frequency vector for each bin
3. Compute genome background vector (mean over all bins)
4. Score each bin by cosine distance from the background
5. Flag bins whose score exceeds mean + Z * std (Z-score threshold)
6. Merge adjacent flagged bins into islands
7. Report coordinates, size, GC deviation, and anomaly score

USAGE
-----
# Basic — detect anomalous regions in one genome
python kali_hgt.py -g genome.fasta -k 4 -b 5000 -o results/ecoli_hgt

# With multiple k values (more sensitive)
python kali_hgt.py -g genome.fasta -k 3 4 5 -b 5000 -o results/ecoli_hgt

# Adjust sensitivity (lower Z = more islands, higher Z = fewer but stronger)
python kali_hgt.py -g genome.fasta -k 4 -b 5000 -z 2.0 -o results/ecoli_hgt

# Multiple genomes at once
python kali_hgt.py -g genomes/*.fna -k 4 -b 5000 -o results/batch_hgt
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.stats import zscore


# ── K-mer helpers ─────────────────────────────────────────────────────────────

BASE_MAP = {'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3}

def kmer_vector(seq: str, k: int, normalize: bool = True) -> np.ndarray:
    """Compute normalised k-mer frequency vector for a sequence."""
    size = 4 ** k
    vec  = np.zeros(size, dtype=np.float64)
    for i in range(len(seq) - k + 1):
        idx, valid = 0, True
        for j in range(k):
            b = BASE_MAP.get(seq[i+j], -1)
            if b == -1: valid = False; break
            idx = idx * 4 + b
        if valid:
            vec[idx] += 1
    if normalize and vec.sum() > 0:
        vec /= vec.sum()
    return vec


def gc_content(seq: str) -> float:
    """GC content as fraction 0–1."""
    seq = seq.upper()
    gc  = seq.count('G') + seq.count('C')
    at  = seq.count('A') + seq.count('T')
    total = gc + at
    return gc / total if total > 0 else 0.0


# ── FASTA reader ──────────────────────────────────────────────────────────────

def read_fasta(path: str) -> dict[str, str]:
    """Return {contig_id: sequence} from a FASTA file."""
    records, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name: records[name] = ''.join(buf)
                name = line[1:].split()[0]
                buf  = []
            elif name:
                buf.append(line.upper().replace(' ', ''))
    if name: records[name] = ''.join(buf)
    return records


# ── Core detection ────────────────────────────────────────────────────────────

def detect_islands(seq: str, k_list: list[int], bin_size: int,
                   z_threshold: float, min_island_bp: int,
                   verbose: bool = False) -> dict:
    """
    Detect genomic islands in a single sequence.

    Returns a dict with:
      - bins     : DataFrame with per-bin stats
      - islands    : DataFrame with merged island regions
      - background_gc : genome-wide GC content
    """
    n = len(seq)
    if n < bin_size * 3:
        if verbose:
            print(f"  Sequence too short ({n} bp) — need at least {bin_size*3} bp, skipping")
        return {"bins":        pd.DataFrame(),
                "islands":       pd.DataFrame(
                    columns=['start','end','size_bp','n_bins',
                             'max_zscore','mean_score','gc','gc_dev']),
                "background_gc": gc_content(seq),
                "n_bins":      0,
                "mean_score":    0.0,
                "std_score":     0.0}

    # Build bin vectors for each k, concatenate
    starts = list(range(0, n - bin_size + 1, bin_size))
    bin_seqs = [seq[s:s+bin_size] for s in starts]
    n_bins   = len(bin_seqs)

    if verbose:
        print(f"  {n:,} bp | {n_bins} bins of {bin_size} bp | k={k_list}")

    # Feature matrix — concat vectors across all k values
    feature_cols = []
    for k in k_list:
        vecs = np.array([kmer_vector(bs, k) for bs in bin_seqs])
        feature_cols.append(vecs)
    X = np.hstack(feature_cols)   # shape: (n_bins, sum(4^k))

    # Genome background = mean vector
    background = X.mean(axis=0)
    bg_norm = np.linalg.norm(background)

    # Score each bin: cosine distance from background
    scores = np.array([
        distance.cosine(X[i], background) if np.linalg.norm(X[i]) > 0 else 0.0
        for i in range(n_bins)
    ])

    # Z-score normalise
    if scores.std() > 0:
        z_scores = (scores - scores.mean()) / scores.std()
    else:
        z_scores = np.zeros(n_bins)

    # GC per bin
    gc_vals = np.array([gc_content(bs) for bs in bin_seqs])
    bg_gc   = gc_content(seq)

    # Build bins DataFrame
    blocks_df = pd.DataFrame({
        'bin_idx':   range(n_bins),
        'start':       starts,
        'end':         [s + bin_size - 1 for s in starts],
        'gc':          gc_vals,
        'gc_dev':      gc_vals - bg_gc,
        'score':       scores,
        'zscore':      z_scores,
        'flagged':     z_scores > z_threshold,
    })

    # Merge adjacent flagged bins into islands
    islands = []
    in_island = False
    island_start = None
    island_blocks = []

    for _, row in blocks_df.iterrows():
        if row['flagged']:
            if not in_island:
                in_island    = True
                island_start = row['start']
                island_blocks = []
            island_blocks.append(row)
        else:
            if in_island:
                island_end = island_blocks[-1]['end']
                size = island_end - island_start + 1
                if size >= min_island_bp:
                    islands.append({
                        'start':      island_start,
                        'end':        island_end,
                        'size_bp':    size,
                        'n_bins':   len(island_blocks),
                        'max_zscore': max(b['zscore'] for b in island_blocks),
                        'mean_score': np.mean([b['score'] for b in island_blocks]),
                        'gc':         np.mean([b['gc'] for b in island_blocks]),
                        'gc_dev':     np.mean([b['gc_dev'] for b in island_blocks]),
                    })
                in_island = False

    # Handle island at end of sequence
    if in_island and island_blocks:
        island_end = island_blocks[-1]['end']
        size = island_end - island_start + 1
        if size >= min_island_bp:
            islands.append({
                'start':      island_start,
                'end':        island_end,
                'size_bp':    size,
                'n_bins':   len(island_blocks),
                'max_zscore': max(b['zscore'] for b in island_blocks),
                'mean_score': np.mean([b['score'] for b in island_blocks]),
                'gc':         np.mean([b['gc'] for b in island_blocks]),
                'gc_dev':     np.mean([b['gc_dev'] for b in island_blocks]),
            })

    islands_df = pd.DataFrame(islands) if islands else pd.DataFrame(
        columns=['start','end','size_bp','n_bins',
                 'max_zscore','mean_score','gc','gc_dev'])

    return {
        'bins':        blocks_df,
        'islands':       islands_df,
        'background_gc': bg_gc,
        'n_bins':      n_bins,
        'mean_score':    float(scores.mean()),
        'std_score':     float(scores.std()),
    }


# ── HTML report ───────────────────────────────────────────────────────────────

def save_html_report(results: dict, out_path: str,
                     genome_name: str, args_info: dict,
                     amr_genes: list = None) -> None:
    """Generate self-contained HTML report — pure SVG, no CDN dependencies."""
    import json as _json

    amr_genes = amr_genes or []

    data = _json.dumps({
        'genome':    genome_name,
        'args':      args_info,
        'contigs':   results,
        'amr_genes': amr_genes,
    }, ensure_ascii=True)

    k_str   = ','.join(str(k) for k in args_info['k_list'])
    meta    = f"k={k_str} | bin={args_info['bin_size']} bp | z-threshold={args_info['z_threshold']} | min-island={args_info['min_island_bp']} bp"

    total_islands = sum(len(v.get('islands',[])) for v in results.values())
    total_bp      = sum(i['size_bp'] for v in results.values() for i in v.get('islands',[]))
    max_z         = max((i['max_zscore'] for v in results.values() for i in v.get('islands',[])), default=0)
    amr_in_isl    = sum(1 for g in amr_genes if g.get('max_island_z',0) > 0)

    cards = [
        ('Contigs',        len(results),             '#1dc9a0'),
        ('Islands Found',  total_islands,             '#f05050'),
        ('Total Island Bp',f'{total_bp:,}',           '#f5a623'),
        ('Max Z-Score',    f'{max_z:.2f}',            '#9b72f8'),
    ]
    if amr_genes:
        cards.append(('AMR Genes',    len(amr_genes), '#f5a623'))
        cards.append(('AMR in Island', amr_in_isl,    '#1dc9a0'))

    cards_html = ''.join(
        f'<div class="card"><div class="cl">{c[0]}</div><div class="cv" style="color:{c[2]}">{c[1]}</div></div>'
        for c in cards
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>KALI HGT &#8212; {genome_name}</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0;}}
body{{font-family:'IBM Plex Sans',Arial,sans-serif;background:#0d0f14;color:#e2e8f8;padding:1.5rem;}}
h1{{font-size:1rem;font-weight:500;color:#1dc9a0;font-family:'IBM Plex Mono',monospace;margin-bottom:.2rem;}}
.meta{{font-size:11px;color:#6b7a9e;font-family:'IBM Plex Mono',monospace;margin-bottom:1.2rem;}}
.cards{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:8px;margin-bottom:1.5rem;}}
.card{{background:#1c2133;border:1px solid #2a3050;border-radius:8px;padding:.6rem .9rem;}}
.cl{{font-size:10px;color:#6b7a9e;margin-bottom:2px;text-transform:uppercase;letter-spacing:.08em;}}
.cv{{font-size:1.2rem;font-weight:600;}}
.sec{{margin-bottom:1.5rem;}}
.sec h2{{font-size:.75rem;color:#6b7a9e;font-family:'IBM Plex Mono',monospace;text-transform:uppercase;
  letter-spacing:.1em;margin-bottom:.75rem;border-bottom:1px solid #2a3050;padding-bottom:4px;}}
svg{{display:block;width:100%;}}
.ctabs{{display:flex;gap:4px;flex-wrap:wrap;margin-bottom:.75rem;}}
.ctab{{padding:4px 10px;border-radius:4px;font-size:10px;font-family:'IBM Plex Mono',monospace;
  border:1px solid #2a3050;background:#1c2133;color:#6b7a9e;cursor:pointer;}}
.ctab.active{{border-color:#5b8df8;color:#5b8df8;background:rgba(91,141,248,.1);}}
table{{width:100%;border-collapse:collapse;font-size:11px;font-family:'IBM Plex Mono',monospace;}}
th,td{{border:.5px solid #2a3050;padding:5px 10px;text-align:left;}}
th{{background:#161a24;color:#6b7a9e;}}
tr:nth-child(even) td{{background:#161a24;}}
.badge{{display:inline-block;padding:2px 7px;border-radius:4px;font-size:10px;font-weight:500;margin-right:2px;}}
.hgt-b{{background:rgba(240,80,80,.2);color:#f05050;}}
.cln-b{{background:rgba(91,141,248,.2);color:#5b8df8;}}
.isl-b{{background:rgba(29,201,160,.2);color:#1dc9a0;}}
.str-b{{background:rgba(240,80,80,.2);color:#f05050;}}
.mod-b{{background:rgba(245,166,35,.2);color:#f5a623;}}
.wek-b{{background:rgba(91,141,248,.2);color:#5b8df8;}}
.leg{{display:flex;gap:14px;margin-bottom:6px;align-items:center;font-size:10px;color:#6b7a9e;flex-wrap:wrap;}}
#tip{{position:fixed;background:#1c2133;border:1px solid #2a3050;border-radius:6px;
  padding:8px 12px;font-size:11px;font-family:'IBM Plex Mono',monospace;pointer-events:none;
  z-index:999;display:none;line-height:1.7;box-shadow:0 4px 16px rgba(0,0,0,.5);max-width:280px;}}
</style>
</head>
<body>
<script id="D" type="application/json">{data}</script>
<h1>KALI HGT detector &#8212; {genome_name}</h1>
<p class="meta">{meta}</p>
<div class="cards">{cards_html}</div>
<div id="tabs"></div>
<div id="view"></div>
<div id="tip"></div>
<script>
var D=JSON.parse(document.getElementById("D").textContent);
var CONTIGS=Object.keys(D.contigs);
var AMR=D.amr_genes||[];
var COL={{HGT:"#f05050",CLONAL:"#5b8df8",UNKNOWN:"#6b7a9e"}};
var tip=document.getElementById("tip");

function showTip(html,e){{tip.innerHTML=html;tip.style.display="block";tip.style.left=(e.clientX+14)+"px";tip.style.top=(e.clientY-10)+"px";}}
function hideTip(){{tip.style.display="none";}}

// Store AMR data by key for safe tooltip access
var _amr={{}};
AMR.forEach(function(g,i){{_amr["a"+i]=g;}});

function amrTip(k,e){{
  var g=_amr[k]; if(!g) return;
  var c=COL[g.gene_class]||"#6b7a9e";
  var inIsl=g.max_island_z>0;
  showTip("<b style=\"color:"+c+"\">"+g.gene_symbol+"</b><br>Type: "+g.gene_class+"<br>Class: "+g.amr_class+"<br>Pos: "+g.start.toLocaleString()+" - "+g.stop.toLocaleString()+"<br>"+(inIsl?"<span style=\"color:#1dc9a0\">&#10003; Z="+g.max_island_z.toFixed(2)+"</span>":"<span style=\"color:#6b7a9e\">Not in island</span>"),e);
}}

// Tabs
if(CONTIGS.length>1){{
  var th="<div class=\"ctabs\">";
  CONTIGS.forEach(function(c,i){{th+="<button class=\"ctab"+(i===0?" active":"")+"\" onclick=\"show(\'"+c+"\')\">"+c+"</button>";}}  );
  th+="</div>";
  document.getElementById("tabs").innerHTML=th;
}}

function show(cid){{
  document.querySelectorAll(".ctab").forEach(function(b){{b.classList.remove("active");}});
  var tb=document.querySelector("[onclick=\"show(\'"+cid+"\')\"");
  if(tb) tb.classList.add("active");

  var cd=D.contigs[cid], bins=cd.bins||[], islands=cd.islands||[], bgGc=cd.background_gc||0;
  var myAmr=AMR.map(function(g,i){{return {{g:g,k:"a"+i}};}}).filter(function(x){{return x.g.contig===cid;}});
  var out="";

  // Legend
  out+="<div class=\"sec\"><h2>Genome Map &#8212; "+cid+"</h2>";
  if(myAmr.length>0){{
    out+="<div class=\"leg\">";
    out+="<span><span style=\"display:inline-block;width:10px;height:8px;background:rgba(240,80,80,0.3);border:1px solid #f05050;margin-right:4px;\"></span>Island</span>";
    out+="<span><span style=\"display:inline-block;width:3px;height:14px;background:#f05050;margin-right:4px;vertical-align:middle;\"></span>HGT gene</span>";
    out+="<span><span style=\"display:inline-block;width:3px;height:14px;background:#5b8df8;margin-right:4px;vertical-align:middle;\"></span>Clonal gene</span>";
    out+="<span><span style=\"display:inline-block;width:8px;height:8px;border-radius:50%;background:#1dc9a0;margin-right:4px;vertical-align:middle;\"></span>In island</span>";
    out+="</div>";
  }}
  out+="<div id=\"gm-"+cid+"\"></div></div>";
  out+="<div class=\"sec\"><h2>Bin Score Profile</h2><div id=\"bs-"+cid+"\"></div></div>";

  // AMR table
  if(myAmr.length>0){{
    out+="<div class=\"sec\"><h2>AMR Genes on this Contig ("+myAmr.length+")</h2>";
    out+="<table><thead><tr><th>Gene</th><th>Drug Class</th><th>Type</th><th>Start</th><th>Stop</th><th>In Island</th><th>Island Z</th></tr></thead><tbody>";
    myAmr.forEach(function(x){{
      var g=x.g, inIsl=g.max_island_z>0;
      var bc=g.gene_class==="HGT"?"hgt-b":g.gene_class==="CLONAL"?"cln-b":"";
      out+="<tr><td><b>"+g.gene_symbol+"</b></td><td>"+g.amr_class+"</td>";
      out+="<td><span class=\"badge "+bc+"\">"+g.gene_class+"</span></td>";
      out+="<td>"+g.start.toLocaleString()+"</td><td>"+g.stop.toLocaleString()+"</td>";
      out+="<td>"+(inIsl?"<span class=\"badge isl-b\">&#10003; Yes</span>":"<span style=\"color:#6b7a9e\">&#8212;</span>")+"</td>";
      out+="<td>"+(inIsl?g.max_island_z.toFixed(2):"&#8212;")+"</td></tr>";
    }});
    out+="</tbody></table></div>";
  }}

  // Islands table
  out+="<div class=\"sec\"><h2>Detected Islands ("+islands.length+")</h2>";
  if(islands.length===0){{
    out+="<div style=\"color:#6b7a9e;padding:1rem;background:#1c2133;border-radius:8px;\">No islands at z="+D.args.z_threshold+"</div>";
  }}else{{
    out+="<table><thead><tr><th>#</th><th>Start</th><th>End</th><th>Size (bp)</th><th>GC%</th><th>GC dev</th><th>Max Z</th><th>AMR genes</th><th>Strength</th></tr></thead><tbody>";
    islands.forEach(function(isl,i){{
      var sc=isl.max_zscore>=4?"str-b":isl.max_zscore>=2.5?"mod-b":"wek-b";
      var lb=isl.max_zscore>=4?"Strong":isl.max_zscore>=2.5?"Moderate":"Weak";
      var ia=myAmr.filter(function(x){{return x.g.start<=isl.end&&x.g.stop>=isl.start;}});
      var ab=ia.length>0?ia.map(function(x){{return "<span class=\"badge "+(x.g.gene_class==="HGT"?"hgt-b":"cln-b")+"\">"+x.g.gene_symbol+"</span>";}}).join(""):"<span style=\"color:#6b7a9e\">&#8212;</span>";
      out+="<tr><td>GI-"+(i+1)+"</td><td>"+isl.start.toLocaleString()+"</td><td>"+isl.end.toLocaleString()+"</td>";
      out+="<td>"+isl.size_bp.toLocaleString()+"</td><td>"+(isl.gc*100).toFixed(1)+"%</td>";
      out+="<td style=\"color:"+(isl.gc_dev>0?"#f5a623":"#5b8df8")+"\">"+(isl.gc_dev>0?"+":"")++(isl.gc_dev*100).toFixed(1)+"%</td>";
      out+="<td>"+isl.max_zscore.toFixed(2)+"</td><td>"+ab+"</td><td><span class=\"badge "+sc+"\">"+lb+"</span></td></tr>";
    }});
    out+="</tbody></table>";
  }}
  out+="</div>";

  document.getElementById("view").innerHTML=out;
  drawMap("gm-"+cid, bins, islands, bgGc, myAmr, D.args);
  drawScore("bs-"+cid, bins, D.args.z_threshold);
}}

function drawMap(id, bins, islands, bgGc, myAmr, args){{
  var el=document.getElementById(id); if(!el||!bins.length) return;
  var W=el.clientWidth||900, ML=50,MR=20,MT=20,MB=35, IH=100, H=IH+MT+MB, IW=W-ML-MR;
  var maxEnd=Math.max.apply(null,bins.map(function(b){{return b.end;}}));
  var gcArr=bins.map(function(b){{return b.gc;}});
  var gcMin=Math.min.apply(null,gcArr), gcMax=Math.max.apply(null,gcArr);
  function xS(v){{return ML+(v/maxEnd)*IW;}}
  function yS(v){{return MT+IH-((v-gcMin)/(gcMax-gcMin||1))*IH;}}
  var s="<svg viewBox=\"0 0 "+W+" "+H+"\" xmlns=\"http://www.w3.org/2000/svg\">";
  s+="<rect x=\""+ML+"\" y=\""+MT+"\" width=\""+IW+"\" height=\""+IH+"\" fill=\"#161a24\" rx=\"3\"/>";
  islands.forEach(function(isl){{
    var x1=xS(isl.start),x2=xS(isl.end);
    s+="<rect x=\""+x1.toFixed(1)+"\" y=\""+MT+"\" width=\""+Math.max(x2-x1,2).toFixed(1)+"\" height=\""+IH+"\" fill=\"rgba(240,80,80,0.2)\" stroke=\"#f05050\" stroke-width=\"1\"/>";
  }});
  var pts=bins.map(function(b){{return xS((b.start+b.end)/2).toFixed(1)+","+yS(b.gc).toFixed(1);}}).join(" ");
  s+="<polyline points=\""+pts+"\" fill=\"none\" stroke=\"#5b8df8\" stroke-width=\"1.5\"/>";
  var bgY=yS(bgGc).toFixed(1);
  s+="<line x1=\""+ML+"\" x2=\""+(ML+IW)+"\" y1=\""+bgY+"\" y2=\""+bgY+"\" stroke=\"#6b7a9e\" stroke-width=\"0.5\" stroke-dasharray=\"4 4\"/>";
  if(myAmr.length>0){{
    var TT=MT+IH-18, TB=MT+IH-2;
    s+="<line x1=\""+ML+"\" x2=\""+(ML+IW)+"\" y1=\""+(TT-2)+"\" y2=\""+(TT-2)+"\" stroke=\"#2a3050\" stroke-width=\"0.5\" stroke-dasharray=\"2 4\"/>";
    s+="<text x=\""+(ML-5)+"\" y=\""+(((TT+TB)/2).toFixed(1))+"\" text-anchor=\"end\" dominant-baseline=\"middle\" font-size=\"9\" fill=\"#6b7a9e\">AMR</text>";
    myAmr.forEach(function(x){{
      var g=x.g, k=x.k;
      var gx=xS((g.start+g.stop)/2).toFixed(1);
      var c=COL[g.gene_class]||"#6b7a9e";
      var inIsl=g.max_island_z>0;
      s+="<line x1=\""+gx+"\" x2=\""+gx+"\" y1=\""+TT+"\" y2=\""+TB+"\" stroke=\""+c+"\" stroke-width=\""+(inIsl?3:2)+"\" opacity=\""+(inIsl?1:0.6)+"\" style=\"cursor:pointer\" onmouseover=\"amrTip(\'"+k+"\',event)\" onmouseout=\"hideTip()\"/>";
      if(inIsl) s+="<circle cx=\""+gx+"\" cy=\""+TT+"\" r=\"4\" fill=\"#1dc9a0\" style=\"cursor:pointer\" onmouseover=\"amrTip(\'"+k+"\',event)\" onmouseout=\"hideTip()\"/>";
    }});
  }}
  var step=Math.pow(10,Math.floor(Math.log10(maxEnd/5)));
  for(var p=0;p<=maxEnd;p+=step){{
    var tx=xS(p).toFixed(1);
    s+="<line x1=\""+tx+"\" x2=\""+tx+"\" y1=\""+(MT+IH)+"\" y2=\""+(MT+IH+4)+"\" stroke=\"#6b7a9e\" stroke-width=\"0.5\"/>";
    s+="<text x=\""+tx+"\" y=\""+(MT+IH+14)+"\" text-anchor=\"middle\" font-size=\"9\" fill=\"#6b7a9e\">"+(p/1000).toFixed(0)+"k</text>";
  }}
  s+="<text x=\""+(ML-5)+"\" y=\""+(MT)+"\" text-anchor=\"end\" dominant-baseline=\"hanging\" font-size=\"9\" fill=\"#6b7a9e\">GC%</text>";
  s+="</svg>";
  el.innerHTML=s;
}}

function drawScore(id, bins, thr){{
  var el=document.getElementById(id); if(!el||!bins.length) return;
  var W=el.clientWidth||900, ML=50,MR=20,MT=10,MB=25, IH=70, H=IH+MT+MB, IW=W-ML-MR;
  var maxZ=Math.max.apply(null,bins.map(function(b){{return Math.abs(b.zscore);}}));
  maxZ=Math.max(maxZ,thr+0.5);
  function xS(i){{return ML+(i/(bins.length-1||1))*IW;}}
  function yS(v){{return MT+IH-(v/maxZ)*IH;}}
  var s="<svg viewBox=\"0 0 "+W+" "+H+"\" xmlns=\"http://www.w3.org/2000/svg\">";
  s+="<line x1=\""+ML+"\" x2=\""+(ML+IW)+"\" y1=\""+yS(thr).toFixed(1)+"\" y2=\""+yS(thr).toFixed(1)+"\" stroke=\"#f05050\" stroke-width=\"0.8\" stroke-dasharray=\"4 3\"/>";
  var bw=Math.max(1,IW/bins.length-0.3);
  bins.forEach(function(b,i){{
    var z=Math.abs(b.zscore), bh=(z/maxZ)*IH;
    s+="<rect x=\""+xS(i).toFixed(1)+"\" y=\""+(MT+IH-bh).toFixed(1)+"\" width=\""+bw.toFixed(1)+"\" height=\""+bh.toFixed(1)+"\" fill=\""+(b.flagged?"#f05050":"#2a3050")+"\"/>";
  }});
  for(var i=0;i<bins.length;i+=Math.floor(bins.length/5)||1){{
    s+="<text x=\""+xS(i).toFixed(1)+"\" y=\""+(MT+IH+14)+"\" text-anchor=\"middle\" font-size=\"9\" fill=\"#6b7a9e\">B"+i+"</text>";
  }}
  s+="</svg>";
  el.innerHTML=s;
}}

show(CONTIGS[0]);
</script>
</body>
</html>"""

    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(html)


def batch_summary(all_results: list[dict], out_base: str,
                  args_info: dict) -> None:
    """
    Cross-genome summary — find genomic positions that are anomalous
    in multiple genomes. Takes the list of per-genome result dicts and
    produces a summary CSV and HTML showing recurrent hotspots.

    Position is expressed as fractional genome position (0.0–1.0) so
    genomes of different lengths can be compared.
    """
    BINS = 20   # divide each genome into 20 positional bins

    # Collect all islands with fractional positions
    records = []
    for entry in all_results:
        genome_name  = entry['genome_name']
        genome_len   = entry['genome_len']
        for contig_id, data in entry['results'].items():
            for isl in data.get('islands', []):
                mid      = (isl['start'] + isl['end']) / 2
                frac_pos = mid / genome_len if genome_len > 0 else 0
                bin_idx  = min(int(frac_pos * BINS), BINS - 1)
                records.append({
                    'genome':     genome_name,
                    'contig':     contig_id,
                    'start':      isl['start'],
                    'end':        isl['end'],
                    'size_bp':    isl['size_bp'],
                    'gc_dev':     isl['gc_dev'],
                    'max_zscore': isl['max_zscore'],
                    'frac_pos':   round(frac_pos, 4),
                    'bin':        bin_idx,
                })

    if not records:
        print("  No islands found across any genome — no summary generated")
        return

    df = pd.DataFrame(records)

    # Count how many genomes have an island in each positional bin
    n_genomes = len(all_results)
    bin_counts = df.groupby('bin').agg(
        n_genomes_with_island=('genome', 'nunique'),
        mean_zscore=('max_zscore', 'mean'),
        mean_gc_dev=('gc_dev', 'mean'),
        total_islands=('genome', 'count'),
    ).reset_index()
    bin_counts['frac_start'] = bin_counts['bin'] / BINS
    bin_counts['frac_end']   = (bin_counts['bin'] + 1) / BINS
    bin_counts['pct_genomes'] = (bin_counts['n_genomes_with_island'] / n_genomes * 100).round(1)

    # Save summary CSV
    csv_path = str(out_base) + '_batch_summary.csv'
    df.to_csv(csv_path, index=False)

    # Save hotspot CSV
    hot_path = str(out_base) + '_hotspots.csv'
    bin_counts.sort_values('pct_genomes', ascending=False).to_csv(hot_path, index=False)

    # Print top hotspots
    print(f"\n── Batch summary ({n_genomes} genomes, {len(records)} total islands) ──")
    print(f"{'Position':>12}  {'Genomes':>8}  {'%':>6}  {'Mean z':>8}  {'Mean GC dev':>12}")
    for _, row in bin_counts.sort_values('pct_genomes', ascending=False).head(10).iterrows():
        pos_label = f"{row['frac_start']*100:.0f}–{row['frac_end']*100:.0f}%"
        print(f"  {pos_label:>10}  {int(row['n_genomes_with_island']):>8}  "
              f"{row['pct_genomes']:>5.1f}%  {row['mean_zscore']:>8.2f}  "
              f"{row['mean_gc_dev']*100:>+11.1f}%")

    # Save HTML summary report
    html_path = str(out_base) + '_batch_summary.html'
    _save_batch_html(df, bin_counts, n_genomes, html_path, args_info)

    print(f"\n  Saved all islands CSV: {csv_path}")
    print(f"  Saved hotspots CSV:    {hot_path}")
    print(f"  Saved summary HTML:    {html_path}")


def _save_batch_html(df: pd.DataFrame, bin_counts: pd.DataFrame,
                     n_genomes: int, out_path: str, args_info: dict) -> None:
    """Generate React HTML batch summary report."""

    payload = json.dumps({
        'n_genomes':  n_genomes,
        'args':       args_info,
        'islands':    df.to_dict(orient='records'),
        'hotspots':   bin_counts.to_dict(orient='records'),
        'genomes':    sorted(df['genome'].unique().tolist()),
    })

    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>KALI HGT Batch Summary</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/react/18.2.0/umd/react.production.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/react-dom/18.2.0/umd/react-dom.production.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/babel-standalone/7.23.2/babel.min.js"></script>
<script id="kali-batch" type="application/json">""" + payload + """</script>
<style>
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:'IBM Plex Sans',Arial,sans-serif;background:#0d0f14;color:#e2e8f8;padding:1.5rem;}
h1{font-size:1rem;font-weight:500;color:#f5a623;font-family:'IBM Plex Mono',monospace;margin-bottom:.25rem;}
.meta{font-size:11px;color:#6b7a9e;font-family:'IBM Plex Mono',monospace;margin-bottom:1rem;}
.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:1.5rem;}
.card{background:#1c2133;border:1px solid #2a3050;border-radius:8px;padding:.6rem .9rem;}
.cl{font-size:10px;color:#6b7a9e;margin-bottom:2px;text-transform:uppercase;letter-spacing:.08em;}
.cv{font-size:1.2rem;font-weight:600;}
h2{font-size:.8rem;font-weight:500;color:#6b7a9e;font-family:'IBM Plex Mono',monospace;
   text-transform:uppercase;letter-spacing:.1em;margin:1.5rem 0 .75rem;}
table{width:100%;border-collapse:collapse;font-size:11px;font-family:'IBM Plex Mono',monospace;}
th,td{border:.5px solid #2a3050;padding:5px 10px;text-align:left;}
th{background:#161a24;color:#6b7a9e;font-weight:500;}
tr:nth-child(even) td{background:#161a24;}
.bar-cell{display:flex;align-items:center;gap:8px;}
.bar-bg{height:8px;background:#2a3050;border-radius:4px;flex:1;}
.bar-fill{height:8px;border-radius:4px;background:#f5a623;}
.badge{display:inline-bin;padding:2px 7px;border-radius:4px;font-size:10px;font-weight:500;}
.hot{background:rgba(240,80,80,.2);color:#f05050;}
.warm{background:rgba(245,166,35,.2);color:#f5a623;}
.cool{background:rgba(91,141,248,.2);color:#5b8df8;}
</style>
</head>
<body>
<div id="root"></div>
<script type="text/babel">
const D = JSON.parse(document.getElementById('kali-batch').textContent);

function App() {
  const totalIslands = D.islands.length;
  const hotspots     = D.hotspots.sort((a,b) => b.pct_genomes - a.pct_genomes);
  const topHot       = hotspots[0] || {};

  return (
    <div>
      <h1>KALI HGT batch summary</h1>
      <p className="meta">
        {D.n_genomes} genomes | k={D.args.k_list.join(',')} |
        bin={D.args.bin_size} bp | z={D.args.z_threshold}
      </p>

      <div className="cards">
        <div className="card">
          <div className="cl">Genomes</div>
          <div className="cv" style={{color:'#1dc9a0'}}>{D.n_genomes}</div>
        </div>
        <div className="card">
          <div className="cl">Total islands</div>
          <div className="cv" style={{color:'#f05050'}}>{totalIslands}</div>
        </div>
        <div className="card">
          <div className="cl">Avg islands/genome</div>
          <div className="cv" style={{color:'#f5a623'}}>
            {(totalIslands/D.n_genomes).toFixed(1)}
          </div>
        </div>
        <div className="card">
          <div className="cl">Top hotspot</div>
          <div className="cv" style={{color:'#9b72f8'}}>
            {topHot.pct_genomes ? topHot.pct_genomes.toFixed(0)+'%' : 'N/A'}
          </div>
        </div>
      </div>

      <h2>Recurrent hotspots by genome position</h2>
      <p style={{fontSize:11,color:'#6b7a9e',marginBottom:12,fontFamily:'monospace'}}>
        Genome divided into 20 positional bins (0–100% of genome length).
        High % = this region is anomalous in many genomes = strong HGT signal.
      </p>
      <HotspotChart hotspots={hotspots} nGenomes={D.n_genomes} />

      <h2>Hotspot table</h2>
      <table>
        <thead>
          <tr>
            <th>Position</th>
            <th>Genomes affected</th>
            <th>% of all genomes</th>
            <th>Mean z-score</th>
            <th>Mean GC dev</th>
            <th>Signal</th>
          </tr>
        </thead>
        <tbody>
          {hotspots.filter(h => h.n_genomes_with_island > 0).map((h,i) => {
            const pct      = h.pct_genomes;
            const strength = pct >= 50 ? 'hot' : pct >= 25 ? 'warm' : 'cool';
            const label    = pct >= 50 ? 'Recurrent' : pct >= 25 ? 'Moderate' : 'Sporadic';
            return (
              <tr key={i}>
                <td>{(h.frac_start*100).toFixed(0)}–{(h.frac_end*100).toFixed(0)}%</td>
                <td>{h.n_genomes_with_island} / {D.n_genomes}</td>
                <td>
                  <div className="bar-cell">
                    <div className="bar-bg">
                      <div className="bar-fill" style={{width:pct+'%'}} />
                    </div>
                    <span>{pct.toFixed(1)}%</span>
                  </div>
                </td>
                <td>{h.mean_zscore.toFixed(2)}</td>
                <td style={{color: h.mean_gc_dev < 0 ? '#5b8df8' : '#f5a623'}}>
                  {h.mean_gc_dev > 0 ? '+' : ''}{(h.mean_gc_dev*100).toFixed(1)}%
                </td>
                <td><span className={'badge ' + strength}>{label}</span></td>
              </tr>
            );
          })}
        </tbody>
      </table>

      <h2>Per-genome island count</h2>
      <table>
        <thead>
          <tr><th>Genome</th><th>Islands</th><th>Positions (% of genome)</th></tr>
        </thead>
        <tbody>
          {D.genomes.map(g => {
            const gi = D.islands.filter(r => r.genome === g);
            return (
              <tr key={g}>
                <td style={{fontFamily:'monospace',fontSize:10}}>{g}</td>
                <td>{gi.length}</td>
                <td style={{fontSize:10,color:'#6b7a9e'}}>
                  {gi.map(r => (r.frac_pos*100).toFixed(0)+'%').join(', ') || '—'}
                </td>
              </tr>
            );
          })}
        </tbody>
      </table>
    </div>
  );
}

function HotspotChart({hotspots, nGenomes}) {
  const ref = React.useRef(null);
  React.useEffect(() => {
    if (!ref.current) return;
    const W = ref.current.clientWidth || 700;
    const H = 140;
    const mL = 50, mR = 20, mT = 10, mB = 30;
    const iW = W - mL - mR, iH = H - mT - mB;

    d3.select(ref.current).selectAll('*').remove();
    const svg = d3.select(ref.current).append('svg').attr('width',W).attr('height',H);
    const g   = svg.append('g').attr('transform',`translate(${mL},${mT})`);

    const xScale = d3.scaleBand()
      .domain(hotspots.map(h => h.bin))
      .range([0, iW]).padding(0.1);
    const yScale = d3.scaleLinear([0, 100], [iH, 0]);

    // 50% reference line
    g.append('line')
     .attr('x1',0).attr('x2',iW)
     .attr('y1',yScale(50)).attr('y2',yScale(50))
     .attr('stroke','#f05050').attr('stroke-width',.8)
     .attr('stroke-dasharray','4 3');
    g.append('text').attr('x',iW+4).attr('y',yScale(50)+4)
     .style('font-size','9px').style('fill','#f05050').text('50%');

    // Bars
    hotspots.forEach(h => {
      const pct = h.pct_genomes;
      const col = pct >= 50 ? '#f05050' : pct >= 25 ? '#f5a623' : '#2a3050';
      g.append('rect')
       .attr('x', xScale(h.bin))
       .attr('y', yScale(pct))
       .attr('width', xScale.bandwidth())
       .attr('height', iH - yScale(pct))
       .attr('fill', col)
       .attr('rx', 2);
    });

    // X axis — show 0%, 25%, 50%, 75%, 100% labels
    const ticks = [0, 5, 10, 15, 19];
    ticks.forEach(bin => {
      const x = xScale(bin) + xScale.bandwidth()/2;
      if (x === undefined) return;
      g.append('text').attr('x', x || 0)
       .attr('y', iH + 16)
       .attr('text-anchor','middle')
       .style('font-size','9px').style('fill','#6b7a9e')
       .text(Math.round(bin/19*100) + '%');
    });

    g.append('text').attr('x',-5).attr('y',0)
     .attr('text-anchor','end').attr('dominant-baseline','hanging')
     .style('font-size','9px').style('fill','#6b7a9e').text('%');

    // Y axis
    g.append('g').call(d3.axisLeft(yScale).ticks(4)
      .tickFormat(d => d+'%'))
     .selectAll('text').style('fill','#6b7a9e').style('font-size','9px');

    svg.append('text').attr('x', mL + iW/2).attr('y', H-2)
     .attr('text-anchor','middle')
     .style('font-size','9px').style('fill','#6b7a9e')
     .text('Relative genome position');

  }, [hotspots]);

  return <div ref={ref} style={{width:'100%'}} />;
}

ReactDOM.createRoot(document.getElementById('root')).render(<App/>);
</script>
</body>
</html>"""

    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(html)


def save_csv(results: dict, out_path: str, genome_name: str) -> None:
    """Save islands table as CSV."""
    rows = []
    for contig_id, data in results.items():
        for i, isl in enumerate(data.get('islands', [])):
            rows.append({
                'genome':      genome_name,
                'contig':      contig_id,
                'island_id':   f'GI-{i+1}',
                'start':       isl['start'],
                'end':         isl['end'],
                'size_bp':     isl['size_bp'],
                'n_bins':    isl['n_bins'],
                'gc':          round(isl['gc'], 4),
                'gc_dev':      round(isl['gc_dev'], 4),
                'max_zscore':  round(isl['max_zscore'], 4),
                'mean_score':  round(isl['mean_score'], 6),
            })
    df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=['genome','contig','island_id','start','end',
                 'size_bp','n_bins','gc','gc_dev','max_zscore','mean_score'])
    df.to_csv(out_path, index=False)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="KALI Phase 3 — Genomic island and HGT detector"
    )
    p.add_argument('-g', '--genomes', nargs='+', required=True,
                   help='FASTA file(s) to analyse')
    p.add_argument('-k', '--k', nargs='+', type=int, default=[4],
                   help='k-mer size(s) (default: 4)')
    p.add_argument('-b', '--bin-size', type=int, default=5000,
                   help='Bin size in bp (default: 5000)')
    p.add_argument('-z', '--z-threshold', type=float, default=2.5,
                   help='Z-score threshold for flagging (default: 2.5)')
    p.add_argument('--min-island', type=int, default=5000,
                   help='Minimum island size in bp (default: 5000)')
    p.add_argument('-o', '--output', required=True,
                   help='Output base path')
    p.add_argument('--amr', type=str, default='',
                   help='Path to AMRFinderPlus TSV — adds AMR gene track to HTML reports')
    p.add_argument('--summary', action='store_true',
                   help='Generate batch summary showing recurrent hotspots across all genomes')
    p.add_argument('-v', '--verbose', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    out_base = Path(args.output)
    out_base.parent.mkdir(parents=True, exist_ok=True)

    args_info = {
        'k_list':       args.k,
        'bin_size':   args.bin_size,
        'z_threshold':  args.z_threshold,
        'min_island_bp': args.min_island,
    }

    # Load AMR data if provided
    amr_lookup = {}
    if hasattr(args, 'amr') and args.amr and Path(args.amr).exists():
        import csv as _csv
        print(f"  Loading AMR data from: {args.amr}")
        with open(args.amr, newline='') as f:
            reader = _csv.DictReader(f, delimiter='\t')
            for row in reader:
                genome = row.get('genome','').strip()
                try:
                    sym = row.get('Gene symbol','').strip()
                    clonal_p = ('gyrA_','parC_','parE_','uhpT_','ptsI_','cyaA_','cirA_','ompC_','ompK')
                    hgt_p    = ('bla','aadA','aac(','aph(','sul','tet(','dfrA','cat','mph','sat','mrx','qac','flo','arr','rmtF','fosA','qep','ble')
                    if any(sym.startswith(p) for p in clonal_p) or ('_' in sym and len(sym.split('_'))>1 and sym.split('_')[1][0].isupper()):
                        gc = 'CLONAL'
                    elif any(sym.startswith(p) for p in hgt_p):
                        gc = 'HGT'
                    else:
                        gc = 'UNKNOWN'
                    gene = {
                        'contig':       row.get('Contig id','').strip(),
                        'start':        int(row.get('Start', 0)),
                        'stop':         int(row.get('Stop',  0)),
                        'gene_symbol':  sym,
                        'amr_class':    row.get('Class','').strip(),
                        'gene_class':   gc,
                        'max_island_z': 0.0,
                    }
                    amr_lookup.setdefault(genome, []).append(gene)
                except (ValueError, KeyError):
                    continue
        total = sum(len(v) for v in amr_lookup.values())
        print(f"  Loaded {total} AMR genes across {len(amr_lookup)} genomes")

    # Expand glob patterns
    from glob import glob
    genome_paths = []
    for g in args.genomes:
        expanded = glob(g)
        genome_paths.extend(expanded if expanded else [g])

    print(f"\n── KALI HGT Detector ──────────────────────────────────")
    print(f"   Genomes:    {len(genome_paths)}")
    print(f"   k:          {args.k}")
    print(f"   Bin size: {args.bin_size:,} bp")
    print(f"   Z-threshold:{args.z_threshold}")
    print(f"   Min island: {args.min_island:,} bp")
    print()

    all_island_count = 0
    all_results = []   # for batch summary

    for gpath in genome_paths:
        if not Path(gpath).exists():
            print(f"  Not found: {gpath}")
            continue

        genome_name = Path(gpath).stem
        print(f"Processing: {genome_name}")

        contigs = read_fasta(gpath)
        results = {}

        for contig_id, seq in contigs.items():
            if args.verbose:
                print(f"  Contig: {contig_id} ({len(seq):,} bp)")
            res = detect_islands(
                seq, args.k, args.bin_size,
                args.z_threshold, args.min_island,
                verbose=args.verbose
            )
            # Skip contigs too short to analyse
            if res.get('n_bins', 0) == 0:
                print(f"  {contig_id}: skipped (too short — need {args.bin_size*3:,} bp minimum)")
                continue
            # Convert DataFrames to dicts for JSON serialisation
            results[contig_id] = {
                'background_gc': res['background_gc'],
                'n_bins':      res['n_bins'],
                'mean_score':    res['mean_score'],
                'std_score':     res['std_score'],
                'bins':  res['bins'].to_dict(orient='records'),
                'islands': res['islands'].to_dict(orient='records'),
            }
            n_isl = len(res['islands'])
            all_island_count += n_isl
            print(f"  {contig_id}: {n_isl} island(s) detected  "
                  f"(GC={res['background_gc']*100:.1f}%  "
                  f"bins={res['n_bins']})")

        # Save outputs
        suffix = '_' + genome_name if len(genome_paths) > 1 else ''
        html_path = str(out_base) + suffix + '_hgt.html'
        csv_path  = str(out_base) + suffix + '_islands.csv'

        # Annotate AMR genes with island Z-scores
        amr_genes = []
        for gene in amr_lookup.get(genome_name, []):
            g = dict(gene)
            for isl in results.get(g['contig'], {}).get('islands', []):
                if isl['start'] <= g['stop'] and isl['end'] >= g['start']:
                    g['max_island_z'] = round(max(g['max_island_z'], isl['max_zscore']), 3)
            amr_genes.append(g)
        save_html_report(results, html_path, genome_name, args_info, amr_genes=amr_genes)
        save_csv(results, csv_path, genome_name)

        print(f"  Saved HTML: {html_path}")
        print(f"  Saved CSV:  {csv_path}")
        print()

        # Collect for batch summary
        genome_len = sum(len(s) for s in read_fasta(gpath).values())
        all_results.append({
            'genome_name': genome_name,
            'genome_len':  genome_len,
            'results':     results,
        })

    print(f"Total islands detected: {all_island_count}")

    # Batch summary across all genomes
    if args.summary and all_results:
        batch_summary(all_results, str(out_base), args_info)

    print(f"Done.\n")


if __name__ == '__main__':
    main()
