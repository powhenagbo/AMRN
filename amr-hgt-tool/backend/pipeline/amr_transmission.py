"""
AMR Transmission Inference Module
===================================
Novel approach: Uses k-mer phylogenetic distances (Pykali) to statistically
distinguish CLONAL spread vs HORIZONTAL GENE TRANSFER (HGT) of AMR genes
— entirely alignment-free.

Core Logic:
-----------
If an AMR gene spread CLONALLY → isolates sharing the gene should be
phylogenetically CLOSE (low k-mer distance). High Mantel r = clonal signal.
If an AMR gene spread via HGT → isolates sharing the gene are
phylogenetically DISTANT. Low Mantel r = HGT signal.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings("ignore")


def _upper_triangle_values(matrix):
    arr = np.asarray(matrix, dtype=float)
    idx = np.triu_indices_from(arr, k=1)
    return arr[idx]


def mantel_test(matrix_a, matrix_b, permutations=999, seed=42):
    """
    Lightweight Mantel test using Pearson correlation on the upper triangles
    of two symmetric distance matrices. The permutation step shuffles sample
    labels on matrix_b, preserving its distance structure while breaking the
    correspondence to matrix_a.

    Returns (r, p_value, permutations_run), matching the values used by the
    previous scikit-bio call closely enough for the surrounding classification
    logic to remain unchanged.
    """
    a = np.asarray(matrix_a, dtype=float)
    b = np.asarray(matrix_b, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError("Mantel inputs must be square matrices of the same shape.")

    av = _upper_triangle_values(a)
    bv = _upper_triangle_values(b)
    if av.size < 2:
        raise ValueError("Mantel test requires at least three samples.")

    r_obs, _ = pearsonr(av, bv)
    if np.isnan(r_obs):
        raise ValueError("Mantel correlation is undefined for constant distances.")

    rng = np.random.default_rng(seed)
    extreme = 0
    n = a.shape[0]
    for _ in range(int(permutations)):
        perm = rng.permutation(n)
        bp = b[np.ix_(perm, perm)]
        rp, _ = pearsonr(av, _upper_triangle_values(bp))
        if not np.isnan(rp) and abs(rp) >= abs(r_obs):
            extreme += 1

    p_value = (extreme + 1) / (int(permutations) + 1)
    return float(r_obs), float(p_value), int(permutations)


# ──────────────────────────────────────────────────────────
# Mantel Test: per-gene clonal vs HGT classification
# ──────────────────────────────────────────────────────────

def amr_sharing_matrix(amr_matrix: pd.DataFrame) -> pd.DataFrame:
    """Pairwise AMR Jaccard dissimilarity matrix."""
    samples = amr_matrix.index.tolist()
    n = len(samples)
    sharing = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        a = amr_matrix.iloc[i].values.astype(bool)
        b = amr_matrix.iloc[j].values.astype(bool)
        intersection = np.sum(a & b)
        union = np.sum(a | b)
        jaccard_sim = intersection / union if union > 0 else 0.0
        sharing[i][j] = 1 - jaccard_sim
        sharing[j][i] = 1 - jaccard_sim
    return pd.DataFrame(sharing, index=samples, columns=samples)


def per_gene_mantel(
    dist_df: pd.DataFrame,
    amr_matrix: pd.DataFrame,
    permutations: int = 999,
    clonal_r: float = 0.4,
    hgt_r: float = 0.2
) -> pd.DataFrame:
    """
    Run Mantel test for each AMR gene independently.
    permutations, clonal_r, hgt_r are all configurable from CLI.
    """
    common = sorted(set(dist_df.index) & set(amr_matrix.index))
    if not common:
        raise ValueError("No overlapping samples between distance matrix and AMR matrix.")

    dist_aligned = dist_df.loc[common, common]
    amr_aligned  = amr_matrix.loc[common]
    n = len(common)

    results = []
    total_genes = amr_aligned.shape[1]

    for idx, gene in enumerate(amr_aligned.columns):
        presence = amr_aligned[gene].values.astype(int)
        n_carriers = int(presence.sum())

        if n_carriers < 2 or n_carriers == n:
            continue

        gene_dist = np.zeros((n, n))
        for i, j in combinations(range(n), 2):
            val = 0 if (presence[i] == 1 and presence[j] == 1) else 1
            gene_dist[i][j] = val
            gene_dist[j][i] = val

        try:
            r, p, _ = mantel_test(
                dist_aligned.values,
                gene_dist,
                permutations=permutations,
                seed=42,
            )
        except Exception as e:
            print(f"  [WARN] Mantel failed for {gene}: {e}")
            continue

        if r >= clonal_r and p <= 0.05:
            classification = "CLONAL"
        elif r < hgt_r or p > 0.05:
            classification = "HGT"
        else:
            classification = "AMBIGUOUS"

        results.append({
            "gene": gene,
            "mantel_r": round(r, 4),
            "p_value": round(p, 4),
            "n_carriers": n_carriers,
            "classification": classification
        })

        if (idx + 1) % 10 == 0:
            print(f"  Processed {idx+1}/{total_genes} genes...")

    return pd.DataFrame(results).sort_values("mantel_r", ascending=False)


def plot_transmission_landscape(results_df: pd.DataFrame, outdir: str):
    """Mantel r vs -log10(p) per gene — core novel figure."""
    if results_df.empty:
        print("[WARN] No results to plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    color_map = {"CLONAL": "#e74c3c", "HGT": "#3498db", "AMBIGUOUS": "#95a5a6"}

    ax = axes[0]
    colors = results_df["classification"].map(color_map)
    sizes  = results_df["n_carriers"] * 30
    p_clipped = results_df["p_value"].clip(lower=1e-10)

    ax.scatter(
        results_df["mantel_r"],
        -np.log10(p_clipped),
        c=colors, s=sizes, alpha=0.75,
        edgecolors="white", linewidth=0.5
    )

    ax.axvline(x=0.4, color="#e74c3c", linestyle="--", alpha=0.5, label="Clonal r=0.4")
    ax.axvline(x=0.2, color="#3498db", linestyle="--", alpha=0.5, label="HGT r=0.2")
    ax.axhline(y=-np.log10(0.05), color="gray", linestyle=":", alpha=0.5, label="p=0.05")

    for _, row in pd.concat([
        results_df.nlargest(5, "mantel_r"),
        results_df.nsmallest(5, "mantel_r")
    ]).iterrows():
        ax.annotate(
            row["gene"][:20],
            (row["mantel_r"], -np.log10(max(row["p_value"], 1e-10))),
            fontsize=7, alpha=0.8,
            xytext=(5, 5), textcoords="offset points"
        )

    from matplotlib.patches import Patch
    legend_patches = [Patch(color=c, label=k) for k, c in color_map.items()]
    ax.legend(handles=legend_patches, fontsize=8)
    ax.set_xlabel("Mantel r (phylogenetic-AMR correlation)", fontsize=11)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=11)
    ax.set_title("AMR Transmission Mode per Gene\n(Alignment-Free Mantel Test)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.2)

    ax2 = axes[1]
    counts = results_df["classification"].value_counts()
    bar_colors = [color_map.get(c, "gray") for c in counts.index]
    bars = ax2.bar(counts.index, counts.values,
                   color=bar_colors, edgecolor="white", linewidth=1.5)
    for bar, count in zip(bars, counts.values):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.2,
            str(count), ha="center", va="bottom",
            fontsize=12, fontweight="bold"
        )
    ax2.set_title("AMR Gene Transmission Classification",
                  fontsize=12, fontweight="bold")
    ax2.set_xlabel("Transmission Mode", fontsize=11)
    ax2.set_ylabel("Number of AMR Genes", fontsize=11)
    ax2.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    out_path = f"{outdir}/amr_transmission_landscape.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[Output] Transmission landscape → {out_path}")
    plt.close()


def plot_distance_vs_sharing(dist_df: pd.DataFrame, amr_matrix: pd.DataFrame, outdir: str):
    """Pairwise k-mer distance vs AMR gene sharing — clonal signal check."""
    common = sorted(set(dist_df.index) & set(amr_matrix.index))
    dist_aligned = dist_df.loc[common, common]
    sharing_df   = amr_sharing_matrix(amr_matrix.loc[common])

    dists, sharings = [], []
    for i, j in combinations(range(len(common)), 2):
        dists.append(dist_aligned.iloc[i, j])
        sharings.append(1 - sharing_df.iloc[i, j])

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(dists, sharings, alpha=0.3, s=15, color="#2c3e50")

    z = np.polyfit(dists, sharings, 1)
    p_fit = np.poly1d(z)
    x_line = np.linspace(min(dists), max(dists), 100)
    ax.plot(x_line, p_fit(x_line), "r--", linewidth=2, label=f"Trend (slope={z[0]:.3f})")

    r, pval = pearsonr(dists, sharings)
    ax.set_title(
        f"Phylogenetic Distance vs AMR Gene Sharing\nr = {r:.3f}, p = {pval:.4f}",
        fontsize=12, fontweight="bold"
    )
    ax.set_xlabel("Pykali k-mer Distance (Bray-Curtis)", fontsize=11)
    ax.set_ylabel("AMR Gene Sharing (Jaccard Similarity)", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = f"{outdir}/distance_vs_amr_sharing.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[Output] Distance vs sharing → {out_path}")
    plt.close()


def run_transmission_analysis(
    dist_df: pd.DataFrame,
    amr_matrix: pd.DataFrame,
    outdir: str,
    permutations: int = 999,
    clonal_r: float = 0.4,
    hgt_r: float = 0.2
):
    """Full transmission inference pipeline."""
    import os
    os.makedirs(outdir, exist_ok=True)

    print(f"[Module 1] Running per-gene Mantel tests ({permutations} permutations)...")
    results = per_gene_mantel(
        dist_df, amr_matrix,
        permutations=permutations,
        clonal_r=clonal_r,
        hgt_r=hgt_r
    )

    if results.empty:
        print("[WARN] No informative genes found.")
        return results

    results.to_csv(f"{outdir}/transmission_classification.csv", index=False)
    print(f"[Output] Classification table → {outdir}/transmission_classification.csv")

    print("\n[Summary]")
    print(results["classification"].value_counts().to_string())
    print("\nTop CLONAL genes:")
    clonal = results[results["classification"] == "CLONAL"]
    if not clonal.empty:
        print(clonal.head(5)[["gene", "mantel_r", "p_value"]].to_string(index=False))
    print("\nTop HGT genes:")
    hgt = results[results["classification"] == "HGT"]
    if not hgt.empty:
        print(hgt.head(5)[["gene", "mantel_r", "p_value"]].to_string(index=False))

    plot_transmission_landscape(results, outdir)
    plot_distance_vs_sharing(dist_df, amr_matrix, outdir)

    return results
