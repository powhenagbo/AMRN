"""
AMR Co-occurrence Network Module
==================================
Novel approach: Builds a network of AMR genes that co-occur across isolates,
then tests whether each co-occurring gene PAIR is phylogenetically CONSERVED
(same clade = co-transferred as a unit) or SCATTERED (independent acquisition).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations
from scipy.stats import fisher_exact
import warnings
warnings.filterwarnings("ignore")


# ──────────────────────────────────────────────────────────
# Gene-gene co-occurrence statistics
# ──────────────────────────────────────────────────────────

def phi_coefficient(a: np.ndarray, b: np.ndarray) -> tuple:
    """Phi coefficient + Fisher's exact p-value between two binary vectors."""
    n11 = np.sum(a & b)
    n10 = np.sum(a & ~b)
    n01 = np.sum(~a & b)
    n00 = np.sum(~a & ~b)
    table = np.array([[n11, n10], [n01, n00]])
    denom = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    phi = (n11 * n00 - n10 * n01) / denom if denom > 0 else 0.0
    _, p = fisher_exact(table, alternative="greater")
    return phi, p, int(n11)


def build_cooccurrence_network(
    amr_matrix: pd.DataFrame,
    phi_threshold: float = 0.3,
    p_threshold: float = 0.05,
    min_cooccurrences: int = 2      # FIX: now actually applied
) -> nx.Graph:
    """Build gene co-occurrence network. Applies all three filters."""
    genes = amr_matrix.columns.tolist()
    G = nx.Graph()

    for gene in genes:
        freq = amr_matrix[gene].sum() / len(amr_matrix)
        G.add_node(gene, frequency=round(freq, 3), n_carriers=int(amr_matrix[gene].sum()))

    for g1, g2 in combinations(genes, 2):
        a = amr_matrix[g1].values.astype(bool)
        b = amr_matrix[g2].values.astype(bool)
        phi, p, n_both = phi_coefficient(a, b)

        # FIX: min_cooccurrences now actually filters edges
        if phi >= phi_threshold and p <= p_threshold and n_both >= min_cooccurrences:
            G.add_edge(g1, g2, phi=round(phi, 4), p_value=round(p, 4), n_cooccurrences=n_both)

    return G


# ──────────────────────────────────────────────────────────
# Phylogenetic scatter scoring
# ──────────────────────────────────────────────────────────

def phylogenetic_scatter_score(
    g1: str, g2: str,
    amr_matrix: pd.DataFrame,
    dist_df: pd.DataFrame,
    mean_all_dist: float        # FIX: pre-computed baseline passed in
) -> float:
    """
    Scatter score for a gene pair using pre-computed baseline distance.
    Score ~0 = co-carriers clustered (CONSERVED).
    Score ~1 = co-carriers dispersed (SCATTERED).
    """
    common = sorted(set(dist_df.index) & set(amr_matrix.index))
    amr    = amr_matrix.loc[common]
    dist   = dist_df.loc[common, common]

    both = amr[(amr[g1] == 1) & (amr[g2] == 1)].index.tolist()
    if len(both) < 2:
        return None

    co_dists = [dist.loc[i, j] for i, j in combinations(both, 2)]
    mean_co_dist = np.mean(co_dists)
    scatter = mean_co_dist / mean_all_dist if mean_all_dist > 0 else 0.0
    return round(scatter, 4)


def classify_cooccurrences(
    G: nx.Graph,
    amr_matrix: pd.DataFrame,
    dist_df: pd.DataFrame,
    scatter_threshold: float = 0.7
) -> pd.DataFrame:
    """Classify all edges as CONSERVED or SCATTERED."""

    # FIX: compute baseline mean distance ONCE here, not inside the loop
    common = sorted(set(dist_df.index) & set(amr_matrix.index))
    dist   = dist_df.loc[common, common]
    all_dists = [dist.iloc[i, j] for i, j in combinations(range(len(common)), 2)]
    mean_all_dist = np.mean(all_dists) if all_dists else 1.0
    print(f"  Baseline mean pairwise distance: {mean_all_dist:.4f}")

    records = []
    for g1, g2, data in G.edges(data=True):
        score = phylogenetic_scatter_score(g1, g2, amr_matrix, dist_df, mean_all_dist)

        if score is None:
            classification = "INSUFFICIENT_DATA"
        elif score < scatter_threshold:
            classification = "CONSERVED"
        else:
            classification = "SCATTERED"

        records.append({
            "gene_1":         g1,
            "gene_2":         g2,
            "phi":            data["phi"],
            "p_value":        data["p_value"],
            "n_cooccurrences":data["n_cooccurrences"],
            "scatter_score":  score,
            "classification": classification
        })

        G[g1][g2]["scatter_score"]  = score
        G[g1][g2]["classification"] = classification

    return pd.DataFrame(records).sort_values("scatter_score")


# ──────────────────────────────────────────────────────────
# Visualization
# ──────────────────────────────────────────────────────────

def plot_cooccurrence_network(G: nx.Graph, drug_class_map: dict, outdir: str):
    """Draw AMR co-occurrence network."""
    if G.number_of_edges() == 0:
        print("[WARN] No significant co-occurrences found. Try lowering phi_threshold.")
        return

    DRUG_CLASS_COLORS = {
        "aminoglycoside": "#e74c3c", "beta-lactam":    "#3498db",
        "fluoroquinolone":"#2ecc71", "tetracycline":   "#f39c12",
        "macrolide":      "#9b59b6", "sulfonamide":    "#1abc9c",
        "glycopeptide":   "#e67e22", "colistin":       "#c0392b",
    }

    def get_color(gene):
        dc = drug_class_map.get(gene, "").lower()
        for key, color in DRUG_CLASS_COLORS.items():
            if key in dc:
                return color
        return "#bdc3c7"

    fig, ax = plt.subplots(figsize=(14, 10))
    pos = nx.spring_layout(G, k=2.5, seed=42, weight="phi")

    node_colors = [get_color(n) for n in G.nodes()]
    node_sizes  = [G.nodes[n].get("n_carriers", 1) * 80 + 100 for n in G.nodes()]
    edge_colors = []
    edge_widths = []
    for u, v, data in G.edges(data=True):
        cls = data.get("classification", "INSUFFICIENT_DATA")
        edge_colors.append(
            "#e74c3c" if cls == "CONSERVED" else
            "#3498db" if cls == "SCATTERED" else "#bdc3c7"
        )
        edge_widths.append(max(1, data.get("phi", 0.3) * 6))

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, width=edge_widths, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=node_sizes, alpha=0.9)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7)

    from matplotlib.patches import Patch
    drug_patches = [Patch(color=c, label=k) for k, c in DRUG_CLASS_COLORS.items()]
    edge_patches = [
        Patch(color="#e74c3c", label="CONSERVED (co-transfer)"),
        Patch(color="#3498db", label="SCATTERED (independent)")
    ]
    ax.legend(handles=drug_patches + edge_patches, loc="lower left",
              fontsize=8, title="Node=Drug Class | Edge=Transfer Mode", ncol=2)
    ax.set_title("AMR Gene Co-occurrence Network\nEdge color = phylogenetic scatter classification",
                 fontsize=13, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    out_path = f"{outdir}/amr_cooccurrence_network.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[Output] Co-occurrence network → {out_path}")
    plt.close()


def plot_scatter_score_distribution(results_df: pd.DataFrame, outdir: str):
    if results_df.empty or results_df["scatter_score"].isna().all():
        print("[WARN] No scatter scores to plot.")
        return

    valid = results_df.dropna(subset=["scatter_score"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.hist(valid["scatter_score"], bins=20, color="#2c3e50", edgecolor="white", alpha=0.8)
    ax.axvline(x=0.7, color="#e74c3c", linestyle="--", linewidth=2, label="Scatter threshold (0.7)")
    ax.set_xlabel("Phylogenetic Scatter Score", fontsize=11)
    ax.set_ylabel("Number of Gene Pairs", fontsize=11)
    ax.set_title("Distribution of AMR Co-occurrence\nPhylogenetic Scatter Scores",
                 fontsize=12, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.2)

    ax2 = axes[1]
    color_map = {"CONSERVED": "#e74c3c", "SCATTERED": "#3498db", "INSUFFICIENT_DATA": "#95a5a6"}
    colors = valid["classification"].map(color_map)
    ax2.scatter(valid["phi"], valid["scatter_score"],
                c=colors, alpha=0.7, s=60, edgecolors="white")
    ax2.axhline(y=0.7, color="gray", linestyle=":", alpha=0.7)
    ax2.set_xlabel("Phi Coefficient (Co-occurrence Strength)", fontsize=11)
    ax2.set_ylabel("Phylogenetic Scatter Score", fontsize=11)
    ax2.set_title("Co-occurrence Strength vs Phylogenetic Scatter",
                  fontsize=12, fontweight="bold")
    from matplotlib.patches import Patch
    ax2.legend(handles=[Patch(color=c, label=k) for k, c in color_map.items()], fontsize=9)
    ax2.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = f"{outdir}/scatter_score_distribution.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"[Output] Scatter score distribution → {out_path}")
    plt.close()


# ──────────────────────────────────────────────────────────
# Main runner
# ──────────────────────────────────────────────────────────

def run_cooccurrence_analysis(
    dist_df: pd.DataFrame,
    amr_matrix: pd.DataFrame,
    drug_class_map: dict,
    outdir: str,
    phi_threshold: float = 0.3,
    scatter_threshold: float = 0.7
):
    import os
    os.makedirs(outdir, exist_ok=True)

    print("[Module 2] Building AMR co-occurrence network...")
    G = build_cooccurrence_network(amr_matrix, phi_threshold=phi_threshold)
    print(f"  Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

    if G.number_of_edges() == 0:
        print("[WARN] No co-occurring pairs found. Lower phi_threshold.")
        return None, pd.DataFrame()

    print("[Module 2] Classifying co-occurrences by phylogenetic scatter...")
    results = classify_cooccurrences(G, amr_matrix, dist_df, scatter_threshold)

    results.to_csv(f"{outdir}/cooccurrence_classification.csv", index=False)
    print(f"[Output] Co-occurrence table → {outdir}/cooccurrence_classification.csv")

    print(f"\n[Summary]")
    print(results["classification"].value_counts().to_string())

    conserved = results[results["classification"] == "CONSERVED"]
    if not conserved.empty:
        print(f"\nTop CONSERVED pairs (likely mobile element co-transfer):")
        print(conserved.head(5)[["gene_1", "gene_2", "phi", "scatter_score"]].to_string(index=False))

    scattered = results[results["classification"] == "SCATTERED"]
    if not scattered.empty:
        print(f"\nTop SCATTERED pairs (independent acquisition):")
        print(scattered.head(5)[["gene_1", "gene_2", "phi", "scatter_score"]].to_string(index=False))

    plot_cooccurrence_network(G, drug_class_map, outdir)
    plot_scatter_score_distribution(results, outdir)

    return G, results
