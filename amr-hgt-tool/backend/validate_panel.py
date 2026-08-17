"""
validate_panel.py
====================
Runs the curated known-HGT / known-clonal gene set (validation_gene_set.py)
against your actual reference panel's Mantel classifications, and reports
matches, mismatches, and coverage gaps — a structured validation pass,
not an opportunistic one-off check.

Usage:
    cd ~/AMRN/amr-hgt-tool/backend
    python validate_panel.py                          # uses panel/v1
    python validate_panel.py --panel-dir panel/v1_amrfinder
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validation_gene_set import KNOWN_HGT_GENES, KNOWN_CLONAL_GENES, check_gene_name


def run_validation(panel_dir: str):
    csv_path = os.path.join(panel_dir, "transmission_classification.csv")
    if not os.path.exists(csv_path):
        print(f"No transmission_classification.csv found at {csv_path}")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    print(f"Loaded panel: {csv_path} ({len(df)} genes)\n")

    def check_group(curated_list, expected_label):
        results = []
        for patterns, justification in curated_list:
            matches = df[df["gene"].apply(lambda g: check_gene_name(str(g), patterns))]
            if matches.empty:
                results.append({
                    "patterns": patterns, "found": False, "expected": expected_label,
                    "justification": justification,
                })
            else:
                for _, row in matches.iterrows():
                    observed = row["classification"]
                    results.append({
                        "patterns": patterns, "found": True, "expected": expected_label,
                        "gene": row["gene"], "observed": observed,
                        "mantel_r": row["mantel_r"], "p_value": row["p_value"],
                        "n_carriers": row["n_carriers"],
                        "match": observed == expected_label,
                        "justification": justification,
                    })
        return results

    hgt_results = check_group(KNOWN_HGT_GENES, "HGT")
    clonal_results = check_group(KNOWN_CLONAL_GENES, "CLONAL")
    all_results = hgt_results + clonal_results

    found = [r for r in all_results if r["found"]]
    not_found = [r for r in all_results if not r["found"]]
    matches = [r for r in found if r["match"]]
    mismatches = [r for r in found if not r["match"]]

    print("=" * 70)
    print(f"COVERAGE: {len(found)} of {len(all_results)} curated genes found in this panel")
    print(f"ACCURACY (of those found): {len(matches)}/{len(found)} matched literature expectation")
    print("=" * 70)

    if matches:
        print(f"\n✓ MATCHES ({len(matches)}):")
        for r in matches:
            print(f"  {r['gene']:40s} expected={r['expected']:7s} observed={r['observed']:7s} "
                  f"r={r['mantel_r']:+.3f} p={r['p_value']:.4f} n={r['n_carriers']}")

    if mismatches:
        print(f"\n✗ MISMATCHES ({len(mismatches)}) — worth investigating each of these:")
        for r in mismatches:
            print(f"  {r['gene']:40s} expected={r['expected']:7s} observed={r['observed']:7s} "
                  f"r={r['mantel_r']:+.3f} p={r['p_value']:.4f} n={r['n_carriers']}")
            print(f"    why expected {r['expected']}: {r['justification']}")

    if not_found:
        print(f"\n— NOT IN PANEL ({len(not_found)}) — no coverage to validate against:")
        for r in not_found:
            print(f"  {'/'.join(r['patterns']):40s} (expected {r['expected']})")

    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-dir", default="panel/v1")
    args = parser.parse_args()
    run_validation(args.panel_dir)
