"""
check_islands.py
==================
Diagnostic: calls run_kali_islands() directly against a genome FASTA and
reports how many islands were found, genome-wide — bypassing the
gene-overlap step entirely.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "pipeline"))

from wrapper import run_kali_islands

if len(sys.argv) < 2:
    print("Usage: python check_islands.py /path/to/genome.fasta")
    sys.exit(1)

fasta_path = sys.argv[1]
print(f"Running island detection on: {fasta_path}")

islands_df = run_kali_islands(fasta_path)

print(f"\n{len(islands_df)} islands detected, genome-wide.\n")

if len(islands_df) > 0:
    print("Top 10 by Z-score:")
    top = islands_df.sort_values("max_zscore", ascending=False).head(10)
    print(top.to_string(index=False))
else:
    print("Zero islands — the detection logic itself isn't flagging anything on this genome.")
