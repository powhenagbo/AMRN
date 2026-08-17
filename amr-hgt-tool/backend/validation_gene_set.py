"""
validation_gene_set.py
========================
A deliberately curated (not opportunistically discovered) set of AMR
genes with well-established literature consensus on their transmission
mode — used to systematically check the tool's Mantel-panel HGT/CLONAL
calls against ground truth, rather than relying on findings stumbled
into while doing something else (as gyrA/parC and aadA5 were).

Sources for each classification are standard, widely-cited AMR biology
(the mobile-genetic-element literature for the HGT list; the target-
site-mutation and intrinsic-resistance literature for the clonal list),
summarized here rather than fabricated — but this file should be treated
as a starting point to expand/correct, not a final authority. If your
own panel or literature review disagrees with an entry, that disagreement
is itself a useful validation finding.
"""

# (gene_name_patterns, expected_classification, justification)
# gene_name_patterns: list of substrings to match against panel gene
# names (case-insensitive) — some genes appear under multiple naming
# conventions (e.g. RGI's CARD-style long names vs. short symbols).
KNOWN_HGT_GENES = [
    (["CTX-M"], "Plasmid-borne ESBL — among the most extensively documented mobile AMR genes; spread globally via conjugative plasmids across Enterobacteriaceae."),
    (["NDM-1", "NDM-"], "Plasmid-borne metallo-beta-lactamase — textbook case of rapid global dissemination via horizontal transfer, first reported 2008, spread to 100+ countries within a decade."),
    (["KPC-"], "Plasmid-borne carbapenemase — classic mobile resistance determinant, spread via conjugative plasmids and transposons (Tn4401)."),
    (["tet(A)", "tet(B)", "tetA", "tetB"], "Transposon-borne tetracycline efflux — canonical mobile resistance gene, associated with Tn10 and related transposons."),
    (["qnrA", "qnrB", "qnrS"], "Plasmid-mediated quinolone resistance (PMQR) — by definition a horizontally-acquired mechanism, distinct from chromosomal gyrA/parC target-site mutations."),
    (["sul1", "sul2"], "Sulfonamide resistance — sul1 is a class 1 integron cassette gene, canonically mobile; sul2 is plasmid-associated."),
    (["dfrA"], "Trimethoprim resistance — dfrA-family genes are classic integron gene cassettes, horizontally mobile by definition."),
    (["aadA"], "Aminoglycoside adenylyltransferase — canonical integron cassette gene family (streptomycin/spectinomycin resistance)."),
    (["aac(6", "aac6"], "Aminoglycoside acetyltransferase — frequently plasmid/integron-associated, well-documented mobile determinant."),
    (["qacE"], "Quaternary ammonium compound resistance — class 1 integron cassette gene, near-universally co-located with sul1/other cassette genes on the same mobile element."),
]

KNOWN_CLONAL_GENES = [
    (["gyrA"], "DNA gyrase target-site mutation — fluoroquinolone resistance arises via spontaneous point mutation in the resident chromosomal gene, then clonal expansion; not horizontally transferred as a whole gene in Enterobacteriaceae."),
    (["gyrB"], "Same mechanism class as gyrA — chromosomal target-site mutation, not acquisition."),
    (["parC"], "Topoisomerase IV target-site mutation — same biology as gyrA, chromosomal point mutation and clonal spread."),
    (["parE"], "Same mechanism class as parC — chromosomal target-site mutation."),
    (["rpoB"], "RNA polymerase target-site mutation (rifampin resistance) — classic chromosomal point-mutation resistance mechanism, not horizontally transferred."),
    (["marA"], "Chromosomal multidrug-resistance regulator — intrinsic regulatory gene present in essentially all E. coli, not an acquired determinant."),
    (["soxS"], "Chromosomal regulon regulator — same intrinsic-regulatory category as marA."),
    (["acrA", "acrB"], "Core structural genes of the AcrAB-TolC intrinsic efflux pump — present as baseline biology in essentially all Enterobacteriaceae, not a horizontally acquired trait."),
    (["tolC"], "Outer-membrane channel shared by multiple intrinsic efflux systems — core chromosomal gene, not acquired."),
]


def check_gene_name(panel_gene_name: str, patterns: list) -> bool:
    """Case-insensitive substring match against any of the given patterns."""
    name_lower = panel_gene_name.lower()
    return any(p.lower() in name_lower for p in patterns)
