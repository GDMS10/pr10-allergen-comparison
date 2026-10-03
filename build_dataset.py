"""Rebuild pairs.csv from the downloaded structures.

    python3 fetch_data.py     # once: 133 entries into pdb_all/
    python3 build_dataset.py  # ~20 s: 86 proteins, 3616 pairs

Two deviations from allergen_comparison.csv, both deliberate:

  * NMR entries are averaged over their whole ensemble instead of using model 0,
    so the target is not one arbitrary conformer.
  * The IgE epitope patch comes from 9Y0A (Bet v 1 + Fab 2H22), which is what
    turns a plain identity number into a cross-reactivity feature.
"""

from __future__ import annotations

import pandas as pd

from pipeline import epitope_index_per_protein, stage1_proteins, stage2_epitope, stage3_pairs

EPITOPE_CUTOFF = 4.5  # A, all-atom contact with the Fab
BETV1_UNIPROT = "P15494"


def main():
    proteins = stage1_proteins()
    print(f"{len(proteins)} unique PR-10 proteins with coordinates")
    by_method = {}
    for p in proteins:
        by_method[p.method] = by_method.get(p.method, 0) + 1
    print(f"  {by_method}")

    epitope, _ = stage2_epitope()
    residues = epitope[EPITOPE_CUTOFF]
    print(f"{len(residues)} Bet v 1 residues within {EPITOPE_CUTOFF} A of the IgE Fab: "
          f"{[r[0] for r in residues]}")

    reference = next(p for p in proteins
                     if p.uniprot == BETV1_UNIPROT and p.method == "X-RAY DIFFRACTION")
    index = epitope_index_per_protein(residues, reference.seq, proteins)
    print(f"epitope reference: {reference.pid} ({reference.entry})")

    rows = stage3_pairs(proteins, index)
    df = pd.DataFrame(rows)
    df.to_csv("pairs.csv", index=False)
    print(f"\nwrote pairs.csv: {len(df)} pairs, {df.shape[1]} columns")
    print(df[["identity_short", "epi_identity", "rmsd", "rmsd_pruned"]].describe().round(2))


if __name__ == "__main__":
    main()