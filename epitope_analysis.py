"""The README's "Next" step: are the IgE-contact positions more conserved than the
rest of the protein? Answered with 3616 pairs instead of one.

The README proposes using 9Y0A (Bet v 1.0101 bound to the human IgE Fab 2H22) to
measure whether the residues the antibody actually contacts are more conserved in
the food allergens than the protein as a whole, on the theory that recognition
tracks the epitope rather than average sequence similarity.

The epitope is extracted here, then compared against the rest of the same
alignment with a paired test, plus random patches of equal size as a control.

Result: the epitope positions are significantly LESS conserved than background,
not more. See the module docstring in the README for what that implies.

    python3 epitope_analysis.py
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy import stats
from Bio.Align import substitution_matrices

from pipeline import ALIGNER, _alignment_blocks, epitope_index_per_protein, stage1_proteins, stage2_epitope

CUTOFF = 4.5          # A, all-atom contact between Bet v 1 and the Fab
BETV1 = "P15494"
N_RANDOM_PATCHES = 20


def conservation_scores(proteins, epitope_idx):
    """Per pair: mean BLOSUM score at the epitope positions vs everywhere else.

    BLOSUM rather than raw identity because 16 positions give a very coarse
    identity fraction; the substitution matrix keeps resolution.
    """
    sub = substitution_matrices.load("BLOSUM62")
    epi, bg = [], []
    for a, b in combinations(proteins, 2):
        if abs(a.length - b.length) > 40:
            continue
        ia, ib = _alignment_blocks(ALIGNER.align(a.seq, b.seq)[0])
        if len(ia) < 60:
            continue
        on_epi = [(i, j) for i, j in zip(ia, ib) if i in epitope_idx[a.pid]]
        off_epi = [(i, j) for i, j in zip(ia, ib) if i not in epitope_idx[a.pid]]
        if not on_epi or len(off_epi) < 50:
            continue
        epi.append(np.mean([sub[a.seq[i], b.seq[j]] for i, j in on_epi]))
        bg.append(np.mean([sub[a.seq[i], b.seq[j]] for i, j in off_epi]))
    return np.array(epi), np.array(bg)


def report(label, epi, bg):
    diff = epi - bg
    t, p = stats.ttest_rel(epi, bg)
    w = stats.wilcoxon(epi, bg).pvalue
    print(f"  {label:26} n={len(epi):5}  epitope {epi.mean():.3f}  "
          f"background {bg.mean():.3f}  diff {diff.mean():+.3f} "
          f"+-{diff.std(ddof=1)/np.sqrt(len(diff)):.3f}  t-test p={p:.1e}  Wilcoxon p={w:.1e}")
    return diff.mean()


def main():
    proteins = stage1_proteins()
    epitope, _ = stage2_epitope()
    residues = epitope[CUTOFF]
    ref = next(p for p in proteins
               if p.uniprot == BETV1 and p.method == "X-RAY DIFFRACTION")
    print(f"epitope: {len(residues)} Bet v 1 residues within {CUTOFF} A of the IgE Fab "
          f"(9Y0A, chain C), from {ref.entry}")
    print(f"positions: {[r[0] for r in residues]}")
    print(f"residues:  {''.join(r[1] for r in residues)}\n")

    idx = epitope_index_per_protein(residues, ref.seq, proteins)
    print("Is the epitope more conserved than its own background?")
    for cutoff in sorted(epitope):
        i = epitope_index_per_protein(epitope[cutoff], ref.seq, proteins)
        report(f"{cutoff} A contact", *conservation_scores(proteins, i))

    print("\nControl: random patches of 16 residues, same size as the real epitope")
    rng = np.random.default_rng(0)
    diffs = []
    for _ in range(N_RANDOM_PATCHES):
        pick = sorted(rng.choice(len(ref.seq), len(residues), replace=False).tolist())
        fake = [(n + 1, ref.seq[n], 0.0) for n in pick]
        i = epitope_index_per_protein(fake, ref.seq, proteins)
        e, b = conservation_scores(proteins, i)
        diffs.append((e - b).mean())
    diffs = np.array(diffs)
    print(f"  {'random patch, mean':26} {diffs.mean():+.3f} "
          f"(range {diffs.min():+.3f}..{diffs.max():+.3f} over {N_RANDOM_PATCHES} draws)")

    print("\nInterpretation: the real epitope is a genuine outlier against random "
          "patches of the same size, but in the opposite direction to the "
          "hypothesis. The IgE-contact residues sit in the most variable part of "
          "the fold, which is consistent with IgE binding being specific to "
          "allergenicity rather than to general sequence conservation.")


if __name__ == "__main__":
    main()