# PR-10 allergen comparison

Does sequence similarity predict structural similarity among PR-10 plant allergens?

This compares five PR-10 allergens pairwise — percent sequence identity against backbone
RMSD — and finds that the two measures do not track each other closely.

![Sequence identity vs RMSD](identity_vs_rmsd.png)

## Background

People with birch pollen allergy often react to raw apple, cherry, celery and carrot. The
reaction is usually limited to the mouth and lips, and it disappears when the food is cooked.
This is oral allergy syndrome (also called pollen-food allergy syndrome).

The cause is cross-reactivity. IgE antibodies raised against the birch pollen allergen
**Bet v 1** also bind structurally similar proteins in those foods. All of them belong to the
**PR-10** protein family and share the same fold.

The epitopes involved are **conformational**: the antibody recognises a patch of surface
residues that are far apart in the sequence but adjacent in the folded structure. That is why
cooking abolishes the reaction — heat unfolds the protein and the patch stops existing.

Because recognition depends on the folded surface rather than the sequence, sequence identity
may be a poor proxy for whether two allergens cross-react. This project tests how well the two
measures agree.

## Data

Five experimental structures from the RCSB PDB:

| PDB | Allergen | Source | Method |
|---|---|---|---|
| 1BV1 | Bet v 1 | birch pollen | X-ray |
| 5MMU | Mal d 1 | apple | NMR |
| 1E09 | Pru av 1 | cherry | NMR |
| 2BK0 | Api g 1 | celery | X-ray |
| 2WQL | Dau c 1 | carrot | X-ray |

## Method

For each of the 10 pairs:

1. Extract the amino acid sequence and the CA atoms of one chain (model 0).
2. Align the two sequences globally with **BLOSUM62**, gap open −11, gap extend −1.
3. Walk the aligned blocks to build two equal-length lists of paired CA atoms, one per
   aligned residue. Residues opposite a gap are excluded.
4. **Percent identity** = identical aligned pairs ÷ length of the shorter sequence × 100.
5. **RMSD** = superimpose the paired CA atoms (Biopython `Superimposer`) and take the
   root-mean-square deviation.

Method choices worth stating: identity is computed over the shorter sequence; the
superposition uses every aligned pair with no outlier pruning, so values are slightly higher
than tools that prune iteratively (e.g. ChimeraX `matchmaker`); where residues had alternate
conformations, the highest-occupancy position was used.

## Result

The full table is in [`allergen_comparison.csv`](allergen_comparison.csv). Two pairs make the
point:

| Pair | Sequence identity | RMSD |
|---|---|---|
| Api g 1 – Dau c 1 (celery/carrot) | 80.3% | **0.69 Å** — most similar structure |
| Mal d 1 – Pru av 1 (apple/cherry) | **84.8%** — most similar sequence | 2.11 Å |

The most sequence-similar pair is not the most structurally similar. The same disagreement
appears from birch: Bet v 1 – Api g 1 has lower sequence identity than Bet v 1 – Mal d 1
(42.5% vs 56.3%) but a lower RMSD (1.49 Å vs 2.11 Å).

## Limitations

- **Structure determination method is a confounder.** The lowest RMSD pair (celery/carrot) is
  the only X-ray/X-ray comparison. Every pair involving the NMR entries (apple, cherry) lands
  near 2 Å. NMR ensembles are systematically less tightly defined, so part of this spread is
  likely methodological rather than biological. Regenerating all five structures with
  AlphaFold would give one uniform method and test this.
- **Five proteins, one family, ten pairs.** A small sample.
- **No clinical data.** This shows that two measures disagree; it does not show which one
  predicts allergy. Validating that needs cross-reactivity rates from clinical studies.
- **One isoform each.** Bet v 1 in particular exists as many isoforms that differ in IgE
  reactivity, so these numbers describe specific isoforms, not the allergens in general.

## Next

Epitope-level analysis. `9Y0A` is a crystal structure of Bet v 1 bound to a human IgE Fab, so
the residues the antibody actually contacts can be measured rather than predicted. The question
is whether those specific positions are more conserved in the food allergens than the rest of
the protein is.

## Files

| File | |
|---|---|
| `allergen_functions.py` | the analysis functions |
| `data_creation.ipynb` | builds `allergen_comparison.csv` from the PDB files |
| `figure.ipynb` | produces `identity_vs_rmsd.png` |
| `pdb_files/` | the five structures, as downloaded |

## Running it

```
python3 -m venv .venv
source .venv/bin/activate
pip install biopython pandas matplotlib seaborn jupyter
```

Then run `data_creation.ipynb` followed by `figure.ipynb`.

## Notes

Written by MS as a self-directed project.
