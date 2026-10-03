"""Screen a pair of PR-10 sequences for predicted structural deviation.

    python3 screen.py seq1.fasta seq2.fasta
    python3 screen.py --betv1 seq2.fasta
    cat two.fasta | python3 screen.py -

Needs pdb_all/ (run fetch_data.py) so it can find the epitope patch and load the
trained model. Trains on the fly from pairs.csv if that is present.

This estimates how far apart two structures will lie. It does NOT predict
whether IgE will bind them; that needs clinical data.
"""

from __future__ import annotations

import sys

import numpy as np

from model import FEATURES, TARGET, _gbr, features_for_pair, load
from pipeline import epitope_index_per_protein, stage1_proteins, stage2_epitope

BETV1_UNIPROT = "P15494"
CUTOFF = 4.5
TRIAGE_HIGH, TRIAGE_LOW = 2.0, 3.5  # A; below = close, above = far


def read_fasta(path):
    if path == "-":
        lines = sys.stdin.read().splitlines()
    else:
        lines = open(path).read().splitlines()
    seq, name, chunks = "", "", []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq:
                chunks.append((name, seq))
            name, seq = line[1:], ""
        else:
            seq += line
    if seq:
        chunks.append((name, seq))
    return chunks


def one(chunks):
    """Single-sequence input: FASTA, or a bare sequence on the command line."""
    if chunks:
        name, seq = chunks[0]
        return name or "input", seq.upper()
    return "input", chunks


def triage(rmsd):
    if rmsd < TRIAGE_HIGH:
        return "close - structurally near-identical, prioritise for an IgE assay"
    if rmsd < TRIAGE_LOW:
        return "moderate - same fold, region differences"
    return "far - same family, but architecture has diverged"


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    betv1 = "--betv1" in argv
    if not args:
        print(__doc__)
        return 1

    proteins = stage1_proteins()
    epitope, _ = stage2_epitope()
    reference = next(p for p in proteins
                     if p.uniprot == BETV1_UNIPROT and p.method == "X-RAY DIFFRACTION")
    index = epitope_index_per_protein(epitope[CUTOFF], reference.seq, proteins)

    loaded = []
    for path in args:
        name, seq = one(read_fasta(path))
        loaded.append((name, seq))

    if betv1:
        loaded.insert(0, (f"Bet v 1 ({reference.entry})", reference.seq))

    model = _gbr().fit(load()[FEATURES].values, load()[TARGET].values)

    print(f"epitope: {len(epitope[CUTOFF])} IgE-contact positions from 9Y0A, "
          f"anchored on {reference.entry}")
    for (n1, s1), (n2, s2) in zip(loaded, loaded[1:]):
        f = features_for_pair(s1, s2, index[n1] if n1 in index else
                              index.get(reference.pid, set()))
        pred = float(model.predict(np.array([[f[c] for c in FEATURES]]))[0])
        print(f"\n{n1} ({len(s1)} aa)  vs  {n2} ({len(s2)} aa)")
        print(f"  identity over shorter seq {f['identity_short']*100:.1f}%   "
              f"epitope identity {f['epi_identity']*100:.0f}%   "
              f"gaps {f['gap_frac']*100:.1f}%")
        print(f"  predicted RMSD {pred:.2f} A   ->  {triage(pred)}")
    print("\nPredicted structural deviation only. Not an IgE cross-reactivity call.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))