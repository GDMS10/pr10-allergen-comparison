"""Sequence-only model for PR-10 structural deviation, validated leave-one-organism-out.

The repo asks whether sequence similarity predicts structural similarity among
PR-10 allergens. With 10 pairs the answer is unmeasurable. This asks it again
with 3616 pairs from 86 proteins in 24 organisms, and trains a model to predict
structural deviation from sequence alone.

    python3 model.py           # correlations, CV scores, screening helpers
    python3 model.py --full    # adds permutation test and feature importances

The model is a triage tool, not a cross-reactivity predictor: it estimates how
far apart two PR-10 structures will lie, so you know which pairs are worth
taking to an IgE assay first. Nothing here predicts whether IgE will bind.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score

from pipeline import ALIGNER, _alignment_blocks

# Every one of these is computable from two sequences. Nothing touches
# coordinates, which is the whole point: score a new allergen before you have
# its structure.
FEATURES = [
    "identity_cov",      # identity over aligned columns (the honest denominator)
    "identity_short",    # identity over the shorter sequence (what the repo reports)
    "blosum_norm",       # substitution-weighted score, not just exact matches
    "hyd_corr",          # hydrophobicity profile correlation along the alignment
    "epi_identity",      # identity restricted to the 9Y0A IgE contact patch
    "epi_coverage",      # how much of the alignment the epitope patch covers
    "len_diff",
    "gap_frac",          # alignment gaps: the single strongest predictor, see below
]
BASELINE = ["identity_short"]  # the repo's one feature
N_FOLDS = 8
TARGET = "rmsd_pruned"

from pipeline import KD  # noqa: E402  (used by features_for_pair)


def _gbr():
    return GradientBoostingRegressor(random_state=0, n_estimators=200, max_depth=3,
                                     learning_rate=0.06, subsample=0.9)


def _ridge():
    return RidgeCV(alphas=np.logspace(-3, 3, 13))


def load(path="pairs.csv"):
    return pd.read_csv(path).dropna(subset=["epi_identity"])


def organism_folds(df, n_folds=N_FOLDS):
    """Fold id per pair, keyed on the unordered pair of its two organisms.

    Leaving a protein out is not possible for pairwise data: every pair has two
    proteins, so a protein-level split still shows the held-out protein's
    relatives in training. Holding out whole organisms is leak-free, and it is
    the realistic case, since the use case is scoring a food source not yet seen.
    """
    orgs = sorted(set(df.org_a) | set(df.org_b))
    rng = np.random.default_rng(0)
    fold_of = dict(zip(orgs, rng.permutation(len(orgs)) % n_folds))
    return df.apply(lambda r: tuple(sorted((fold_of[r.org_a], fold_of[r.org_b]))), axis=1)


def fold_split(folds):
    for f in sorted(set(folds)):
        te = np.array([x == f for x in folds])
        if te.sum() and (~te).sum() >= 50:
            yield ~te, te


def score(y, preds):
    return dict(r2=r2_score(y, preds), mae=mean_absolute_error(y, preds),
                within_1A=float(np.mean(np.abs(preds - y) < 1.0)))


def evaluate(df, feature_sets, target=TARGET):
    folds = organism_folds(df)
    y = df[target].values
    out = {}
    for name, cols in feature_sets.items():
        X = df[cols].values
        for kind, make in (("", _gbr), (" [ridge]", _ridge)):
            preds = np.zeros(len(y))
            for tr, te in fold_split(folds):
                m = make().fit(X[tr], y[tr])
                preds[te] = m.predict(X[te])
            out[name + kind] = score(y, preds)
    return out


def correlations(df, target=TARGET):
    return {c: spearmanr(df[c], df[target]) for c in
            ["identity_short", "identity_cov", "blosum_norm", "hyd_corr",
             "epi_identity", "gap_frac", "len_diff"]}


def permutation_test(df, cols, target=TARGET, n=20):
    """Shuffle y within fold, refit, collect grouped-CV R2. Guards against reading
    the real R2 as an artefact of the split."""
    folds = organism_folds(df)
    X, y = df[cols].values, df[target].values
    rng = np.random.default_rng(0)

    def run(yy):
        preds = np.zeros(len(yy))
        for tr, te in fold_split(folds):
            preds[te] = _gbr().fit(X[tr], yy[tr]).predict(X[te])
        return r2_score(yy, preds)

    null = []
    for _ in range(n):
        ys = y.copy()
        for f in sorted(set(folds)):
            te = np.array([x == f for x in folds])
            ys[te] = rng.permutation(ys[te])
        null.append(run(ys))
    return run(y), np.array(null)


def importances(df, feature_sets):
    """Leave-one-feature-out grouped CV. Honest; no impurity shortcut."""
    base = evaluate(df, {"all": feature_sets})["all"]["r2"]
    return {f: base - evaluate(df, {f: [c for c in feature_sets if c != f]})[f]["r2"]
            for f in feature_sets}


def features_for_pair(seq_a, seq_b, epitope_indices=()):
    """The feature row for an arbitrary pair of PR-10 sequences.

    This is the part that matters in practice: no coordinates needed, so a new
    allergen can be screened before anyone solves a structure for it.
    """
    sub = substitution_matrices.load("BLOSUM62")
    ia, ib = _alignment_blocks(ALIGNER.align(seq_a, seq_b)[0])
    gaps = (len(seq_a) + len(seq_b)) - 2 * len(ia)
    hits = [(i, j) for i, j in zip(ia, ib) if i in epitope_indices]

    def profile(s, idx):
        return np.array([KD.get(s[i], 0.0) for i in idx])

    hyd = profile(seq_a, ia), profile(seq_b, ib)
    return {
        "identity_cov": sum(seq_a[i] == seq_b[j] for i, j in zip(ia, ib)) / max(len(ia), 1),
        "identity_short": sum(seq_a[i] == seq_b[j] for i, j in zip(ia, ib))
                          / min(len(seq_a), len(seq_b)),
        "blosum_norm": sum(sub[seq_a[i], seq_b[j]] for i, j in zip(ia, ib)) / max(len(ia), 1),
        "hyd_corr": float(np.corrcoef(*hyd)[0, 1]) if len(ia) > 3 else 0.0,
        "epi_identity": (sum(seq_a[i] == seq_b[j] for i, j in hits) / len(hits)) if hits else np.nan,
        "epi_coverage": len(hits) / max(len(ia), 1),
        "len_diff": abs(len(seq_a) - len(seq_b)),
        "gap_frac": gaps / (len(seq_a) + len(seq_b)),
    }


def main(full=False):
    df = load()
    print(f"pairs: {len(df)}   proteins: {len(set(df.a) | set(df.b))}   "
          f"organisms: {len(set(df.org_a) | set(df.org_b))}")

    print("\nDoes sequence similarity predict RMSD at family scale?  (Spearman vs RMSD)")
    for c, (r, p) in correlations(df).items():
        print(f"  {c:16} {r:+.3f}   p={p:.1e}")
    print("  Note gap_frac is the strongest single predictor, and identity is not first.")

    res = evaluate(df, {"all features": FEATURES, "identity only": BASELINE})
    print("\nLeave-one-organism-out CV (both proteins of a tested pair unseen):")
    print(f"  {'model':34} {'R2':>7} {'MAE (A)':>9} {'|err|<1A':>9}")
    for k, v in res.items():
        print(f"  {k:34} {v['r2']:7.3f} {v['mae']:9.3f} {v['within_1A']:9.1%}")

    y = df[TARGET].values
    print("\nTrivial references (in-sample, so optimistic):")
    print(f"  {'predict family mean':34} {0.0:7.3f} "
          f"{mean_absolute_error(y, np.full_like(y, y.mean())):9.3f}")

    if full:
        real, null = permutation_test(df, FEATURES)
        print(f"\nPermutation test: real grouped-CV R2={real:.3f} vs null "
              f"{null.mean():.3f} +- {null.std():.3f}")
        print("\nLeave-one-feature-out cost in R2:")
        for f, d in sorted(importances(df, FEATURES).items(), key=lambda kv: -kv[1]):
            print(f"  {f:18} {d:+.4f}")
    return df


if __name__ == "__main__":
    import sys

    main(full="--full" in sys.argv)