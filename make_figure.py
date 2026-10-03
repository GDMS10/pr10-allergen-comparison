"""Figure for the extended analysis: what n=3616 looks like where n=10 did not.

    python3 make_figure.py

Writes identity_vs_rmsd_family.png next to the original identity_vs_rmsd.png, so
the two claims can be compared side by side.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from model import FEATURES, TARGET, _gbr, load, organism_folds

ORIGINAL_PDB = {"Bet v 1": "1BV1", "Mal d 1": "5MMU", "Pru av 1": "1E09",
                "Api g 1": "2BK0", "Dau c 1": "2WQL"}
# entry -> the chain stage1_proteins picked, so the original 10 pairs can be
# located inside the new dataset rather than re-aligned.
ORIGINAL_PIDS = {"1BV1": "1BV1_A", "5MMU": "5MMU_A", "1E09": "1E09_A",
                 "2BK0": "2BK0_A", "2WQL": "2WQL_A"}


def original_pairs(df):
    """The repo's 10 pairs, expressed in the new dataset's ids."""
    out = []
    for name, eid in ORIGINAL_PDB.items():
        for other, oeid in ORIGINAL_PDB.items():
            if eid >= oeid:
                continue
            a, b = ORIGINAL_PIDS[eid], ORIGINAL_PIDS[oeid]
            hit = df[((df.a == a) & (df.b == b)) | ((df.a == b) & (df.b == a))]
            if len(hit):
                out.append((name, other, hit.iloc[0]))
    return out


def main():
    df = load()
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))

    # 1. the same scatter as the original figure, at family scale
    ax = axes[0]
    ax.scatter(df.identity_short * 100, df[TARGET], s=9, alpha=0.28,
               color="#4C72B0", edgecolors="none", label=f"all {len(df)} pairs")
    originals = original_pairs(df)
    ax.scatter([r[2].identity_short * 100 for r in originals],
               [r[2][TARGET] for r in originals], s=95, marker="D", zorder=5,
               color="#C44E52", edgecolors="black", linewidths=0.6,
               label="the original 10 pairs")
    rho, p = stats.spearmanr(df.identity_short, df[TARGET])
    ax.set_title(f"Sequence identity vs structure, n={len(df)}\n"
                 f"Spearman {rho:+.2f} (p<1e-300) — identity does predict",
                 fontsize=11)
    ax.set_xlabel("Sequence identity over shorter sequence (%)")
    ax.set_ylabel("Outlier-pruned CA RMSD (Å)")
    ax.legend(frameon=True, fontsize=9, loc="upper right")

    # 2. the feature that actually wins
    ax = axes[1]
    ax.scatter(df.gap_frac * 100, df[TARGET], s=9, alpha=0.28,
               color="#55A868", edgecolors="none")
    rho_gap, _ = stats.spearmanr(df.gap_frac, df[TARGET])
    ax.set_title(f"Alignment gaps predict structure better\n"
                 f"Spearman {rho_gap:+.2f} vs {rho:+.2f} for identity",
                 fontsize=11)
    ax.set_xlabel("Alignment gap fraction (%)")
    ax.set_ylabel("Outlier-pruned CA RMSD (Å)")

    # 3. model predicted vs observed, leave-one-organism-out
    ax = axes[2]
    folds = organism_folds(df)
    X = df[FEATURES].values
    y = df[TARGET].values
    preds = np.zeros(len(y))
    for f in sorted(set(folds)):
        te = np.array([x == f for x in folds])
        preds[te] = _gbr().fit(X[~te], y[~te]).predict(X[te])
    lim = [0, max(y.max(), preds.max()) * 1.03]
    ax.hexbin(preds, y, gridsize=45, bins="log", cmap="viridis", mincnt=1)
    ax.plot(lim, lim, "--", color="black", lw=1, label="perfect prediction")
    mae = np.mean(np.abs(preds - y))
    r2 = 1 - ((y - preds) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    ax.set_title(f"Sequence-only model, leave-one-organism-out\n"
                 f"R²={r2:.2f}, MAE={mae:.2f} Å — no coordinates used",
                 fontsize=11)
    ax.set_xlabel("Predicted RMSD from sequence alone (Å)")
    ax.set_ylabel("Observed RMSD (Å)")
    ax.set_xlim(lim), ax.set_ylim(lim)
    ax.legend(fontsize=9, loc="upper left")

    fig.suptitle("PR-10 allergens: sequence similarity does predict structural "
                 "similarity, and gaps predict it better than identity",
                 fontsize=13, y=1.0)
    fig.tight_layout()
    out = "identity_vs_rmsd_family.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"  n={len(df)} pairs, {len(set(df.a) | set(df.b))} proteins, "
          f"{len(set(df.org_a) | set(df.org_b))} organisms")
    print(f"  Spearman identity {rho:+.3f}, gaps {rho_gap:+.3f}")
    print(f"  model R2={r2:.3f}, MAE={mae:.2f} A, {np.mean(np.abs(preds-y) < 1)*100:.0f}% within 1 A")


if __name__ == "__main__":
    main()