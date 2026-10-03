"""Checks for pipeline.py. Run: python3 test_pipeline.py

No test framework: these are the assertions that caught the real bugs, so they
stay small enough to read in one screen.
"""

from __future__ import annotations

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

from pipeline import ALIGNER, _kabsch, _rmsd_pruned, _alignment_blocks


class Atom:
    """Minimal stand-in so Bio's Superimposer will take plain coordinates."""

    def __init__(self, coord):
        self.coord = np.asarray(coord, dtype=float)

    def get_coord(self):
        return self.coord


def _bio_rmsd(fixed, moving):
    s = SVDSuperimposer()
    s.set(fixed, moving)
    s.run()
    return s.get_rms()


def test_kabsch_matches_biopython():
    """My Kabsch must equal Bio's on random and reflected point sets.

    Two bugs lived here: composing the rotation as V @ U.T instead of U @ Vt, and
    negating a row of U for the reflection case, which leaves the determinant
    negative and the RMSD wrong. Random data plus an explicit reflection catches
    both; a hand-checked pair does not.
    """
    rng = np.random.default_rng(0)
    for trial in range(50):
        fixed = rng.random((rng.integers(20, 120), 3)) * 30
        moving = rng.random((len(fixed), 3)) * 30
        assert abs(_kabsch(fixed, moving)[0] - _bio_rmsd(fixed, moving)) < 1e-9, trial

        # Force the reflection branch and confirm it is still correct.
        if trial % 5 == 0:
            mirrored = fixed * np.array([1.0, 1.0, -1.0])
            rot = _kabsch(fixed, mirrored)[2]
            assert np.linalg.det(rot) > 0, "reflection not corrected"


def test_kabsch_on_a_known_pair():
    """Identity under a known rigid motion, and zero for a perfect copy."""
    rng = np.random.default_rng(1)
    P = rng.random((50, 3)) * 20
    R = np.linalg.qr(rng.random((3, 3)))[0]
    Q = P @ R + np.array([7.0, -3.0, 11.0])
    assert _kabsch(P, Q)[0] < 1e-9
    assert _kabsch(P, P)[0] < 1e-12


def test_rmsd_pruned_never_exceeds_raw():
    """Pruning only removes outliers, so the target can only go down."""
    rng = np.random.default_rng(2)
    for _ in range(20):
        n = 80
        a = rng.random((n, 3)) * 20
        b = a + rng.normal(0, 0.4, (n, 3))
        b[:6] += rng.normal(0, 12, (6, 3))  # a few wild outliers
        ia = list(range(n))
        raw = _kabsch(a, b)[0]
        pruned, kept = _rmsd_pruned(a, b)
        assert pruned <= raw + 1e-9, (pruned, raw)
        assert 30 <= kept <= n, kept


def test_alignment_blocks_are_monotonic():
    """The index pairs feed the coordinate gather, so gaps must be dropped and
    order preserved."""
    a = "ACDEFGHIKLMNPQRSTVWY"
    b = "ACDEFGHIKLMNPQRSTVWY"
    ia, ib = _alignment_blocks(ALIGNER.align(a, b)[0])
    assert ia == ib == list(range(20))

    ia, ib = _alignment_blocks(ALIGNER.align(a, "ACDEFGHIKLMNPQRSTVWY")[0])
    assert ia == sorted(ia) and ib == sorted(ib)
    assert len(ia) == len(ib)


def demo():
    """End-to-end on the repo's own five proteins, if the structures are present."""
    import os

    if not os.path.isdir("pdb_all"):
        print("pdb_all/ absent, run fetch_data.py first; skipping demo")
        return

    from pipeline import stage1_proteins, stage2_epitope

    proteins = stage1_proteins()
    by_pid = {p.pid: p for p in proteins}
    assert "1BV1_A" in by_pid, "1BV1 (Bet v 1) should be in the family set"

    # The repo's headline numbers must still reproduce from the new code path,
    # which is also the check that the NMR change did not move them.
    for pid_a, pid_b, expected in [("1BV1_A", "5MMU_A", 2.11),
                                   ("1BV1_A", "2BK0_A", 1.49),
                                   ("2BK0_A", "2WQL_A", 0.69)]:
        a, b = by_pid[pid_a], by_pid[pid_b]
        ia, ib = _alignment_blocks(ALIGNER.align(a.seq, b.seq)[0])
        got = _kabsch(a.ca[ia], b.ca[ib])[0]
        assert abs(got - expected) < 0.05, (pid_a, pid_b, got, expected)

    epitope, _ = stage2_epitope()
    assert len(epitope[4.5]) == 16, len(epitope[4.5])
    assert len(epitope[4.0]) == 15, len(epitope[4.0])
    print(f"{len(proteins)} proteins, epitope 4.5 A = {len(epitope[4.5])} residues, "
          "repo RMSD values reproduce")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
    demo()
    print("ok  demo")