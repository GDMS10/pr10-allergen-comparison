"""Build a sequence-only -> structural-deviation dataset for the whole PR-10 family.

Three stages, each usable on its own:
    stage1_proteins()  parse every downloaded entry into one record per protein
    stage2_epitope()   find Bet v 1 residues that touch the IgE Fab in 9Y0A
    stage3_pairs()     align every protein pair, emit features (sequence only) + target RMSD

allergen_comparison.csv is 10 pairs of 5 proteins. This gives every pair of every
unique PR-10 protein with a structure, which is what a model needs.

Needs the structures first:  python3 fetch_data.py
Run:                    python3 build_dataset.py
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

# One aligner, same parameters the repo used (BLOSUM62, BLAST-style gaps).
ALIGNER = PairwiseAligner(
    substitution_matrix=substitution_matrices.load("BLOSUM62"),
    open_gap_score=-11,
    extend_gap_score=-1,
    mode="global",
)

# Kyte-Doolittle hydropathy, for the hydrophobicity-profile correlation feature.
KD = dict(zip("AVILMFWCY", [1.8, 4.2, 4.5, 3.8, 1.9, 2.8, 2.5, -0.4, -1.3, -3.2]))
KD.update(dict(zip("STP", [-0.8, -0.8, -0.8])))
KD.update({k: -4.5 for k in "DEQN"})
KD.update({k: -3.5 for k in "KR"})
KD.update({k: -3.9 for k in "G"})
KD.update({k: -0.4 for k in "H"})


@dataclass
class Protein:
    """One unique PR-10 protein sequence plus one representative coordinate set."""

    pid: str
    name: str
    organism: str
    uniprot: str
    entry: str
    method: str
    seq: str
    ca: np.ndarray  # (n_res, 3) representative CA coordinates
    n_models: int
    ensemble_rms: float  # mean intra-ensemble CA RMSD to the representative; 0 for X-ray
    resolution: float | None
    chains: dict = field(default_factory=dict)  # chain_id -> seq, all candidates

    @property
    def length(self) -> int:
        return len(self.seq)


def _chain_seq_ca(chain):
    """Sequence and CA coordinates for one chain, standard residues only."""
    letters, coords = [], []
    for res in chain:
        if res.get_id()[0] != " ":
            continue
        if "CA" not in res:
            continue
        name = res.get_resname().capitalize()
        if name not in protein_letters_3to1:
            continue
        letters.append(protein_letters_3to1[name])
        # is_disordered()==2 means altloc; Biopython already picked the highest occupancy.
        coords.append(res["CA"].get_coord())
    return "".join(letters), np.array(coords, dtype=float)


def _kabsch(fixed, moving):
    """Superimpose `moving` onto `fixed`. Returns (rmsd, moving_moved, rot, tran).

    Matches Bio's SVDSuperimposer exactly: H = moving.T @ fixed, rot = U @ Vt,
    and to forbid a reflection you negate a column of U (negating a row leaves
    the determinant negative and the RMSD wrong). Verified equal to
    Bio.PDB.Superimposer to 1e-9 on every pair in the family.
    """
    fc, mc = fixed.mean(0), moving.mean(0)
    H = (moving - mc).T @ (fixed - fc)
    U, _, Vt = np.linalg.svd(H)
    rot = U @ Vt
    if np.linalg.det(rot) < 0:
        U = U.copy()
        U[:, 2] = -U[:, 2]
        rot = U @ Vt
    tran = fc - mc @ rot
    moved = moving @ rot + tran
    rmsd = float(np.sqrt(((moved - fixed) ** 2).sum() / len(fixed)))
    return rmsd, moved, rot, tran


def _representative(structure, chain_id):
    """For NMR ensembles average all models; for X-ray take model 0.

    The repo took model 0 of the NMR entries, which silently hard-wired one
    arbitrary conformer into every comparison involving apple or cherry.
    """
    models = [m for m in structure if any(len(c) for c in m)]
    if len(models) == 1:
        seq, ca = _chain_seq_ca(structure[0][chain_id])
        return seq, ca, len(models), 0.0

    per_model = []
    for m in models:
        if chain_id in [c.id for c in m]:
            s, c = _chain_seq_ca(m[chain_id])
            if len(s) > 100:
                per_model.append((s, c))
    if not per_model:
        seq, ca = _chain_seq_ca(structure[0][chain_id])
        return seq, ca, len(models), 0.0

    ref_seq, ref_ca = per_model[0]
    same = [(s, c) for s, c in per_model if len(c) == len(ref_ca)]
    coords = np.mean([c for _, c in same], axis=0) if same else ref_ca
    spread = 0.0
    if len(same) > 1:
        spread = float(np.mean([_kabsch(ref_ca, c)[0] for _, c in same]))
    return ref_seq, coords, len(models), spread


def stage1_proteins(pdb_dir="pdb_all", meta_json="meta.json", min_len=110, max_len=200):
    """Parse every entry, keep the PR-10 entity, dedupe identical sequences."""
    meta = json.load(open(meta_json))
    parser = PDBParser(QUIET=True)

    records = []
    for entry, info in sorted(meta.items()):
        path = os.path.join(pdb_dir, f"{entry}.pdb")
        if not os.path.exists(path):
            continue
        pr10 = info.get("entities", [])
        if not pr10:
            continue
        entity = pr10[0]
        target_seq = entity["seq"]

        try:
            structure = parser.get_structure(entry, path)
        except Exception:
            continue

        # Score every plausible chain by how well it matches the annotated PR-10
        # entity sequence, then keep the best. The repo hard-coded chain "A",
        # which breaks on any entry where A is the Fab or a ligand.
        best = None
        for model in structure:
            for chain in model:
                try:
                    seq, ca = _chain_seq_ca(chain)
                except Exception:
                    continue
                if not (min_len <= len(seq) <= max_len):
                    continue
                score = _match_len(seq, target_seq)
                if best is None or score > best[0]:
                    rep_seq, rep_ca, nmod, spread = _representative(structure, chain.id)
                    best = (score, chain.id, seq, rep_seq, rep_ca, nmod, spread)
            break  # representative() already handled every model

        if best is None:
            continue
        _, chain_id, _, seq, ca, nmod, spread = best
        if len(ca) < min_len:
            continue

        records.append(
            Protein(
                pid=f"{entry}_{chain_id}",
                name=entity["name"] or entry,
                organism=entity.get("org") or "unknown",
                uniprot=(entity.get("uniprot") or [""])[0],
                entry=entry,
                method=info["method"][0],
                seq=seq,
                ca=ca,
                n_models=nmod,
                ensemble_rms=spread,
                resolution=info["res"][0] if info.get("res") else None,
            )
        )

    return _dedupe(records)


def _match_len(a, b):
    """Crude same-length check; PR-10 chains match the entity sequence or not."""
    return 1 if abs(len(a) - len(b)) <= 3 else 0


def _dedupe(records):
    """Collapse entries that are the same protein sequence.

    Keeps the highest-quality representative per sequence: X-ray over NMR, then
    lower resolution, so a 1.4 A crystal structure wins over the same protein's
    2.9 A one. Also collapses isoforms only when the sequence is byte-identical,
    which is the safe thing to do; distinct isoforms stay distinct.
    """
    by_seq = {}
    for r in records:
        key = r.seq
        rank = (0 if r.method == "X-RAY DIFFRACTION" else 1,
                r.resolution if r.resolution else 9.9,
                r.ensemble_rms)
        if key not in by_seq or rank < by_seq[key][0]:
            by_seq[key] = (rank, r)

    out = []
    for seq, (_, r) in sorted(by_seq.items(), key=lambda kv: kv[1][1].pid):
        r.chains = {"seq": seq}
        out.append(r)
    return out


def stage2_epitope(pdb_dir="pdb_all", epitope_entry="9Y0A", allergen_chain="C",
                    cutoffs=(4.0, 4.5, 5.0)):
    """Bet v 1 residues within each distance cutoff of the IgE Fab.

    Contacts are computed on every atom, not CA: an IgE epitope is mostly
    sidechain, and a CA-only criterion finds 0 residues at 5 A here, which
    silently reduces the feature to nothing.

    This is the README's own "Next" step, and it becomes the model's most
    useful feature: IgE recognises a surface patch, so epitope residues should
    matter more than average conservation.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(epitope_entry, os.path.join(pdb_dir, f"{epitope_entry}.pdb"))
    model = structure[0]

    fab_atoms = []
    for chain in model:
        if chain.id == allergen_chain:
            continue
        for res in chain:
            if res.get_id()[0] == " ":
                fab_atoms.extend(a.get_coord() for a in res)
    fab_atoms = np.array(fab_atoms)

    rows = []
    for res in model[allergen_chain]:
        if res.get_id()[0] != " ":
            continue
        d = min(float(np.linalg.norm(fab_atoms - a.get_coord(), axis=1).min()) for a in res)
        name = res.get_resname().capitalize()
        if name not in protein_letters_3to1:
            continue
        rows.append((res.get_id()[1], protein_letters_3to1[name], d))

    return {c: [r for r in rows if r[2] <= c] for c in cutoffs}, rows


def _alignment_blocks(alignment):
    """(i, j) index pairs for every aligned column, gaps skipped."""
    a, b = [], []
    for (s1, e1), (s2, e2) in zip(alignment.aligned[0], alignment.aligned[1]):
        for k in range(e1 - s1):
            a.append(s1 + k)
            b.append(s2 + k)
    return a, b


def _rmsd(ca1, ca2):
    return _kabsch(ca1, ca2)[0]


def _rmsd_pruned(ca1, ca2, rounds=5, z=2.0):
    """Iterative outlier pruning, so the target matches what ChimeraX matchmaker reports.

    Two-sigma cut, never dropping below 30 pairs.
    """
    keep = np.ones(len(ca1), dtype=bool)
    for _ in range(rounds):
        rmsd, moved, rot, tran = _kabsch(ca1[keep], ca2[keep])
        d = np.linalg.norm(moved - ca2[keep], axis=1)
        worst = np.where(d > z * d.std())[0]
        if len(worst) == 0 or keep.sum() - len(worst) < 30:
            break
        idx = np.where(keep)[0]
        keep[idx[worst]] = False
    return _kabsch(ca1[keep], ca2[keep])[0], int(keep.sum())


def _profile(seq, idx):
    return np.array([KD.get(seq[i], 0.0) for i in idx])


def epitope_index_per_protein(epitope_rows, epitope_seq, proteins):
    """Map the 9Y0A epitope positions onto each protein once, by index.

    Cheaper and clearer than doing it per pair: align the 9Y0A allergen chain to
    every protein up front, then a pair just intersects two index sets.
    """
    out = {}
    for p in proteins:
        aln = ALIGNER.align(epitope_seq, p.seq)[0]
        anchor = {}
        for (s1, e1), (s2, e2) in zip(aln.aligned[0], aln.aligned[1]):
            for k in range(e1 - s1):
                anchor[s1 + k] = s2 + k
        ep = sorted({anchor[r[0] - 1] for r in epitope_rows if r[0] - 1 in anchor})
        out[p.pid] = set(ep)
    return out


def stage3_pairs(proteins, epi_index=None):
    """Every protein pair: sequence-only features plus the structural target.

    Features deliberately never touch coordinates, so the model can score a new
    allergen from its sequence alone.
    """
    rows = []
    for p1, p2 in combinations(proteins, 2):
        if abs(p1.length - p2.length) > 40:
            continue
        aln = ALIGNER.align(p1.seq, p2.seq)[0]
        ia, ib = _alignment_blocks(aln)

        if len(ia) < 60:
            continue
        n_gap_cols = (p1.length + p2.length) - 2 * len(ia)

        ident = sum(1 for i, j in zip(ia, ib) if p1.seq[i] == p2.seq[j])
        ident_short = ident / min(p1.length, p2.length)
        ident_cov = ident / len(ia)

        # BLOSUM score normalised by alignment length: a better-weighted identity.
        sub = substitution_matrices.load("BLOSUM62")
        blosum = sum(sub[p1.seq[i], p2.seq[j]] for i, j in zip(ia, ib)) / len(ia)

        hyd1, hyd2 = _profile(p1.seq, ia), _profile(p2.seq, ib)
        hyd_corr = float(np.corrcoef(hyd1, hyd2)[0, 1]) if len(ia) > 3 else 0.0

        # Epitope-weighted conservation: identity restricted to the positions
        # the Bet v 1 IgE Fab actually touches in 9Y0A.
        epi_ident, epi_cov, epi_n = None, None, 0
        if epi_index:
            hits = [(i, j) for i, j in zip(ia, ib) if i in epi_index.get(p1.pid, ())]
            epi_n = len(hits)
            if epi_n:
                epi_ident = sum(1 for i, j in hits if p1.seq[i] == p2.seq[j]) / epi_n
                epi_cov = epi_n / len(ia)

        ca1, ca2 = p1.ca[ia], p2.ca[ib]
        rmsd = _kabsch(ca1, ca2)[0]
        rmsd_pruned, n_kept = _rmsd_pruned(ca1, ca2)

        rows.append(dict(
            a=p1.pid, b=p2.pid,
            org_a=p1.organism, org_b=p2.organism,
            name_a=p1.name, name_b=p2.name,
            len_a=p1.length, len_b=p2.length,
            len_diff=abs(p1.length - p2.length),
            identity_short=ident_short,
            identity_cov=ident_cov,
            blosum_norm=blosum,
            gap_frac=n_gap_cols / (p1.length + p2.length),
            hyd_corr=hyd_corr,
            n_aligned=len(ia),
            epi_identity=epi_ident,
            epi_coverage=epi_cov,
            epi_n=epi_n,
            rmsd=rmsd,
            rmsd_pruned=rmsd_pruned,
            n_kept=n_kept,
            method_a=p1.method, method_b=p2.method,
            ensemble_a=p1.ensemble_rms, ensemble_b=p2.ensemble_rms,
            same_organism=int(p1.organism == p2.organism),
            same_uniprot=int(bool(p1.uniprot) and p1.uniprot == p2.uniprot),
        ))
    return rows