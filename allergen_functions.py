from itertools import combinations

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser, Superimposer

aligner = PairwiseAligner(substitution_matrix = substitution_matrices.load("BLOSUM62"), open_gap_score = -11, extend_gap_score = -1, mode = "global")

def ids(structure, model_id=0):
    """Return the chain IDs in model"""
    ls = []
    for chain in structure[model_id]:
        ls.append(chain.id)
    return ls

def sep_aminos(residues):
    """Seperate amino acids"""
    amino_acid_residues = [residue for residue in residues if residue.get_id()[0] == ' ']
    return amino_acid_residues

def get_seq_and_ca(structure, chain_id, model_id=0):
    """Gives sequence of amino acids and list of CA atoms for a single chain"""
    residues = sep_aminos(structure[model_id][chain_id])
    letters = []
    list_of_CA_atoms = []

    for residue in residues:
        if "CA" in residue:
            letters.append(protein_letters_3to1[residue.get_resname().capitalize()])
            list_of_CA_atoms.append(residue["CA"])

    seq_string = "".join(letters)

    return seq_string, list_of_CA_atoms

def align_seq(seq1, seq2):
    """Returns the best global alignment of two sequences"""
    alignments = aligner.align(seq1, seq2)
    return alignments[0]

def pair_ca_atoms(algmt, caA, caB):
    """Walk the aligned blocks and return two equal-length lists of paired CA atoms."""
    pairedA = []
    pairedB = []

    rangeA = algmt.aligned[0]
    rangeB = algmt.aligned[1]
    zipped = zip(rangeA, rangeB)
    for ((startA, endA), (startB, endB)) in zipped:
        for i in range(endA - startA):
            a = startA + i
            b = startB + i
            pairedA.append(caA[a])
            pairedB.append(caB[b])
    return pairedA, pairedB

def pair_identicals(algmt, strucA, strucB, chainIdA, chainIdB):
    """Return (a, b) position pairs where the aligned amino acids are identical."""
    identical_pairs = []
    seqA, caA = get_seq_and_ca(strucA, chainIdA)
    seqB, caB = get_seq_and_ca(strucB, chainIdB)
    rangeA = algmt.aligned[0]
    rangeB = algmt.aligned[1]
    zipped = zip(rangeA, rangeB)
    for ((startA, endA), (startB, endB)) in zipped:
        for i in range(endA - startA):
            a = startA + i
            b = startB + i
            if seqA[a] == seqB[b]:
                identical_pairs.append((a, b))
    return identical_pairs

def percent_identity(algmt, strucA, strucB, chainIdA, chainIdB):
    """Identical aligned pairs as a % of the shorter sequence"""
    identical_pair_count = len(pair_identicals(algmt, strucA, strucB, chainIdA, chainIdB))
    seqA, caA = get_seq_and_ca(strucA, chainIdA)
    seqB, caB = get_seq_and_ca(strucB, chainIdB)
    shorterSeq = min(len(seqA), len(seqB))
    percent = identical_pair_count * 100 / shorterSeq
    return round(percent, 1)

def compute_rmsd(pairedA, pairedB):
    """Superimpose two paired CA atom lists and return the RMSD in Å."""
    rmsd_imposer = Superimposer()
    rmsd_imposer.set_atoms(pairedA, pairedB)
    return round(rmsd_imposer.rms, 2)

def allergen_comparison(proteinList):
    """Compare every unique pair of proteins; return one row per pair."""
    compareTable = []
    for p1, p2 in combinations(proteinList, 2):
        structure1 = p1["structure"]
        structure2 = p2["structure"]
        seq1, ca1 = get_seq_and_ca(structure1, p1["chainId"])
        seq2, ca2 = get_seq_and_ca(structure2, p2["chainId"])
        aligned = align_seq(seq1, seq2)
        paired1, paired2 = pair_ca_atoms(aligned, ca1, ca2)
        compareTable.append({
            "protein_a": p1["name"],
            "source_a": structure1.id,
            "protein_b": p2["name"],
            "source_b": structure2.id,
            "seq_identity_pct": percent_identity(aligned, structure1, structure2, p1["chainId"], p2["chainId"]),
            "rmsd_angstrom": compute_rmsd(paired1, paired2),
            "residues_aligned": len(paired1)
        })
    return pd.DataFrame(compareTable)