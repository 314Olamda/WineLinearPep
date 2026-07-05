#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mass_utils.py — Pyteomics-backed mass/ion layer for WineCycloPep
══════════════════════════════════════════════════════════════════

Drop-in replacement for the hand-rolled RESIDUE_MASS dict, bn-ion
generator, residue-loss calculator, and immonium-ion formula in
wine_cyclopep.py Step 3 (cyclic peptide de novo detection).

Why this exists
────────────────
The original functions in wine_cyclopep.py (bn_ions_all_rotations,
residue_loss_ions, and the immonium-ion inline formula in
score_spectrum_vs_cyclic) all recompute standard peptide mass
chemistry from a manually transcribed monoisotopic mass table. That
table is correct today, but every future edit — adding a PTM,
extending to isotope-labelled residues, correcting a mass to more
decimal places — is a manual, unverified change feeding directly
into your FDR/Storey-q statistics.

Pyteomics ships this table (and PTM/isotope handling) as a maintained,
peer-reviewed reference implementation. This module is a thin
compatibility layer: same function names and return types as the
original, so integration into wine_cyclopep.py is a one-line import
swap, not a rewrite of the scoring logic.

Requirements
────────────
pip install pyteomics

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations

from pyteomics import mass

# ══════════════════════════════════════════════════════════════════════════
# CONSTANTS — now sourced from Pyteomics instead of hand-transcribed
# ══════════════════════════════════════════════════════════════════════════

PROTON = mass.calculate_mass(formula="H") - mass.nist_mass["e*"][0][0]
# Equivalent to the original PROTON = 1.007276; kept as a fallback below
# for environments where the electron-mass lookup differs by version.
if not (1.007 < PROTON < 1.008):
    PROTON = 1.007276

H2O = mass.calculate_mass(formula="H2O")

# Monoisotopic residue masses — pulled directly from Pyteomics' std_aa_mass,
# which is validated against the same IUPAC values as the original dict.
# Proline (P) is included here (unlike AA_SIDECHAIN in the 3D builder,
# which correctly excludes it — that limitation is about SMILES ring
# construction, not mass, and is unaffected by this change).
#
# IMPORTANT: Pyteomics' std_aa_mass also includes non-canonical single-letter
# codes — J (Leu/Ile ambiguity, identical mass to I/L), U (selenocysteine),
# O (pyrrolysine). Left unfiltered, these silently expand a de novo
# composition search with a duplicate-mass ambiguity code (J) and two
# residues (U, O) that don't occur in yeast/grapevine peptidomes relevant
# to wine lees work. Restricting explicitly to the 20 canonical residues
# is what the original hardcoded dict did implicitly — this preserves
# that behaviour rather than silently changing it.
_CANONICAL_20 = set("ACDEFGHIKLMNPQRSTVWY")
RESIDUE_MASS: dict[str, float] = {
    aa: m for aa, m in mass.std_aa_mass.items()
    if len(aa) == 1 and aa in _CANONICAL_20
}

# CO mass, used for immonium ion calculation (residue - CO + H)
_CO_MASS = mass.calculate_mass(formula="CO")


def cyclic_neutral_mass(sequence: str) -> float:
    """Neutral mass of a head-to-tail cyclic peptide (no H2O). Unchanged API."""
    return sum(RESIDUE_MASS.get(aa, 0.0) for aa in sequence)


def linear_neutral_mass(sequence: str) -> float:
    """Neutral mass of a linear peptide (+H2O). Unchanged API."""
    return sum(RESIDUE_MASS.get(aa, 0.0) for aa in sequence) + H2O


def precursor_neutral(mz: float, charge: int) -> float:
    return mz * charge - charge * PROTON


# ══════════════════════════════════════════════════════════════════════════
# ION CALCULATIONS — replaces bn_ions_all_rotations / residue_loss_ions /
# the inline immonium formula in score_spectrum_vs_cyclic
# ══════════════════════════════════════════════════════════════════════════

def bn_ions_all_rotations(sequence: str, charge: int = 1) -> list[float]:
    """
    All bn ions across every ring-opening rotation of a cyclic peptide.
    Same signature and return type as the original — safe drop-in.
    """
    n = len(sequence)
    mz_set: set[float] = set()
    for start in range(n):
        rot = sequence[start:] + sequence[:start]
        b_mass = sum(RESIDUE_MASS.get(aa, 0.0) for aa in rot[:-1])
        mz = round((b_mass + PROTON * charge) / charge, 5)
        mz_set.add(mz)
    return sorted(mz_set)


def residue_loss_ions(sequence: str, cyclic_mh: float) -> dict[str, float]:
    """[M+H]+ − residue_mass for each unique residue. Same API as original."""
    return {
        f"loss_{aa}": round(cyclic_mh - RESIDUE_MASS[aa], 5)
        for aa in set(sequence)
    }


def immonium_ions(sequence: str) -> dict[str, float]:
    """
    Immonium ion m/z per unique residue: residue_mass - CO + H+.
    Previously inlined in score_spectrum_vs_cyclic using a hardcoded
    27.9949 for CO — now sourced from Pyteomics' formula mass calculator.
    """
    return {
        aa: round(RESIDUE_MASS[aa] - _CO_MASS + PROTON, 5)
        for aa in set(sequence)
    }


def y1_ion(aa: str) -> float:
    """y1 ion for a single C-terminal residue: residue_mass + H2O + H+."""
    return RESIDUE_MASS[aa] + H2O + PROTON


# ══════════════════════════════════════════════════════════════════════════
# OPTIONAL: MSCI-based decoy indistinguishability check
# ══════════════════════════════════════════════════════════════════════════
#
# generate_decoys() in wine_cyclopep.py already handles the short-sequence
# degeneracy problem with an ad hoc I<->L isobaric substitution. MSCI
# formalizes this: instead of only patching I/L, it scores fragmentation
# *similarity* between any two same-composition candidates and flags pairs
# that are genuinely indistinguishable by MS2 alone — useful for sequences
# like AAI vs AIA vs IAA where composition is shared but rotation differs.
#
# This is presented separately (not wired into generate_decoys) because it
# changes what "valid decoy" means — worth validating against a known
# spectrum before it feeds your FDR pipeline.
#
# pip install msci   (see MSCI docs for current package name/import path)
#
# def flag_indistinguishable(sequence: str, decoy_candidates: list[str]) -> list[str]:
#     """Return the subset of decoy_candidates whose theoretical fragmentation
#     is statistically indistinguishable from `sequence` per MSCI scoring."""
#     from msci import compare_fragmentation  # exact import per MSCI version
#     flagged = []
#     for cand in decoy_candidates:
#         sim = compare_fragmentation(sequence, cand, ion_types=("b", "y"))
#         if sim.is_indistinguishable:
#             flagged.append(cand)
#     return flagged


if __name__ == "__main__":
    # Quick self-check against the original hardcoded values
    test_seq = "IAA"
    print(f"Sequence: {test_seq}")
    print(f"Cyclic neutral mass: {cyclic_neutral_mass(test_seq):.5f} Da")
    print(f"bn ions (all rotations): {bn_ions_all_rotations(test_seq)}")
    print(f"Residue loss ions: {residue_loss_ions(test_seq, cyclic_neutral_mass(test_seq) + PROTON)}")
    print(f"Immonium ions: {immonium_ions(test_seq)}")
