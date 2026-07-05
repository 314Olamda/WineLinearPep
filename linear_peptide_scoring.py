#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
linear_peptide_scoring.py — a/b/y + immonium scoring for linear peptides
══════════════════════════════════════════════════════════════════════════

Sibling module to WineCycloPep's score_spectrum_vs_cyclic, adapted for
LINEAR peptides (2–50 AA, optimized for the 2–10 AA range discussed
earlier). Built on mass_utils.py (Pyteomics-backed) rather than a
hand-rolled mass table, so it inherits the same fix already applied to
the cyclic pipeline.

Why the scoring logic differs from the cyclic version
────────────────────────────────────────────────────
Cyclic peptides (WineCycloPep):
  • n possible ring-opening positions → n overlapping bn-ion ladders
  • NO free C-terminus → y1 ABSENCE is diagnostic of cyclicity
  • Immonium ions: composition check only (low weight, 0.10)

Linear peptides (this module):
  • ONE fragmentation series → single b-ion ladder, single y-ion ladder
  • Free C-terminus → y-ion PRESENCE is diagnostic (opposite logic)
  • a-ions (b - CO) are a genuine third series worth scoring separately —
    for cyclic peptides a-ions are largely redundant with the multi-
    rotation b-ion ladder, but for linear peptides they're an independent
    confirmation, so they get their own scoring dimension here rather
    than being folded into "loss ions" as in the cyclic script.
  • Immonium ions play the same confirmatory (low-weight) role as in
    WineCycloPep — composition validation, not sequence-order evidence.

This module intentionally covers the SCORING function only — not a full
Bruker .d reading / de novo composition search pipeline. If you want the
full WineCycloPep-style pipeline (raw .d ingestion, composition search,
FDR/decoy validation, 3D structures) built out for linear peptides, that's
a separate, larger project — this gives you the core scoring engine to
drop into whatever spectrum source you're already using (e.g. your
MS-DIAL / GNPS exports).

Requirements
────────────
pip install pyteomics numpy
mass_utils.py must be in the same directory (or on PYTHONPATH)

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations

import numpy as np

from mass_utils import (
    RESIDUE_MASS, PROTON, H2O,
    linear_neutral_mass, precursor_neutral,
    immonium_ions,
)

# ══════════════════════════════════════════════════════════════════════════
# ION SERIES — a, b, y for a LINEAR peptide (single series, not rotated)
# ══════════════════════════════════════════════════════════════════════════

def b_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    b-ion series for a linear peptide: N-terminal fragments.
    b_i = sum(residues[0:i]) + proton, for i = 1 .. n-1
    """
    n = len(sequence)
    mz = []
    running = 0.0
    for i in range(n - 1):  # b_n (full length) is not a real fragment
        running += RESIDUE_MASS.get(sequence[i], 0.0)
        mz.append(round((running + PROTON * charge) / charge, 5))
    return mz


def a_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    a-ion series: b_i - CO. Independent confirmatory series for linear
    peptides (unlike the cyclic case, where a-ions are largely redundant
    with the multi-rotation b-ion ladder).
    """
    from pyteomics import mass as _mass
    co = _mass.calculate_mass(formula="CO")
    return [round(b - co / charge, 5) for b in b_ions(sequence, charge)]


def y_ions(sequence: str, charge: int = 1) -> list[float]:
    """
    y-ion series: C-terminal fragments, INCLUDING y1 (unlike the cyclic
    scorer, where y1 presence is treated as evidence AGAINST cyclicity).
    y_i = sum(residues[n-i:n]) + H2O + proton, for i = 1 .. n-1
    """
    n = len(sequence)
    mz = []
    running = 0.0
    for i in range(1, n):
        running += RESIDUE_MASS.get(sequence[n - i], 0.0)
        mz.append(round((running + H2O + PROTON * charge) / charge, 5))
    return mz


# ══════════════════════════════════════════════════════════════════════════
# SCORING
# ══════════════════════════════════════════════════════════════════════════

# Weights sum to 1.0. b/y carry the most sequence-order evidence; a-ions
# and immonium ions are confirmatory. Intensity rewards spectra where the
# matched ions are actually the dominant peaks, not just present above noise.
SCORE_WEIGHTS = {
    "b_coverage": 0.30,
    "y_coverage": 0.30,
    "a_coverage": 0.15,
    "intensity_score": 0.15,
    "immonium_score": 0.10,
}


def score_spectrum_vs_linear(
    mz_obs: np.ndarray,
    int_obs: np.ndarray,
    sequence: str,
    tol: float = 0.02,
    charge: int = 1,
) -> dict:
    """
    Score an MS2 spectrum against a linear peptide candidate sequence
    using a/b/y ion coverage + immonium confirmation.

    Mirrors score_spectrum_vs_cyclic's structure and return-dict shape
    (composite + per-dimension scores + matched/total counts) so it can
    be dropped into the same downstream FDR/decoy pipeline pattern used
    in WineCycloPep, if you extend this to a full de novo search later.
    """
    if len(mz_obs) == 0:
        return {"composite": 0.0}

    total_int = float(np.sum(int_obs))

    b_theo = b_ions(sequence, charge)
    a_theo = a_ions(sequence, charge)
    y_theo = y_ions(sequence, charge)
    imm_theo = list(immonium_ions(sequence).values())

    def _coverage(theo: list[float]) -> tuple[int, int, float]:
        if not theo:
            return 0, 0, 0.0
        matched = sum(1 for mz in theo if np.any(np.abs(mz_obs - mz) <= tol))
        return matched, len(theo), matched / len(theo)

    b_matched, b_total, b_coverage = _coverage(b_theo)
    a_matched, a_total, a_coverage = _coverage(a_theo)
    y_matched, y_total, y_coverage = _coverage(y_theo)
    imm_matched, imm_total, immonium_score = _coverage(imm_theo)

    # Intensity score: matched peaks' share of total ion current
    all_theo = b_theo + a_theo + y_theo
    matched_int = 0.0
    for mz in all_theo:
        hits = np.where(np.abs(mz_obs - mz) <= tol)[0]
        if len(hits) > 0:
            matched_int += float(np.max(int_obs[hits]))
    intensity_score = matched_int / total_int if total_int > 0 else 0.0

    composite = (
        SCORE_WEIGHTS["b_coverage"] * b_coverage +
        SCORE_WEIGHTS["y_coverage"] * y_coverage +
        SCORE_WEIGHTS["a_coverage"] * a_coverage +
        SCORE_WEIGHTS["intensity_score"] * intensity_score +
        SCORE_WEIGHTS["immonium_score"] * immonium_score
    )

    return {
        "composite": round(composite, 4),
        "b_coverage": round(b_coverage, 3),
        "y_coverage": round(y_coverage, 3),
        "a_coverage": round(a_coverage, 3),
        "intensity_score": round(intensity_score, 3),
        "immonium_score": round(immonium_score, 3),
        "b_matched": b_matched, "b_total": b_total,
        "y_matched": y_matched, "y_total": y_total,
        "a_matched": a_matched, "a_total": a_total,
        "immonium_matched": imm_matched, "immonium_total": imm_total,
    }


def confidence_tier(score: float) -> str:
    """Same tier boundaries as WineCycloPep, for consistency across the series."""
    if score >= 0.60:
        return "HIGH"
    if score >= 0.35:
        return "MEDIUM"
    if score >= 0.15:
        return "LOW"
    return "VERY_LOW"


# ══════════════════════════════════════════════════════════════════════════
# SELF-TEST
# ══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Synthetic test: perfect spectrum for linear tripeptide "IAA"
    seq = "IAA"
    b_theo = b_ions(seq)
    y_theo = y_ions(seq)
    imm_theo = list(immonium_ions(seq).values())

    print(f"Sequence: {seq} (linear)")
    print(f"Linear neutral mass: {linear_neutral_mass(seq):.5f} Da")
    print(f"b-ions: {b_theo}")
    print(f"a-ions: {a_ions(seq)}")
    print(f"y-ions: {y_theo}")
    print(f"Immonium ions: {imm_theo}")

    # Build a synthetic "perfect" spectrum from theoretical ions
    mz_obs = np.array(b_theo + y_theo + imm_theo)
    int_obs = np.ones_like(mz_obs) * 1000.0

    result = score_spectrum_vs_linear(mz_obs, int_obs, seq)
    print(f"\nScore: {result['composite']} → {confidence_tier(result['composite'])}")
    print(result)
