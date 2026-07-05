#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
linear_denovo.py — De novo linear peptide identification (2–50 AA)
════════════════════════════════════════════════════════════════════

Full pipeline counterpart to WineCycloPep's Step 3, adapted for LINEAR
peptides: mass-directed composition search -> permutation candidates ->
a/b/y + immonium scoring (linear_peptide_scoring.py) -> target-decoy FDR
validation (Benjamini-Hochberg + Storey q-value + Hotelling T²/F).

What's ported unchanged from WineCycloPep's statistical framework
────────────────────────────────────────────────────────────────
The FDR machinery (empirical p-value, BH correction, Storey q-value,
Hotelling T²) is architecture-agnostic — it operates on score vectors,
not on cyclic-vs-linear fragment logic. Only three things change:
  1. Neutral mass uses +H2O (linear_neutral_mass), not the cyclic
     no-water rule
  2. Decoys exclude only the exact candidate string (no rotational
     equivalents to exclude, since linear sequences aren't cyclic)
  3. Score dimensions are b/y/a/immonium (linear_peptide_scoring.py),
     not bn/loss/y1-absence/immonium (WineCycloPep)

Requirements
────────────
pip install pyteomics numpy scipy
mass_utils.py and linear_peptide_scoring.py must be in the same directory.

Author : Pol Giménez-Gil — ISVV, Université de Bordeaux
ORCID  : 0000-0002-7720-3733
"""

from __future__ import annotations

import itertools
import random
from dataclasses import dataclass, field

import numpy as np
from scipy import stats as _stats

from mass_utils import RESIDUE_MASS, linear_neutral_mass, precursor_neutral
from linear_peptide_scoring import (
    score_spectrum_vs_linear, confidence_tier, SCORE_WEIGHTS,
)

SCORE_DIMS = ["b_coverage", "y_coverage", "a_coverage",
              "intensity_score", "immonium_score"]
MIN_DECOYS_FOR_PVAL = 20


# ══════════════════════════════════════════════════════════════════════════
# 1. MASS-DIRECTED COMPOSITION SEARCH  (ported from WineCycloPep,
#    target mass swapped from cyclic to linear)
# ══════════════════════════════════════════════════════════════════════════

def _all_compositions(
    target_mass: float,
    residues: list[tuple[str, float]],
    n_min: int,
    n_max: int,
    tol: float,
    current: list[str],
    current_mass: float,
    results: list[list[str]],
) -> None:
    """Recursive mass-directed composition finder (branch-and-bound).
    Identical algorithm to WineCycloPep — the branch-and-bound logic
    doesn't depend on cyclic vs. linear, only the target_mass passed in
    does (see compositions_for_linear_mass, which subtracts H2O upfront)."""
    n = len(current)
    if n >= n_min and abs(current_mass - target_mass) <= tol:
        results.append(current[:])
    if n >= n_max or current_mass > target_mass + tol:
        return
    last = current[-1] if current else None
    for aa, mass_val in residues:
        if last and aa < last:
            continue
        new_mass = current_mass + mass_val
        if new_mass > target_mass + tol:
            continue
        current.append(aa)
        _all_compositions(target_mass, residues, n_min, n_max, tol,
                           current, new_mass, results)
        current.pop()


def compositions_for_linear_mass(
    precursor_neutral_mass: float,
    n_min: int = 2,
    n_max: int = 10,
    tol: float = 0.02,
    allowed_aa: str | None = None,
) -> list[str]:
    """
    Find all amino acid compositions whose linear neutral mass
    (sum of residues + H2O) matches precursor_neutral_mass within tol.
    Subtracts H2O once here so the branch-and-bound target is a pure
    residue-mass sum, same as the cyclic search.
    """
    from mass_utils import H2O
    target = precursor_neutral_mass - H2O
    if allowed_aa:
        residues = [(aa, RESIDUE_MASS[aa]) for aa in allowed_aa if aa in RESIDUE_MASS]
    else:
        residues = sorted(RESIDUE_MASS.items(), key=lambda x: x[0])
    results: list[list[str]] = []
    _all_compositions(target, residues, n_min, n_max, tol, [], 0.0, results)
    return ["".join(r) for r in results]


def all_permutations_linear(composition: str, max_perms: int | None = None) -> list[str]:
    """
    All unique orderings of a composition. Unlike the cyclic case, there
    is no rotational equivalence to collapse — cyclo(IAA), cyclo(AAI),
    cyclo(AIA) are the same molecule, but linear IAA, AAI, AIA are three
    genuinely distinct peptides with different b/y ladders. So this
    returns every permutation, deduplicated only for literal repeats
    (relevant when the composition has repeated residues, e.g. "AAI").
    """
    seen: set[str] = set()
    out: list[str] = []
    for perm in itertools.permutations(composition):
        s = "".join(perm)
        if s not in seen:
            seen.add(s)
            out.append(s)
    if max_perms and len(out) > max_perms:
        out = random.sample(out, max_perms)
    return out


# ══════════════════════════════════════════════════════════════════════════
# 2. DECOY GENERATION  (simpler than cyclic: only exclude the exact
#    candidate, no rotational family to exclude)
# ══════════════════════════════════════════════════════════════════════════

def generate_decoys_linear(sequence: str, n_decoys: int = 100) -> list[str]:
    """
    Decoys for a linear candidate: permutations of the same composition,
    excluding the candidate itself. Same isobaric I<->L fallback as
    WineCycloPep for short/degenerate compositions.
    """
    n = len(sequence)
    all_perms = {"".join(p) for p in itertools.permutations(sequence)}
    valid = all_perms - {sequence}

    if len(valid) < MIN_DECOYS_FOR_PVAL:
        pseudo: set[str] = set()
        for i, aa in enumerate(sequence):
            sub = "L" if aa == "I" else ("I" if aa == "L" else None)
            if sub:
                s = sequence[:i] + sub + sequence[i + 1:]
                if s != sequence:
                    pseudo.add(s)
                    for p in itertools.permutations(s):
                        ps = "".join(p)
                        if ps != sequence:
                            pseudo.add(ps)
        valid |= pseudo

    result = list(valid)
    if len(result) > n_decoys:
        result = random.sample(result, n_decoys)
    return result


# ══════════════════════════════════════════════════════════════════════════
# 3. FDR STATISTICS  (identical algorithms to WineCycloPep — reusable
#    verbatim since they operate on generic score vectors)
# ══════════════════════════════════════════════════════════════════════════

def compute_empirical_pvalue(target_score: float, decoy_scores: np.ndarray) -> float:
    """One-tailed empirical p-value with Laplace pseudocount."""
    if len(decoy_scores) == 0:
        return np.nan
    return float((np.sum(decoy_scores >= target_score) + 1) / (len(decoy_scores) + 1))


def bh_fdr(p_values: np.ndarray, alpha: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Benjamini-Hochberg (1995) step-up FDR correction."""
    n = len(p_values)
    if n == 0:
        return np.array([]), np.array([], dtype=bool)
    order = np.argsort(p_values)
    ranked = np.arange(1, n + 1, dtype=float)
    bh_thresh = ranked / n * alpha
    below = p_values[order] <= bh_thresh
    k_max = int(np.max(np.where(below)[0])) if np.any(below) else -1
    rejected = np.zeros(n, dtype=bool)
    if k_max >= 0:
        rejected[order[:k_max + 1]] = True
    adj = np.minimum(1.0, p_values[order] * n / ranked)
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    adj_out = np.zeros(n)
    adj_out[order] = adj
    return adj_out, rejected


def storey_qvalue(p_values: np.ndarray, lambdas: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Storey & Tibshirani (2003) q-value with pi0 estimation."""
    n = len(p_values)
    if n == 0:
        return np.array([]), 1.0
    if lambdas is None:
        lambdas = np.arange(0.05, 0.90, 0.05)
    pi0_hats = np.array([np.sum(p_values > lam) / (n * (1.0 - lam)) for lam in lambdas])
    pi0 = float(np.clip(np.min(pi0_hats[-5:]) if len(pi0_hats) >= 5 else pi0_hats[-1], 0.0, 1.0))
    order = np.argsort(p_values)
    ranked = np.arange(1, n + 1, dtype=float)
    q = p_values[order] * n * pi0 / ranked
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    q_out = np.zeros(n)
    q_out[order] = np.minimum(1.0, q)
    return q_out, pi0


def hotellings_t2(X_target: np.ndarray, X_decoy: np.ndarray) -> dict:
    """
    Multivariate Hotelling's T² / F-test: do the 5 score dimensions
    JOINTLY separate targets from the global decoy pool?

    T² = (n1*n2)/(n1+n2) * (mu1-mu2)^T * Sp^-1 * (mu1-mu2)
    F  = T² * (n1+n2-p-1) / ((n1+n2-2) * p)  ~  F(p, n1+n2-p-1)

    Also reports per-dimension Mann-Whitney U (one-tailed, target > decoy)
    with BH correction, identifying which dimensions drive the separation.
    """
    n1, n2 = len(X_target), len(X_decoy)
    p = X_target.shape[1]
    if n1 < 2 or n2 < 2:
        return {"T2": np.nan, "F": np.nan, "p_value": np.nan, "df1": p, "df2": np.nan}

    mu1, mu2 = X_target.mean(axis=0), X_decoy.mean(axis=0)
    S1, S2 = np.cov(X_target, rowvar=False), np.cov(X_decoy, rowvar=False)
    Sp = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)

    try:
        Sp_inv = np.linalg.pinv(Sp)
    except np.linalg.LinAlgError:
        return {"T2": np.nan, "F": np.nan, "p_value": np.nan, "df1": p, "df2": np.nan}

    diff = (mu1 - mu2).reshape(-1, 1)
    T2 = float((n1 * n2) / (n1 + n2) * (diff.T @ Sp_inv @ diff)[0, 0])
    df1, df2 = p, n1 + n2 - p - 1
    F = T2 * df2 / ((n1 + n2 - 2) * df1) if df2 > 0 else np.nan
    p_val = float(1 - _stats.f.cdf(F, df1, df2)) if df2 > 0 and not np.isnan(F) else np.nan

    per_dim = {}
    for i, dim in enumerate(SCORE_DIMS):
        try:
            u_stat, u_p = _stats.mannwhitneyu(
                X_target[:, i], X_decoy[:, i], alternative="greater"
            )
        except ValueError:
            u_stat, u_p = np.nan, np.nan
        per_dim[dim] = {"U": u_stat, "p_value": u_p}

    dim_p = np.array([per_dim[d]["p_value"] for d in SCORE_DIMS])
    valid = ~np.isnan(dim_p)
    if valid.any():
        adj = np.full_like(dim_p, np.nan)
        adj[valid], _ = bh_fdr(dim_p[valid])
        for i, dim in enumerate(SCORE_DIMS):
            per_dim[dim]["p_adj_BH"] = adj[i]

    return {"T2": T2, "F": F, "p_value": p_val, "df1": df1, "df2": df2,
            "per_dimension": per_dim}


# ══════════════════════════════════════════════════════════════════════════
# 4. MAIN DRIVER — candidate identification for one MS2 spectrum
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class LinearCandidate:
    sequence: str
    composite: float
    scores: dict
    p_value: float = np.nan
    p_adj_BH: float = np.nan
    q_storey: float = np.nan
    rejected_BH: bool = False
    n_decoys_in_pool: int = 0
    stat_note: str = "ok"


def identify_linear_peptide(
    mz_obs: np.ndarray,
    int_obs: np.ndarray,
    precursor_mz: float,
    charge: int,
    n_min: int = 2,
    n_max: int = 10,
    tol: float = 0.02,
    min_score: float = 0.10,
    max_perms_per_composition: int = 200,
    n_decoys: int = 100,
    fdr_alpha: float = 0.05,
    allowed_aa: str | None = None,
) -> list[LinearCandidate]:
    """
    Full de novo identification for one spectrum:
      1. Compute precursor neutral mass, search compositions
      2. Enumerate permutations per composition (capped)
      3. Score every candidate against the spectrum
      4. Build a global decoy pool (all candidates' decoys pooled together,
         matching WineCycloPep's approach) and compute empirical p-values
      5. BH + Storey correction; Hotelling T²/F on the full score matrix

    Returns candidates sorted by q_storey ascending (best first).
    """
    neutral_mass = precursor_neutral(precursor_mz, charge)
    compositions = compositions_for_linear_mass(
        neutral_mass, n_min=n_min, n_max=n_max, tol=tol, allowed_aa=allowed_aa
    )

    candidates: list[dict] = []
    all_decoy_scores: list[float] = []
    candidate_decoy_map: dict[str, np.ndarray] = {}

    for comp in compositions:
        for seq in all_permutations_linear(comp, max_perms=max_perms_per_composition):
            result = score_spectrum_vs_linear(mz_obs, int_obs, seq, tol=tol, charge=1)
            if result.get("composite", 0.0) < min_score:
                continue
            candidates.append({"sequence": seq, **result})

            decoys = generate_decoys_linear(seq, n_decoys=n_decoys)
            decoy_scores = np.array([
                score_spectrum_vs_linear(mz_obs, int_obs, d, tol=tol, charge=1).get("composite", 0.0)
                for d in decoys
            ])
            candidate_decoy_map[seq] = decoy_scores
            all_decoy_scores.extend(decoy_scores.tolist())

    if not candidates:
        return []

    global_decoy_pool = np.array(all_decoy_scores)

    out: list[LinearCandidate] = []
    for cand in candidates:
        seq = cand["sequence"]
        local_decoys = candidate_decoy_map.get(seq, np.array([]))
        pool = global_decoy_pool if len(global_decoy_pool) >= MIN_DECOYS_FOR_PVAL else local_decoys
        note = "ok" if len(local_decoys) >= MIN_DECOYS_FOR_PVAL else "low_decoy_count"

        p_val = compute_empirical_pvalue(cand["composite"], pool)
        out.append(LinearCandidate(
            sequence=seq, composite=cand["composite"], scores=cand,
            p_value=p_val, n_decoys_in_pool=len(pool), stat_note=note,
        ))

    p_arr = np.array([c.p_value for c in out])
    adj_p, rejected = bh_fdr(p_arr, alpha=fdr_alpha)
    q_arr, pi0 = storey_qvalue(p_arr)
    for i, c in enumerate(out):
        c.p_adj_BH = adj_p[i]
        c.rejected_BH = bool(rejected[i])
        c.q_storey = q_arr[i]

    out.sort(key=lambda c: (c.q_storey, -c.composite))
    return out


if __name__ == "__main__":
    # Smoke test with a synthetic perfect spectrum for linear "IAA"
    from linear_peptide_scoring import b_ions, y_ions
    from mass_utils import immonium_ions, linear_neutral_mass, PROTON

    seq = "IAA"
    b_theo, y_theo = b_ions(seq), y_ions(seq)
    imm_theo = list(immonium_ions(seq).values())
    mz_obs = np.array(b_theo + y_theo + imm_theo)
    int_obs = np.ones_like(mz_obs) * 1000.0

    precursor_mz = linear_neutral_mass(seq) + PROTON  # [M+H]+, charge 1

    print(f"Testing de novo identification against synthetic spectrum for {seq}")
    print(f"Precursor m/z: {precursor_mz:.5f} (charge 1)\n")

    results = identify_linear_peptide(
        mz_obs, int_obs, precursor_mz, charge=1,
        n_min=2, n_max=4, tol=0.02, min_score=0.05,
        n_decoys=50,
    )

    print(f"{len(results)} candidate(s) above threshold:\n")
    for c in results[:5]:
        print(f"  {c.sequence:6s}  composite={c.composite:.3f}  "
              f"tier={confidence_tier(c.composite):9s}  "
              f"p={c.p_value:.4f}  q_storey={c.q_storey:.4f}  "
              f"rejected_BH={c.rejected_BH}")
